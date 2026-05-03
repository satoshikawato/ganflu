#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import csv
import json
import os
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from ganflu.launchers.miniprot import MiniprotCommandLine
from ganflu.scripts import gff3_prune, gff3togbk, validate_reference_files


DEFAULT_AUTO_TARGETS = ("IAV", "IBV", "ICV", "IDV")
MINIPROT_PREFIXES = {
    "IAV": "MPIA",
    "IBV": "MPIB",
    "ICV": "MPIC",
    "IDV": "MPID",
}
TSV_COLUMNS = [
    "contig_id",
    "length",
    "call",
    "target",
    "segment",
    "status",
    "flags",
    "confidence",
    "best_identity",
    "best_aa_coverage",
    "best_score",
    "second_target",
    "second_score",
    "score_margin",
    "best_product",
    "best_target_range",
    "best_query_range",
    "strand",
    "cds_count",
    "frameshift_count",
    "internal_stop_count",
    "missing_start_count",
    "missing_stop_count",
    "notes",
]


@dataclass(frozen=True)
class AutoThresholds:
    min_identity: float = 0.55
    min_aa_coverage: float = 0.35
    min_score: float = 0.25
    min_margin: float = 0.10
    complete_aa_coverage: float = 0.90

    @classmethod
    def from_args(cls, args):
        return cls(
            min_identity=float(args.auto_min_identity),
            min_aa_coverage=float(args.auto_min_aa_coverage),
            min_score=float(args.auto_min_score),
            min_margin=float(args.auto_min_margin),
            complete_aa_coverage=float(args.auto_complete_aa_coverage),
        )

    def as_dict(self):
        return {
            "min_identity": self.min_identity,
            "min_aa_coverage": self.min_aa_coverage,
            "min_score": self.min_score,
            "min_margin": self.min_margin,
            "complete_aa_coverage": self.complete_aa_coverage,
        }


@dataclass
class ReferenceBundle:
    target: str
    ref_dir: str
    toml_path: str
    prot_faa: str
    config: dict
    segment_keys: list[str]
    gene_configs: dict
    protein_lengths: dict[str, int]


@dataclass
class CdsRow:
    start: int
    end: int
    strand: str


@dataclass
class CandidateHit:
    contig_id: str
    target: str
    parent_id: str
    product: str
    segment: str | None
    identity: float
    positive: float
    ref_aa_start: int
    ref_aa_end: int
    ref_aa_length: int
    aa_coverage: float
    query_start: int
    query_end: int
    query_coverage: float
    strand: str
    raw_score: float
    normalized_score: float
    frameshift_count: int = 0
    cds_count: int = 0
    internal_stop_count: int = 0
    missing_start_count: int = 0
    missing_stop_count: int = 0
    cds_length_not_mod3: bool = False
    touches_contig_edge: bool = False
    stop_codon_present: bool = False
    flags: list[str] = field(default_factory=list)

    @property
    def target_range(self) -> str:
        if not self.ref_aa_length:
            return "-"
        return f"{self.ref_aa_start}-{self.ref_aa_end}/{self.ref_aa_length}"

    @property
    def query_range(self) -> str:
        return f"{self.query_start}-{self.query_end}"


@dataclass
class AutoCall:
    contig_id: str
    length: int
    call: str
    target: str = "-"
    segment: str = "-"
    status: str = "no_hit"
    flags: list[str] = field(default_factory=list)
    confidence: str = "low"
    best_hit: CandidateHit | None = None
    second_target: str = "-"
    second_score: float = 0.0
    score_margin: float = 0.0
    notes: list[str] = field(default_factory=list)


def parse_auto_targets(value: str | None) -> list[str]:
    if not value:
        return list(DEFAULT_AUTO_TARGETS)
    targets = [part.strip().upper() for part in str(value).split(",") if part.strip()]
    unsupported = [target for target in targets if target not in DEFAULT_AUTO_TARGETS]
    if unsupported:
        raise ValueError(f"Unsupported auto target(s): {', '.join(unsupported)}")
    return list(dict.fromkeys(targets))


def resolve_target_reference_dir(target: str, db_dir: str | None) -> str:
    if db_dir:
        root = os.path.abspath(db_dir)
        direct_toml = os.path.join(root, f"{target}.toml")
        nested = os.path.join(root, target)
        if os.path.isfile(direct_toml):
            return root
        return nested
    return str(resources.files("ganflu").joinpath("db", target).resolve())


def load_reference_bundle(target: str, db_dir: str | None, logger) -> ReferenceBundle:
    ref_dir = resolve_target_reference_dir(target, db_dir)
    ref_toml = validate_reference_files.validate_reference_files(
        target=target, db_dir=ref_dir, logger=logger
    )
    toml_path = os.path.join(ref_dir, f"{target}.toml")
    prot_faa = os.path.join(ref_dir, ref_toml["metadata"]["prot_faa"])
    protein_lengths = {
        record.id: len(record.seq)
        for record in SeqIO.parse(prot_faa, "fasta")
    }
    return ReferenceBundle(
        target=target,
        ref_dir=ref_dir,
        toml_path=toml_path,
        prot_faa=prot_faa,
        config=ref_toml,
        segment_keys=list(ref_toml.get("segments", {}).keys()),
        gene_configs=ref_toml.get("genes", {}),
        protein_lengths=protein_lengths,
    )


def parse_gff3_attributes(attributes: str) -> dict[str, str]:
    if not attributes or attributes == ".":
        return {}
    parsed = {}
    for item in attributes.split(";"):
        if not item:
            continue
        key, value = item.split("=", 1)
        parsed[key] = value
    return parsed


def parse_target_attribute(value: str) -> tuple[str, int, int]:
    parts = value.split()
    product = parts[0]
    start = int(parts[1]) if len(parts) > 1 else 1
    end = int(parts[2]) if len(parts) > 2 else start
    return product, start, end


def parse_float(value: str | None, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_paf_tags(fields: list[str]) -> dict[str, str]:
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_paf_line(line: str) -> tuple[tuple[str, str, int, int], dict[str, int]]:
    fields = line.rstrip("\n").split("\t")
    product = fields[1]
    query_start = int(fields[3]) + 1
    query_end = int(fields[4])
    contig_id = fields[6]
    tags = parse_paf_tags(fields[13:])
    meta = {
        "frameshift_count": int(tags.get("fs", "0")),
        "raw_score": int(tags.get("AS", "0")),
    }
    return (contig_id, product, query_start, query_end), meta


def get_product_segment(product: str, segment_keys) -> str | None:
    return gff3togbk.product_to_segment(product, segment_keys)


def product_gene_name(product: str) -> str:
    return product.split("_")[0]


def build_cds_feature(cds_rows: list[CdsRow], product: str) -> SeqFeature | None:
    if not cds_rows:
        return None
    locations = [
        FeatureLocation(row.start - 1, row.end, strand=1 if row.strand == "+" else -1)
        for row in cds_rows
    ]
    location = locations[0] if len(locations) == 1 else CompoundLocation(locations)
    return SeqFeature(
        location=location,
        type="CDS",
        qualifiers={"gene": product_gene_name(product), "product": product},
    )


def get_translation_qc(
    contig_record: SeqRecord,
    product: str,
    cds_rows: list[CdsRow],
) -> tuple[int, int, int]:
    feature = build_cds_feature(cds_rows, product)
    if feature is None:
        return 0, 0, 0
    record = SeqRecord(
        contig_record.seq,
        id=contig_record.id,
        name=contig_record.name,
        description=contig_record.description,
        features=[feature],
    )
    gff3togbk.add_translations(record)
    checked = record.features[0]
    notes = checked.qualifiers.get("note", [])
    if isinstance(notes, str):
        notes = [notes]
    internal_stop_count = int(
        checked.type == "misc_feature"
        or any("nonfunctional due to mutation" in note for note in notes)
    )
    missing_start_count = int(
        any("start codon not found" in note for note in notes)
    )
    missing_stop_count = int(
        any("stop codon not found" in note for note in notes)
    )
    return internal_stop_count, missing_start_count, missing_stop_count


def collect_candidate_flags(
    candidate: CandidateHit,
    thresholds: AutoThresholds,
) -> list[str]:
    flags = []
    if candidate.identity < thresholds.min_identity:
        flags.append("low_identity")
    if candidate.aa_coverage < thresholds.min_aa_coverage:
        flags.append("low_aa_coverage")
    if candidate.normalized_score < thresholds.min_score:
        flags.append("low_score")
    if candidate.aa_coverage < thresholds.complete_aa_coverage:
        flags.append("partial_coverage")
    if candidate.ref_aa_start > 1:
        flags.append("reference_start_missing")
    if candidate.ref_aa_length and candidate.ref_aa_end < candidate.ref_aa_length:
        flags.append("reference_end_missing")
    if candidate.touches_contig_edge:
        flags.append("contig_edge")
    if candidate.frameshift_count > 0:
        flags.append("frameshift")
    if candidate.cds_length_not_mod3:
        flags.append("cds_length_not_mod3")
    if candidate.internal_stop_count > 0:
        flags.extend(["internal_stop", "nonfunctional"])
    if candidate.missing_start_count > 0:
        flags.append("missing_start")
    if candidate.missing_stop_count > 0 or not candidate.stop_codon_present:
        flags.append("missing_stop")
    return list(dict.fromkeys(flags))


def parse_miniprot_gff3(
    gff3_path: str,
    reference: ReferenceBundle,
    contigs_by_id: dict[str, SeqRecord],
    thresholds: AutoThresholds,
) -> list[CandidateHit]:
    paf_by_key = {}
    mrna_rows = []
    cds_by_parent = defaultdict(list)
    stop_codon_parents = set()

    with open(gff3_path, "r", encoding="utf-8") as handle:
        for line_index, line in enumerate(handle):
            if not line.strip():
                continue
            if line.startswith("##PAF\t"):
                key, meta = parse_paf_line(line)
                paf_by_key[key] = meta
                continue
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                continue
            seqid, _, feature_type, start, end, score, strand, _, attributes = columns
            attrs = parse_gff3_attributes(attributes)
            if feature_type == "mRNA" and "Target" in attrs:
                product, ref_start, ref_end = parse_target_attribute(attrs["Target"])
                mrna_rows.append(
                    {
                        "contig_id": seqid,
                        "parent_id": attrs.get("ID", ""),
                        "product": product,
                        "ref_start": ref_start,
                        "ref_end": ref_end,
                        "query_start": int(start),
                        "query_end": int(end),
                        "raw_score": parse_float(score),
                        "strand": strand,
                        "identity": parse_float(attrs.get("Identity")),
                        "positive": parse_float(attrs.get("Positive")),
                        "line_index": line_index,
                    }
                )
            elif feature_type == "CDS" and "Parent" in attrs:
                cds_by_parent[attrs["Parent"]].append(
                    CdsRow(start=int(start), end=int(end), strand=strand)
                )
            elif feature_type == "stop_codon" and "Parent" in attrs:
                stop_codon_parents.add(attrs["Parent"])

    candidates = []
    for row in mrna_rows:
        contig_record = contigs_by_id.get(row["contig_id"])
        if contig_record is None:
            continue
        product = row["product"]
        ref_length = reference.protein_lengths.get(product, row["ref_end"])
        aa_coverage = 0.0
        if ref_length:
            aa_coverage = max(0, row["ref_end"] - row["ref_start"] + 1) / ref_length
        query_length = max(0, row["query_end"] - row["query_start"] + 1)
        query_coverage = query_length / len(contig_record.seq) if contig_record.seq else 0.0
        paf_meta = paf_by_key.get(
            (row["contig_id"], product, row["ref_start"], row["ref_end"]),
            {},
        )
        raw_score = row["raw_score"] or float(paf_meta.get("raw_score", 0))
        cds_rows = cds_by_parent.get(row["parent_id"], [])
        cds_length = sum(cds.end - cds.start + 1 for cds in cds_rows)
        gene_config = reference.gene_configs.get(product_gene_name(product), {})
        cds_length_not_mod3 = bool(
            cds_length
            and cds_length % 3 != 0
            and not gene_config.get("ribosomal_slippage")
        )
        internal_stop_count, missing_start_count, missing_stop_count = get_translation_qc(
            contig_record, product, cds_rows
        )
        touches_edge = any(
            cds.start <= 1 or cds.end >= len(contig_record.seq)
            for cds in cds_rows
        )
        normalized_score = row["identity"] * aa_coverage
        candidate = CandidateHit(
            contig_id=row["contig_id"],
            target=reference.target,
            parent_id=row["parent_id"],
            product=product,
            segment=get_product_segment(product, reference.segment_keys),
            identity=row["identity"],
            positive=row["positive"],
            ref_aa_start=row["ref_start"],
            ref_aa_end=row["ref_end"],
            ref_aa_length=ref_length,
            aa_coverage=aa_coverage,
            query_start=row["query_start"],
            query_end=row["query_end"],
            query_coverage=query_coverage,
            strand=row["strand"],
            raw_score=raw_score,
            normalized_score=normalized_score,
            frameshift_count=int(paf_meta.get("frameshift_count", 0)),
            cds_count=len(cds_rows),
            internal_stop_count=internal_stop_count,
            missing_start_count=missing_start_count,
            missing_stop_count=missing_stop_count,
            cds_length_not_mod3=cds_length_not_mod3,
            touches_contig_edge=touches_edge,
            stop_codon_present=row["parent_id"] in stop_codon_parents,
        )
        candidate.flags = collect_candidate_flags(candidate, thresholds)
        candidates.append(candidate)
    adjust_ribosomal_slippage_fragment_qc(
        candidates,
        cds_by_parent,
        contigs_by_id,
        reference,
        thresholds,
    )
    return candidates


def is_ribosomal_slippage_fragment(candidate: CandidateHit, reference: ReferenceBundle) -> bool:
    gene_name = product_gene_name(candidate.product)
    gene_config = reference.gene_configs.get(gene_name, {})
    return bool(gene_config.get("ribosomal_slippage") and candidate.product != gene_name)


def adjust_ribosomal_slippage_fragment_qc(
    candidates: list[CandidateHit],
    cds_by_parent: dict[str, list[CdsRow]],
    contigs_by_id: dict[str, SeqRecord],
    reference: ReferenceBundle,
    thresholds: AutoThresholds,
) -> None:
    grouped = defaultdict(list)
    for candidate in candidates:
        if is_ribosomal_slippage_fragment(candidate, reference):
            grouped[(candidate.contig_id, product_gene_name(candidate.product))].append(candidate)

    for (contig_id, gene_name), fragment_candidates in grouped.items():
        if len(fragment_candidates) < 2:
            continue
        contig_record = contigs_by_id.get(contig_id)
        if contig_record is None:
            continue
        group_cds_rows = []
        for candidate in fragment_candidates:
            group_cds_rows.extend(cds_by_parent.get(candidate.parent_id, []))
        if not group_cds_rows:
            continue

        internal_stop_count, missing_start_count, missing_stop_count = get_translation_qc(
            contig_record,
            gene_name,
            group_cds_rows,
        )
        group_has_stop = missing_stop_count == 0
        for candidate in fragment_candidates:
            candidate.internal_stop_count = internal_stop_count
            candidate.missing_start_count = missing_start_count
            candidate.missing_stop_count = missing_stop_count
            candidate.stop_codon_present = group_has_stop
            candidate.flags = collect_candidate_flags(candidate, thresholds)


def candidate_sort_key(candidate: CandidateHit):
    return (
        candidate.normalized_score,
        candidate.aa_coverage,
        candidate.identity,
        candidate.raw_score,
    )


PRIMARY_SEGMENT_PRODUCTS = {
    "M": {"M1", "P42"},
    "NS": {"NS1"},
}


def is_segment_representative_candidate(candidate: CandidateHit) -> bool:
    segment = candidate.segment
    if not segment:
        return False
    if candidate.product == segment:
        return True
    if candidate.product.startswith(f"{segment}_"):
        return True
    return product_gene_name(candidate.product) in PRIMARY_SEGMENT_PRODUCTS.get(segment, set())


def passes_auto_thresholds(candidate: CandidateHit, thresholds: AutoThresholds) -> bool:
    return (
        candidate.identity >= thresholds.min_identity
        and candidate.aa_coverage >= thresholds.min_aa_coverage
        and candidate.normalized_score >= thresholds.min_score
    )


def choose_report_candidate(
    candidates: list[CandidateHit],
    best: CandidateHit,
    thresholds: AutoThresholds,
) -> CandidateHit:
    representatives = [
        candidate
        for candidate in candidates
        if (
            candidate.target == best.target
            and candidate.segment == best.segment
            and is_segment_representative_candidate(candidate)
            and passes_auto_thresholds(candidate, thresholds)
        )
    ]
    if representatives:
        return max(representatives, key=candidate_sort_key)
    return best


def best_by_target_and_segment(candidates: list[CandidateHit]) -> dict[str, dict[str, CandidateHit]]:
    collapsed = defaultdict(dict)
    for candidate in candidates:
        if not candidate.segment:
            continue
        current = collapsed[candidate.target].get(candidate.segment)
        if current is None or candidate_sort_key(candidate) > candidate_sort_key(current):
            collapsed[candidate.target][candidate.segment] = candidate
    return collapsed


def determine_status(candidate: CandidateHit, thresholds: AutoThresholds) -> str:
    if candidate.frameshift_count > 0 or candidate.cds_length_not_mod3:
        return "frameshift"
    if candidate.internal_stop_count > 0:
        return "nonfunctional"
    partial_flags = {
        "partial_coverage",
        "reference_start_missing",
        "reference_end_missing",
        "missing_start",
        "missing_stop",
    }
    if candidate.aa_coverage < thresholds.complete_aa_coverage:
        return "partial"
    if any(flag in partial_flags for flag in candidate.flags):
        return "partial"
    return "complete"


def determine_confidence(
    call: str,
    status: str,
    best_score: float,
    margin: float,
    thresholds: AutoThresholds,
) -> str:
    if call != "accept":
        return "low"
    if (
        status == "complete"
        and best_score >= 0.80
        and margin >= thresholds.min_margin * 2
    ):
        return "high"
    return "medium"


def classify_contig(
    contig: SeqRecord,
    candidates: list[CandidateHit],
    thresholds: AutoThresholds,
) -> AutoCall:
    if not candidates:
        return AutoCall(
            contig_id=contig.id,
            length=len(contig.seq),
            call="reject",
            status="no_hit",
            notes=["no miniprot mRNA hit"],
        )

    collapsed = best_by_target_and_segment(candidates)
    target_best = {}
    target_second_segment = {}
    for target, by_segment in collapsed.items():
        ranked = sorted(by_segment.values(), key=candidate_sort_key, reverse=True)
        if ranked:
            target_best[target] = ranked[0]
        if len(ranked) > 1:
            target_second_segment[target] = ranked[1]

    if not target_best:
        return AutoCall(
            contig_id=contig.id,
            length=len(contig.seq),
            call="reject",
            status="no_hit",
            notes=["no segment-compatible miniprot hit"],
        )

    ranked_targets = sorted(target_best.values(), key=candidate_sort_key, reverse=True)
    best = ranked_targets[0]
    second = ranked_targets[1] if len(ranked_targets) > 1 else None
    second_target = second.target if second else "-"
    second_score = second.normalized_score if second else 0.0
    margin = best.normalized_score - second_score if second else best.normalized_score

    if (
        best.identity < thresholds.min_identity
        or best.aa_coverage < thresholds.min_aa_coverage
        or best.normalized_score < thresholds.min_score
    ):
        return AutoCall(
            contig_id=contig.id,
            length=len(contig.seq),
            call="reject",
            target="-",
            segment="-",
            status="off_target",
            flags=best.flags,
            confidence="low",
            best_hit=best,
            second_target=second_target,
            second_score=second_score,
            score_margin=margin,
            notes=["best hit below auto thresholds"],
        )

    if second and margin < thresholds.min_margin:
        return AutoCall(
            contig_id=contig.id,
            length=len(contig.seq),
            call="reject",
            target=best.target,
            segment=best.segment or "-",
            status="ambiguous",
            flags=list(dict.fromkeys([*best.flags, "target_ambiguous"])),
            confidence="low",
            best_hit=best,
            second_target=second_target,
            second_score=second_score,
            score_margin=margin,
            notes=["best target and second target scores are too close"],
        )

    second_segment = target_second_segment.get(best.target)
    report_best = choose_report_candidate(candidates, best, thresholds)
    status = determine_status(report_best, thresholds)
    flags = report_best.flags
    notes = []
    if second_segment and best.normalized_score - second_segment.normalized_score < thresholds.min_margin:
        flags = list(dict.fromkeys([*flags, "segment_close_score"]))
        notes = ["best segment and second segment scores are close"]
    return AutoCall(
        contig_id=contig.id,
        length=len(contig.seq),
        call="accept",
        target=report_best.target,
        segment=report_best.segment or "-",
        status=status,
        flags=flags,
        confidence=determine_confidence(
            "accept", status, report_best.normalized_score, margin, thresholds
        ),
        best_hit=report_best,
        second_target=second_target,
        second_score=second_score,
        score_margin=margin,
        notes=notes,
    )


def classify_contigs(
    contigs: list[SeqRecord],
    candidates_by_contig: dict[str, list[CandidateHit]],
    thresholds: AutoThresholds,
) -> list[AutoCall]:
    return [
        classify_contig(contig, candidates_by_contig.get(contig.id, []), thresholds)
        for contig in contigs
    ]


def format_float(value: float) -> str:
    return f"{value:.4f}"


def auto_call_to_row(call: AutoCall) -> dict[str, str | int]:
    best = call.best_hit
    return {
        "contig_id": call.contig_id,
        "length": call.length,
        "call": call.call,
        "target": call.target,
        "segment": call.segment,
        "status": call.status,
        "flags": ";".join(call.flags) if call.flags else "-",
        "confidence": call.confidence,
        "best_identity": format_float(best.identity) if best else "0.0000",
        "best_aa_coverage": format_float(best.aa_coverage) if best else "0.0000",
        "best_score": format_float(best.normalized_score) if best else "0.0000",
        "second_target": call.second_target,
        "second_score": format_float(call.second_score),
        "score_margin": format_float(call.score_margin),
        "best_product": best.product if best else "-",
        "best_target_range": best.target_range if best else "-",
        "best_query_range": best.query_range if best else "-",
        "strand": best.strand if best else "-",
        "cds_count": best.cds_count if best else 0,
        "frameshift_count": best.frameshift_count if best else 0,
        "internal_stop_count": best.internal_stop_count if best else 0,
        "missing_start_count": best.missing_start_count if best else 0,
        "missing_stop_count": best.missing_stop_count if best else 0,
        "notes": ";".join(call.notes) if call.notes else "-",
    }


def write_auto_tsv(calls: list[AutoCall], path: str) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        for call in calls:
            writer.writerow(auto_call_to_row(call))


def make_accepted_segments(calls: list[AutoCall]) -> dict[str, dict[str, str]]:
    accepted = defaultdict(dict)
    for call in calls:
        if call.call == "accept" and call.target != "-" and call.segment != "-":
            accepted[call.target][call.contig_id] = call.segment
    return {target: dict(value) for target, value in accepted.items()}


def paf_line_matches(line: str, accepted_segments: dict[str, str], segment_keys) -> bool:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 7:
        return False
    product = fields[1]
    contig_id = fields[6]
    segment = get_product_segment(product, segment_keys)
    return accepted_segments.get(contig_id) == segment


def filter_gff3_for_target(
    scan_gff3: str,
    output_gff3: str,
    accepted_segments: dict[str, str],
    reference: ReferenceBundle,
) -> None:
    gff3_prune.prune_gff3(
        scan_gff3,
        output_gff3,
        protein_lengths=reference.protein_lengths,
        antigen_names=reference.config.get("serotype", {}).keys(),
        accepted_segments=accepted_segments,
        segment_keys=reference.segment_keys,
        product_to_segment=gff3togbk.product_to_segment,
    )


def write_target_fasta(
    contigs: list[SeqRecord],
    accepted_segments: dict[str, str],
    output_fasta: str,
) -> None:
    with open(output_fasta, "w", encoding="utf-8") as handle:
        SeqIO.write(
            [record for record in contigs if record.id in accepted_segments],
            handle,
            "fasta",
        )


def run_annotation_for_targets(
    *,
    contigs: list[SeqRecord],
    accepted_by_target: dict[str, dict[str, str]],
    scan_gff3_by_target: dict[str, str],
    references: dict[str, ReferenceBundle],
    output_stem: str,
    auto_work_dir: str,
    isolate: str,
    preserve_original_id: bool,
    logger,
) -> dict[str, str]:
    outputs = {}
    for target in sorted(accepted_by_target):
        accepted_segments = accepted_by_target[target]
        if not accepted_segments:
            continue
        reference = references[target]
        target_stem = f"{output_stem}.{target}"
        target_fasta = os.path.join(
            auto_work_dir, f"{Path(output_stem).name}.{target}.accepted.fasta"
        )
        target_gff3 = f"{target_stem}.gff3"
        target_gbk = f"{target_stem}.gbk"
        target_cds = f"{target_stem}.cds.fna"
        target_faa = f"{target_stem}.faa"

        logger.info(
            f"Annotating {len(accepted_segments)} accepted contig(s) as {target}"
        )
        write_target_fasta(contigs, accepted_segments, target_fasta)
        filter_gff3_for_target(
            scan_gff3_by_target[target],
            target_gff3,
            accepted_segments,
            reference,
        )

        gff3togbk_args = [
            "-g", target_gff3,
            "-o", target_gbk,
            "-i", target_fasta,
            "--toml", reference.toml_path,
            "--isolate", isolate,
            "--cds-fna", target_cds,
            "--faa", target_faa,
        ]
        if preserve_original_id:
            gff3togbk_args.append("--preserve_original_id")
        gff3togbk.main(gff3togbk_args)

        outputs[f"{target}.gff3"] = target_gff3
        outputs[f"{target}.gbk"] = target_gbk
        outputs[f"{target}.cds_fna"] = target_cds
        outputs[f"{target}.faa"] = target_faa
    return outputs


def write_rejected_fasta(
    contigs: list[SeqRecord],
    calls: list[AutoCall],
    output_path: str,
) -> None:
    rejected_ids = {
        call.contig_id
        for call in calls
        if call.call == "reject"
    }
    with open(output_path, "w", encoding="utf-8") as handle:
        SeqIO.write(
            [record for record in contigs if record.id in rejected_ids],
            handle,
            "fasta",
        )


def build_summary(
    *,
    input_fasta: str,
    output_stem: str,
    targets: list[str],
    thresholds: AutoThresholds,
    calls: list[AutoCall],
    outputs: dict[str, str],
) -> dict:
    call_counts = Counter(call.call for call in calls)
    status_counts = Counter(call.status for call in calls)
    by_target = {}
    for call in calls:
        if call.target == "-":
            continue
        target_summary = by_target.setdefault(
            call.target,
            {
                "accepted": 0,
                "rejected": 0,
                "complete": 0,
                "partial": 0,
                "frameshift": 0,
                "nonfunctional": 0,
                "ambiguous": 0,
                "segments": {},
            },
        )
        if call.call == "accept":
            target_summary["accepted"] += 1
        else:
            target_summary["rejected"] += 1
        if call.status in target_summary:
            target_summary[call.status] += 1
        if call.segment != "-":
            target_summary["segments"][call.segment] = (
                target_summary["segments"].get(call.segment, 0) + 1
            )

    return {
        "input": os.path.abspath(input_fasta),
        "output_stem": os.path.abspath(output_stem),
        "targets_scanned": targets,
        "thresholds": thresholds.as_dict(),
        "counts": {
            "input_contigs": len(calls),
            "accepted": call_counts.get("accept", 0),
            "rejected": call_counts.get("reject", 0),
        },
        "by_target": by_target,
        "by_status": dict(status_counts),
        "outputs": outputs,
    }


def write_summary_json(summary: dict, path: str) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")


def log_auto_summary(calls: list[AutoCall], logger) -> None:
    call_counts = Counter(call.call for call in calls)
    status_counts = Counter(call.status for call in calls)
    logger.info(
        "Auto mode calls: "
        + ", ".join(f"{key}={value}" for key, value in sorted(call_counts.items()))
    )
    logger.info(
        "Auto mode statuses: "
        + ", ".join(f"{key}={value}" for key, value in sorted(status_counts.items()))
    )
    target_counts = Counter(
        call.target for call in calls if call.call == "accept" and call.target != "-"
    )
    if target_counts:
        logger.info(
            "Accepted targets: "
            + ", ".join(f"{key}={value}" for key, value in sorted(target_counts.items()))
        )
    else:
        logger.warning("No contigs were accepted for annotation")


def run_auto(args, output_stem: str, work_dir: str, logger) -> dict:
    thresholds = AutoThresholds.from_args(args)
    targets = parse_auto_targets(args.auto_targets)
    report_stem = os.path.abspath(args.auto_report_prefix) if args.auto_report_prefix else output_stem
    auto_work_dir = f"{report_stem}.auto.work"
    os.makedirs(auto_work_dir, exist_ok=True)

    logger.info("Auto mode started")
    logger.info(f"Auto targets: {', '.join(targets)}")
    logger.info(f"Auto thresholds: {thresholds.as_dict()}")
    logger.info(f"Auto work directory: {auto_work_dir}")

    contigs = list(SeqIO.parse(args.input, "fasta"))
    if not contigs:
        raise ValueError("Input FASTA contains no records")
    contigs_by_id = {record.id: record for record in contigs}
    if len(contigs_by_id) != len(contigs):
        raise ValueError("Input FASTA contains duplicate record IDs, which auto mode cannot disambiguate")

    references = {
        target: load_reference_bundle(target, args.db_dir, logger)
        for target in targets
    }

    scan_gff3_by_target = {}
    candidates_by_contig = defaultdict(list)
    for target in targets:
        reference = references[target]
        scan_gff3 = os.path.join(auto_work_dir, f"{Path(report_stem).name}.{target}.scan.gff3")
        scan_gff3_by_target[target] = scan_gff3
        logger.info(f"Running miniprot auto scan for {target}")
        miniprot = MiniprotCommandLine(
            input=args.input,
            work_dir=auto_work_dir,
            output=scan_gff3,
            prot_faa=reference.prot_faa,
            miniprot_bin="miniprot",
            stderr_filename=f"{target}.miniprot.stderr",
            kmer_size=15,
            prefix=MINIPROT_PREFIXES.get(target, "MP"),
            max_secondary_alignments=gff3_prune.RELAXED_MAX_SECONDARY_ALIGNMENTS,
            secondary_to_primary_ratio=gff3_prune.RELAXED_SECONDARY_TO_PRIMARY_RATIO,
            output_score_ratio=gff3_prune.RELAXED_OUTPUT_SCORE_RATIO,
        )
        miniprot.run_piped_commands()
        candidates = parse_miniprot_gff3(
            scan_gff3,
            reference,
            contigs_by_id,
            thresholds,
        )
        logger.info(f"{target} auto scan candidates: {len(candidates)}")
        for candidate in candidates:
            logger.debug(
                "Auto candidate %s %s %s %s identity=%.4f coverage=%.4f score=%.4f flags=%s",
                candidate.contig_id,
                candidate.target,
                candidate.segment,
                candidate.product,
                candidate.identity,
                candidate.aa_coverage,
                candidate.normalized_score,
                ";".join(candidate.flags) or "-",
            )
            candidates_by_contig[candidate.contig_id].append(candidate)

    calls = classify_contigs(contigs, candidates_by_contig, thresholds)
    accepted_by_target = make_accepted_segments(calls)
    outputs = run_annotation_for_targets(
        contigs=contigs,
        accepted_by_target=accepted_by_target,
        scan_gff3_by_target=scan_gff3_by_target,
        references=references,
        output_stem=output_stem,
        auto_work_dir=auto_work_dir,
        isolate=args.isolate,
        preserve_original_id=args.preserve_original_id,
        logger=logger,
    )

    tsv_path = f"{report_stem}.auto.tsv"
    summary_path = f"{report_stem}.auto.summary.json"
    write_auto_tsv(calls, tsv_path)
    outputs["auto.tsv"] = tsv_path

    if args.auto_write_rejected:
        rejected_path = f"{report_stem}.auto.rejected.fasta"
        write_rejected_fasta(contigs, calls, rejected_path)
        outputs["auto.rejected_fasta"] = rejected_path

    outputs["auto.summary_json"] = summary_path
    summary = build_summary(
        input_fasta=args.input,
        output_stem=output_stem,
        targets=targets,
        thresholds=thresholds,
        calls=calls,
        outputs=outputs,
    )
    write_summary_json(summary, summary_path)

    log_auto_summary(calls, logger)
    logger.info(f"Auto TSV output: {tsv_path}")
    logger.info(f"Auto summary JSON output: {summary_path}")
    return summary
