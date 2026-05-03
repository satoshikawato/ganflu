#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO


RELAXED_MAX_SECONDARY_ALIGNMENTS = 100
RELAXED_SECONDARY_TO_PRIMARY_RATIO = 0.10
RELAXED_OUTPUT_SCORE_RATIO = 0.10


@dataclass
class Gff3Row:
    line_index: int
    line: str
    kind: str
    columns: list[str] = field(default_factory=list)
    attributes: dict[str, str] = field(default_factory=dict)
    paf_key: tuple[str, str, int, int] | None = None
    paf_target_length: int = 0
    paf_raw_score: float = 0.0

    @property
    def seqid(self) -> str:
        return self.columns[0]

    @property
    def feature_type(self) -> str:
        return self.columns[2]

    @property
    def start(self) -> int:
        return int(self.columns[3])

    @property
    def end(self) -> int:
        return int(self.columns[4])

    @property
    def score(self) -> float:
        return parse_float(self.columns[5])

    @property
    def strand(self) -> str:
        return self.columns[6]


@dataclass
class ParentAlignment:
    parent_id: str
    row: Gff3Row
    product: str
    target_start: int
    target_end: int
    target_length: int
    identity: float
    raw_score: float
    rank: int
    children: list[Gff3Row] = field(default_factory=list)

    @property
    def seqid(self) -> str:
        return self.row.seqid

    @property
    def product_name(self) -> str:
        return self.product.split("_")[0]

    @property
    def aa_coverage(self) -> float:
        if not self.target_length:
            return 0.0
        return max(0, self.target_end - self.target_start + 1) / self.target_length

    @property
    def query_span(self) -> int:
        return abs(self.row.end - self.row.start) + 1

    @property
    def paf_key(self) -> tuple[str, str, int, int]:
        return (self.seqid, self.product, self.target_start, self.target_end)

    @property
    def sort_key(self) -> tuple[float, float, float, int, int]:
        return (
            self.identity,
            self.aa_coverage,
            self.raw_score,
            self.query_span,
            -self.rank,
        )


@dataclass
class PruneResult:
    selected_parent_count: int
    kept_feature_count: int
    kept_paf_count: int


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
        if value in {None, "."}:
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value: str | None, default: int = 0) -> int:
    try:
        if value in {None, "."}:
            return default
        return int(value)
    except (TypeError, ValueError):
        return default


def parse_paf_tags(fields: list[str]) -> dict[str, str]:
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_paf_row(line: str) -> tuple[tuple[str, str, int, int], int, float] | None:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 7:
        return None
    try:
        product = fields[1]
        target_length = int(fields[2])
        target_start = int(fields[3]) + 1
        target_end = int(fields[4])
        contig_id = fields[6]
    except (ValueError, IndexError):
        return None
    tags = parse_paf_tags(fields[13:])
    return (contig_id, product, target_start, target_end), target_length, parse_float(tags.get("AS"))


def split_parent_ids(value: str | None) -> list[str]:
    if not value:
        return []
    return [parent_id for parent_id in value.split(",") if parent_id]


def load_protein_lengths(prot_faa: str | Path | None) -> dict[str, int]:
    if not prot_faa:
        return {}
    return {record.id: len(record.seq) for record in SeqIO.parse(str(prot_faa), "fasta")}


def read_gff3_rows(gff3_path: str | Path) -> tuple[list[Gff3Row], dict[tuple[str, str, int, int], Gff3Row]]:
    rows = []
    paf_rows_by_key = {}
    with open(gff3_path, "r", encoding="utf-8") as handle:
        for line_index, line in enumerate(handle):
            if not line.strip():
                continue
            if line.startswith("##PAF\t"):
                row = Gff3Row(line_index=line_index, line=line, kind="paf")
                parsed_paf = parse_paf_row(line)
                if parsed_paf is not None:
                    paf_key, target_length, raw_score = parsed_paf
                    row.paf_key = paf_key
                    row.paf_target_length = target_length
                    row.paf_raw_score = raw_score
                    paf_rows_by_key[paf_key] = row
                rows.append(row)
                continue
            if line.startswith("#"):
                rows.append(Gff3Row(line_index=line_index, line=line, kind="comment"))
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                continue
            rows.append(
                Gff3Row(
                    line_index=line_index,
                    line=line,
                    kind="feature",
                    columns=columns,
                    attributes=parse_gff3_attributes(columns[8]),
                )
            )
    return rows, paf_rows_by_key


def collect_parent_alignments(
    rows: list[Gff3Row],
    paf_rows_by_key: dict[tuple[str, str, int, int], Gff3Row],
    protein_lengths: dict[str, int] | None = None,
) -> list[ParentAlignment]:
    protein_lengths = protein_lengths or {}
    children_by_parent = defaultdict(list)
    mrna_rows = []
    for row in rows:
        if row.kind != "feature":
            continue
        if row.feature_type == "mRNA" and "ID" in row.attributes and "Target" in row.attributes:
            mrna_rows.append(row)
        for parent_id in split_parent_ids(row.attributes.get("Parent")):
            children_by_parent[parent_id].append(row)

    parents = []
    for row in mrna_rows:
        parent_id = row.attributes["ID"]
        product, target_start, target_end = parse_target_attribute(row.attributes["Target"])
        paf_row = paf_rows_by_key.get((row.seqid, product, target_start, target_end))
        target_length = protein_lengths.get(product)
        if not target_length and paf_row is not None:
            target_length = paf_row.paf_target_length
        if not target_length:
            target_length = target_end
        raw_score = row.score
        if not raw_score and paf_row is not None:
            raw_score = paf_row.paf_raw_score
        parents.append(
            ParentAlignment(
                parent_id=parent_id,
                row=row,
                product=product,
                target_start=target_start,
                target_end=target_end,
                target_length=target_length,
                identity=parse_float(row.attributes.get("Identity")),
                raw_score=raw_score,
                rank=parse_int(row.attributes.get("Rank"), default=1_000_000),
                children=children_by_parent.get(parent_id, []),
            )
        )
    return parents


def choose_best(current: ParentAlignment | None, candidate: ParentAlignment) -> ParentAlignment:
    if current is None or candidate.sort_key > current.sort_key:
        return candidate
    return current


def select_parent_alignments(
    parents: list[ParentAlignment],
    *,
    antigen_names: set[str] | None = None,
    accepted_segments: dict[str, str] | None = None,
    segment_keys=None,
    product_to_segment=None,
) -> list[ParentAlignment]:
    antigen_names = antigen_names or set()
    selected_by_group = {}
    for parent in parents:
        if accepted_segments is not None:
            if parent.seqid not in accepted_segments:
                continue
            if product_to_segment is None:
                raise ValueError("product_to_segment is required when accepted_segments is provided")
            segment = product_to_segment(parent.product, segment_keys)
            if segment != accepted_segments[parent.seqid]:
                continue

        if parent.product_name in antigen_names:
            group_key = (parent.seqid, "antigen", parent.product_name)
        else:
            group_key = (parent.seqid, "product", parent.product)
        selected_by_group[group_key] = choose_best(selected_by_group.get(group_key), parent)

    selected = list(selected_by_group.values())
    selected.sort(key=lambda parent: parent.row.line_index)
    return selected


def prune_gff3(
    input_gff3: str | Path,
    output_gff3: str | Path,
    *,
    protein_lengths: dict[str, int] | None = None,
    prot_faa: str | Path | None = None,
    antigen_names=None,
    accepted_segments: dict[str, str] | None = None,
    segment_keys=None,
    product_to_segment=None,
) -> PruneResult:
    if protein_lengths is None:
        protein_lengths = load_protein_lengths(prot_faa)
    antigen_set = set(antigen_names or ())

    rows, paf_rows_by_key = read_gff3_rows(input_gff3)
    parents = collect_parent_alignments(rows, paf_rows_by_key, protein_lengths)
    selected = select_parent_alignments(
        parents,
        antigen_names=antigen_set,
        accepted_segments=accepted_segments,
        segment_keys=segment_keys,
        product_to_segment=product_to_segment,
    )

    selected_parent_ids = {parent.parent_id for parent in selected}
    selected_paf_keys = {parent.paf_key for parent in selected}

    kept_feature_count = 0
    kept_paf_count = 0
    written_paf_keys = set()
    with open(output_gff3, "w", encoding="utf-8") as output:
        output.write("##gff-version 3\n")
        for row in rows:
            if row.kind == "comment":
                continue
            if row.kind == "paf":
                if row.paf_key in selected_paf_keys and row.paf_key not in written_paf_keys:
                    output.write(row.line)
                    written_paf_keys.add(row.paf_key)
                    kept_paf_count += 1
                continue
            if row.kind != "feature":
                continue

            keep = False
            if row.feature_type == "mRNA" and row.attributes.get("ID") in selected_parent_ids:
                keep = True
            elif any(parent_id in selected_parent_ids for parent_id in split_parent_ids(row.attributes.get("Parent"))):
                keep = True

            if keep:
                output.write(row.line)
                kept_feature_count += 1

    return PruneResult(
        selected_parent_count=len(selected_parent_ids),
        kept_feature_count=kept_feature_count,
        kept_paf_count=kept_paf_count,
    )
