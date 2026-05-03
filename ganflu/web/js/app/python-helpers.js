export const PYTHON_HELPERS = `
import contextlib
import glob
import io
import json
import logging
import os
import re
import traceback
from collections import Counter, defaultdict
from importlib import resources
try:
    import tomllib
except ModuleNotFoundError:
    tomllib = None
    import toml

from Bio import SeqIO
from ganflu.scripts import auto_mode, gff3togbk

SUPPORTED_TARGETS = {"IAV", "IBV", "ICV", "IDV"}
AUTO_TARGETS = ("IAV", "IBV", "ICV", "IDV")
WEB_HIT_SETTING_DEFAULTS = {
    "min_identity": 0.70,
    "min_aa_coverage": 0.35,
    "min_score": 0.25,
    "min_margin": 0.10,
    "complete_aa_coverage": 0.90,
    "max_secondary_alignments": 100,
    "output_score_ratio": 0.10,
    "secondary_to_primary_ratio": 0.10,
}
WEB_HIT_SETTING_ALIASES = {
    "min_identity": ("min_identity", "minIdentity"),
    "min_aa_coverage": ("min_aa_coverage", "minAaCoverage"),
    "min_score": ("min_score", "minScore"),
    "min_margin": ("min_margin", "minMargin"),
    "complete_aa_coverage": ("complete_aa_coverage", "completeAaCoverage"),
    "max_secondary_alignments": ("max_secondary_alignments", "maxSecondaryAlignments"),
    "output_score_ratio": ("output_score_ratio", "outputScoreRatio"),
    "secondary_to_primary_ratio": ("secondary_to_primary_ratio", "secondaryToPrimaryRatio"),
}

def _safe_stem(value):
    value = str(value or "ganflu").strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    return value or "ganflu"

def _read_resource_text(path):
    return path.read_text(encoding="utf-8")

def _load_toml_text(text):
    if tomllib is not None:
        return tomllib.loads(text)
    return toml.loads(text)

def _load_hit_settings(value=None):
    if value is None:
        return {}
    if isinstance(value, str):
        if not value.strip():
            return {}
        loaded = json.loads(value)
    elif isinstance(value, dict):
        loaded = value
    else:
        loaded = dict(value or {})
    if loaded is None:
        return {}
    if not isinstance(loaded, dict):
        raise ValueError("Hit settings must be a JSON object.")
    return loaded

def _raw_hit_setting(raw_settings, name):
    for key in WEB_HIT_SETTING_ALIASES[name]:
        if key in raw_settings and raw_settings[key] not in (None, ""):
            return raw_settings[key]
    return WEB_HIT_SETTING_DEFAULTS[name]

def _coerce_float_setting(name, value, min_value=0.0, max_value=1.0):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name.replace('_', ' ')} must be a number")
    if number < min_value or number > max_value:
        raise ValueError(f"{name.replace('_', ' ')} must be between {min_value} and {max_value}")
    return number

def _coerce_int_setting(name, value, min_value=1):
    if isinstance(value, bool):
        raise ValueError(f"{name.replace('_', ' ')} must be a whole number")
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name.replace('_', ' ')} must be a whole number")
    if not number.is_integer():
        raise ValueError(f"{name.replace('_', ' ')} must be a whole number")
    integer = int(number)
    if integer < min_value:
        raise ValueError(f"{name.replace('_', ' ')} must be at least {min_value}")
    return integer

def _normalize_hit_settings(value=None):
    raw_settings = _load_hit_settings(value)
    return {
        "min_identity": _coerce_float_setting("min_identity", _raw_hit_setting(raw_settings, "min_identity")),
        "min_aa_coverage": _coerce_float_setting("min_aa_coverage", _raw_hit_setting(raw_settings, "min_aa_coverage")),
        "min_score": _coerce_float_setting("min_score", _raw_hit_setting(raw_settings, "min_score")),
        "min_margin": _coerce_float_setting("min_margin", _raw_hit_setting(raw_settings, "min_margin")),
        "complete_aa_coverage": _coerce_float_setting("complete_aa_coverage", _raw_hit_setting(raw_settings, "complete_aa_coverage")),
        "max_secondary_alignments": _coerce_int_setting("max_secondary_alignments", _raw_hit_setting(raw_settings, "max_secondary_alignments")),
        "output_score_ratio": _coerce_float_setting("output_score_ratio", _raw_hit_setting(raw_settings, "output_score_ratio")),
        "secondary_to_primary_ratio": _coerce_float_setting("secondary_to_primary_ratio", _raw_hit_setting(raw_settings, "secondary_to_primary_ratio")),
    }

def _thresholds_from_hit_settings(settings):
    return auto_mode.AutoThresholds(
        min_identity=settings["min_identity"],
        min_aa_coverage=settings["min_aa_coverage"],
        min_score=settings["min_score"],
        min_margin=settings["min_margin"],
        complete_aa_coverage=settings["complete_aa_coverage"],
    )

def _miniprot_settings_for_summary(settings):
    return {
        "max_secondary_alignments": settings["max_secondary_alignments"],
        "output_score_ratio": settings["output_score_ratio"],
        "secondary_to_primary_ratio": settings["secondary_to_primary_ratio"],
    }

def _resource_join(root, relative_path):
    return root.joinpath(*str(relative_path).split("/"))

def get_ganflu_reference_json(target):
    target = str(target or "").upper()
    if target not in SUPPORTED_TARGETS:
        return json.dumps({"error": f"Unsupported target: {target}"})
    try:
        target_root = resources.files("ganflu").joinpath("db", target)
        prot_path = target_root.joinpath("prot", f"{target}_proteome_consensus.faa")
        toml_path = target_root.joinpath(f"{target}.toml")
        return json.dumps(
            {
                "target": target,
                "protein_fasta": _read_resource_text(prot_path),
                "toml": _read_resource_text(toml_path),
            }
        )
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

def _count_genbank_outputs(gbk_path):
    records = list(SeqIO.parse(gbk_path, "genbank"))
    cds_count = sum(1 for record in records for feature in record.features if feature.type == "CDS")
    return {
        "record_count": len(records),
        "cds_count": cds_count,
        "record_ids": [record.id for record in records],
    }

def _fixed_target_output_segments(calls):
    accepted_segments = {}
    for call in calls:
        best = call.best_hit
        if call.call == "accept" and best and best.segment and call.contig_id not in accepted_segments:
            accepted_segments[call.contig_id] = best.segment
    return accepted_segments

def _join_log(stdout_buf, stderr_buf):
    return chr(10).join(
        part for part in [stdout_buf.getvalue().strip(), stderr_buf.getvalue().strip()] if part
    )

def _summary_parse_attributes(attributes):
    if not attributes or attributes == ".":
        return {}
    parsed = {}
    for item in str(attributes).split(";"):
        if not item:
            continue
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        parsed[key] = value
    return parsed

def _summary_parse_target(value):
    parts = str(value or "").split()
    product = parts[0] if parts else ""
    start = int(parts[1]) if len(parts) > 1 else 1
    end = int(parts[2]) if len(parts) > 2 else start
    return product, start, end

def _summary_float(value, default=0.0):
    try:
        if value in (None, "", "."):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default

def _summary_int(value, default=0):
    try:
        if value in (None, "", "."):
            return default
        return int(value)
    except (TypeError, ValueError):
        return default

def _parse_paf_tags(fields):
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags

def _feature_metric_keys(product):
    product = str(product or "")
    keys = [product]
    if "_" in product:
        keys.append(product.split("_", 1)[0])
    return [key for key in dict.fromkeys(keys) if key]

def _parse_gff3_feature_metrics(gff3_text):
    paf_meta = {}
    metrics = defaultdict(list)
    for line in str(gff3_text or "").splitlines():
        if not line.strip():
            continue
        if line.startswith(f"##PAF{chr(9)}"):
            fields = line.split(chr(9))
            if len(fields) < 7:
                continue
            try:
                product = fields[1]
                target_length = int(fields[2])
                target_start = int(fields[3]) + 1
                target_end = int(fields[4])
                contig_id = fields[6]
            except (TypeError, ValueError, IndexError):
                continue
            tags = _parse_paf_tags(fields[13:])
            paf_meta[(contig_id, product, target_start, target_end)] = {
                "target_length": target_length,
                "raw_score": _summary_float(tags.get("AS")),
            }
            continue
        if line.startswith("#"):
            continue

        columns = line.split(chr(9))
        if len(columns) != 9:
            continue
        seqid, _, feature_type, start, end, score, strand, _, attributes = columns
        attrs = _summary_parse_attributes(attributes)
        if feature_type != "mRNA" or "Target" not in attrs:
            continue

        product, target_start, target_end = _summary_parse_target(attrs["Target"])
        paf = paf_meta.get((seqid, product, target_start, target_end), {})
        target_length = paf.get("target_length") or target_end
        aa_coverage = 0.0
        if target_length:
            aa_coverage = max(0, target_end - target_start + 1) / target_length
        identity = _summary_float(attrs.get("Identity"))
        normalized_score = identity * aa_coverage
        metric = {
            "product": product,
            "reference": {
                "product": product,
                "target_range": f"{target_start}-{target_end}/{target_length}" if target_length else "-",
            },
            "metrics": {
                "identity": identity,
                "aa_coverage": aa_coverage,
                "score": normalized_score,
            },
            "query_range": f"{_summary_int(start)}-{_summary_int(end)}",
            "target_range": f"{target_start}-{target_end}/{target_length}" if target_length else "-",
            "strand": strand or ".",
            "raw_score": _summary_float(score, paf.get("raw_score", 0.0)),
        }
        for key in _feature_metric_keys(product):
            metrics[(seqid, key)].append(metric)
    return metrics

def _first_qualifier(feature, key, default=""):
    value = feature.qualifiers.get(key, default)
    if isinstance(value, list):
        return value[0] if value else default
    if value is None:
        return default
    return value

def _qualifier_list(feature, key):
    value = feature.qualifiers.get(key, [])
    if value in (None, ""):
        return []
    if isinstance(value, list):
        return [str(item) for item in value if str(item)]
    return [str(value)]

def _strand_symbol(strand):
    if strand == 1:
        return "+"
    if strand == -1:
        return "-"
    return "."

def _feature_location_dict(feature):
    location = feature.location
    parts = list(getattr(location, "parts", [location]))
    ranges = []
    starts = []
    ends = []
    strands = []
    for part in parts:
        start = int(part.start) + 1
        end = int(part.end)
        starts.append(start)
        ends.append(end)
        strands.append(part.strand)
        range_text = f"{start}..{end}" if start != end else str(start)
        if part.strand == -1:
            range_text = f"complement({range_text})"
        ranges.append(range_text)
    strand = location.strand
    if strand is None and strands:
        unique_strands = {value for value in strands if value is not None}
        strand = unique_strands.pop() if len(unique_strands) == 1 else None
    return {
        "location": ranges[0] if len(ranges) == 1 else f"join({','.join(ranges)})",
        "start": min(starts) if starts else 0,
        "end": max(ends) if ends else 0,
        "strand": _strand_symbol(strand),
    }

def _feature_sequence_dict(record, feature):
    cds_nt = ""
    aa = ""
    try:
        cds_nt = str(feature.extract(record.seq))
    except Exception:
        cds_nt = ""
    translation = _first_qualifier(feature, "translation")
    if not translation and feature.type == "CDS":
        translated = gff3togbk.get_feature_translation(feature, record)
        translation = str(translated) if translated else ""
    aa = str(translation or "")
    return {"cds_nt": cds_nt, "aa": aa}

def _sanitize_summary_id(value):
    value = str(value or "feature")
    cleaned = "".join(char if char.isalnum() or char in "._-" else "_" for char in value)
    return cleaned.strip("._") or "feature"

def _best_metric_for_feature(metrics_by_key, input_id, feature):
    gene = str(_first_qualifier(feature, "gene") or "")
    product = str(_first_qualifier(feature, "product") or "")
    candidates = []
    for key in dict.fromkeys([gene, product]):
        if key:
            candidates.extend(metrics_by_key.get((input_id, key), []))
    if not candidates:
        return None
    return max(
        candidates,
        key=lambda item: (
            item.get("metrics", {}).get("score", 0.0),
            item.get("metrics", {}).get("aa_coverage", 0.0),
            item.get("metrics", {}).get("identity", 0.0),
        ),
    )

def _candidate_hit_dict(hit):
    if not hit:
        return None
    return {
        "product": hit.product,
        "reference": hit.product,
        "identity": hit.identity,
        "aa_coverage": hit.aa_coverage,
        "score": hit.normalized_score,
        "query_range": hit.query_range,
        "target_range": hit.target_range,
        "strand": hit.strand,
    }

def _build_genbank_feature_data(gbk_path, calls, filtered_gff3_text, target, accepted_input_ids=None):
    records = list(SeqIO.parse(gbk_path, "genbank"))
    calls_by_id = {call.contig_id: call for call in calls}
    if accepted_input_ids is None:
        accepted_input_ids = [
            call.contig_id
            for call in calls
            if call.call == "accept" and call.target == target
        ]
    metrics_by_key = _parse_gff3_feature_metrics(filtered_gff3_text)
    feature_data = {}
    for input_id, record in zip(accepted_input_ids, records):
        call = calls_by_id.get(input_id)
        feature_summaries = []
        seen = defaultdict(int)
        for feature in record.features:
            if feature.type not in {"CDS", "misc_feature"}:
                continue
            gene = str(_first_qualifier(feature, "gene") or "")
            product = str(_first_qualifier(feature, "product") or gene or feature.type)
            if not gene and not product:
                continue
            location = _feature_location_dict(feature)
            metric = _best_metric_for_feature(metrics_by_key, input_id, feature)
            base_id = _sanitize_summary_id(f"{record.id}_{gene or product or feature.type}")
            seen[base_id] += 1
            feature_id = base_id if seen[base_id] == 1 else f"{base_id}_{seen[base_id]}"
            notes = _qualifier_list(feature, "note")
            if "ribosomal_slippage" in feature.qualifiers:
                notes.append("ribosomal slippage")
            feature_summaries.append(
                {
                    "id": feature_id,
                    "type": feature.type,
                    "gene": gene,
                    "product": product,
                    "segment": call.segment if call and call.segment != "-" else "",
                    "location": location["location"],
                    "start": location["start"],
                    "end": location["end"],
                    "strand": location["strand"],
                    "reference": metric.get("reference", {}) if metric else {},
                    "metrics": metric.get("metrics", {}) if metric else {},
                    "query_range": metric.get("query_range", "") if metric else "",
                    "target_range": metric.get("target_range", "") if metric else "",
                    "sequences": _feature_sequence_dict(record, feature),
                    "notes": notes,
                }
            )
        feature_data[input_id] = {
            "record_id": record.id,
            "length": len(record.seq),
            "features": feature_summaries,
        }
    return feature_data

def _build_contig_summaries(calls, feature_data_by_input):
    contigs = []
    for call in calls:
        feature_data = feature_data_by_input.get(call.contig_id, {})
        contigs.append(
            {
                "input_id": call.contig_id,
                "record_id": feature_data.get("record_id", ""),
                "target": call.target if call.target != "-" else "",
                "segment": call.segment if call.segment != "-" else "",
                "length": feature_data.get("length", call.length),
                "call": call.call,
                "status": call.status,
                "qc_result": call.qc_result,
                "confidence": call.confidence,
                "best_hit": _candidate_hit_dict(call.best_hit),
                "features": feature_data.get("features", []),
                "flags": list(call.flags),
                "notes": list(call.notes),
            }
        )
    return contigs

def _download_kind(name):
    lowered = str(name).lower()
    if lowered.endswith(".summary.json"):
        return "summary_json"
    if lowered.endswith(".gff3"):
        return "gff3"
    if lowered.endswith(".gbk"):
        return "genbank"
    if lowered.endswith(".cds.fna"):
        return "cds_fasta"
    if lowered.endswith(".faa"):
        return "protein_fasta"
    if lowered.endswith(".hits.tsv") or lowered.endswith(".auto.tsv"):
        return "contig_report"
    if lowered.endswith(".log"):
        return "log"
    return "output"

def _download_target(name):
    for target in AUTO_TARGETS:
        if f".{target}." in str(name):
            return target
    return ""

def _build_download_manifest(outputs):
    order = {
        "summary_json": 0,
        "gff3": 1,
        "genbank": 2,
        "cds_fasta": 3,
        "protein_fasta": 4,
        "contig_report": 5,
        "log": 6,
        "output": 7,
    }
    downloads = []
    for name in outputs:
        kind = _download_kind(name)
        downloads.append(
            {
                "name": name,
                "kind": kind,
                "target": _download_target(name),
            }
        )
    return sorted(downloads, key=lambda item: (order.get(item["kind"], 99), item["name"]))

def _run_counts_for_summary(calls, annotated_contigs):
    call_counts = Counter(call.call for call in calls)
    qc_counts = Counter(call.qc_result for call in calls)
    return {
        "input_contigs": len(calls),
        "annotated_contigs": annotated_contigs,
        "accepted": call_counts.get("accept", 0),
        "rejected": call_counts.get("reject", 0),
        "passed": qc_counts.get("pass", 0),
        "failed": qc_counts.get("fail", 0),
    }

def run_ganflu_web(input_fasta, gff3_text, target, isolate, output_stem="ganflu", preserve_original_id=False, hit_settings_json=None):
    target = str(target or "").upper()
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    original_streams = []
    try:
        if target not in SUPPORTED_TARGETS:
            raise ValueError(f"Unsupported target: {target}")
        if not str(input_fasta or "").strip():
            raise ValueError("Input FASTA is empty.")
        if not str(gff3_text or "").strip():
            raise ValueError("Miniprot GFF3 output is empty.")
        isolate_value = str(isolate or "").strip() or _safe_stem(output_stem)

        work_dir = "/tmp/ganflu-web"
        os.makedirs(work_dir, exist_ok=True)
        for old_path in glob.glob(os.path.join(work_dir, "*")):
            try:
                os.remove(old_path)
            except OSError:
                pass

        stem = _safe_stem(output_stem)
        input_path = os.path.join(work_dir, f"{stem}.fasta")
        raw_input_path = os.path.join(work_dir, f"{stem}.input.fasta")
        scan_gff3_path = os.path.join(work_dir, f"{stem}.scan.gff3")
        gff3_path = os.path.join(work_dir, f"{stem}.gff3")
        toml_path = os.path.join(work_dir, f"{target}.toml")
        gbk_path = os.path.join(work_dir, f"{stem}.gbk")
        cds_path = os.path.join(work_dir, f"{stem}.cds.fna")
        faa_path = os.path.join(work_dir, f"{stem}.faa")

        target_root = resources.files("ganflu").joinpath("db", target)
        with open(raw_input_path, "w", encoding="utf-8") as handle:
            handle.write(str(input_fasta))
        with open(scan_gff3_path, "w", encoding="utf-8") as handle:
            handle.write(str(gff3_text))
        with open(toml_path, "w", encoding="utf-8") as handle:
            handle.write(_read_resource_text(target_root.joinpath(f"{target}.toml")))

        contigs = list(SeqIO.parse(io.StringIO(str(input_fasta)), "fasta"))
        if not contigs:
            raise ValueError("Input FASTA contains no records.")
        contigs_by_id = {record.id: record for record in contigs}
        if len(contigs_by_id) != len(contigs):
            raise ValueError("Input FASTA contains duplicate record IDs.")

        reference = _load_web_reference_bundle(target, work_dir)
        hit_settings = _normalize_hit_settings(hit_settings_json)
        thresholds = _thresholds_from_hit_settings(hit_settings)
        candidates = auto_mode.parse_miniprot_gff3(
            scan_gff3_path,
            reference,
            contigs_by_id,
            thresholds,
        )
        candidates_by_contig = defaultdict(list)
        for candidate in candidates:
            candidates_by_contig[candidate.contig_id].append(candidate)
        calls = auto_mode.classify_contigs(contigs, candidates_by_contig, thresholds)
        for call in calls:
            best = call.best_hit
            if best and best.segment:
                call.target = target
                call.segment = best.segment
        accepted_segments = _fixed_target_output_segments(calls)
        if not accepted_segments:
            raise ValueError(f"No segment-compatible miniprot hits were found for target {target}.")

        auto_mode.write_target_fasta(contigs, accepted_segments, input_path)
        auto_mode.filter_gff3_for_target(
            scan_gff3_path,
            gff3_path,
            accepted_segments,
            reference,
        )
        hit_tsv_path = os.path.join(work_dir, f"{stem}.hits.tsv")
        auto_mode.write_auto_tsv(calls, hit_tsv_path)

        root_logger = logging.getLogger()
        for handler in root_logger.handlers:
            if isinstance(handler, logging.StreamHandler):
                original_streams.append((handler, handler.stream))
                handler.setStream(stdout_buf)

        args = [
            "-i", input_path,
            "-g", gff3_path,
            "--toml", toml_path,
            "-o", gbk_path,
            "--isolate", isolate_value,
            "--cds-fna", cds_path,
            "--faa", faa_path,
        ]
        if bool(preserve_original_id):
            args.append("--preserve_original_id")

        with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
            gff3togbk.main(args)

        with open(gbk_path, "r", encoding="utf-8") as handle:
            gbk_text = handle.read()
        with open(cds_path, "r", encoding="utf-8") as handle:
            cds_text = handle.read()
        with open(faa_path, "r", encoding="utf-8") as handle:
            faa_text = handle.read()
        filtered_gff3_text = _read_text_file(gff3_path)
        hit_tsv_text = _read_text_file(hit_tsv_path)
        log_text = _join_log(stdout_buf, stderr_buf)

        output_names = {
            f"{stem}.gff3": f"{stem}.gff3",
            f"{stem}.gbk": f"{stem}.gbk",
            f"{stem}.cds.fna": f"{stem}.cds.fna",
            f"{stem}.faa": f"{stem}.faa",
            f"{stem}.hits.tsv": f"{stem}.hits.tsv",
            f"{stem}.summary.json": f"{stem}.summary.json",
            f"{stem}.log": f"{stem}.log",
        }
        call_counts = Counter(call.call for call in calls)
        status_counts = Counter(call.status for call in calls)
        qc_counts = Counter(call.qc_result for call in calls)
        segment_counts = Counter(accepted_segments.values())
        feature_data = _build_genbank_feature_data(
            gbk_path,
            calls,
            filtered_gff3_text,
            target,
            accepted_input_ids=list(accepted_segments),
        )
        contig_summaries = _build_contig_summaries(calls, feature_data)
        run_counts = _run_counts_for_summary(calls, len(accepted_segments))
        summary = {
            **_count_genbank_outputs(gbk_path),
            "schema_version": 1,
            "mode": "target",
            "target": target,
            "isolate": isolate_value,
            "output_stem": stem,
            "input_contigs": len(contigs),
            "annotated_contigs": len(accepted_segments),
            "skipped_contigs": len(contigs) - len(accepted_segments),
            "gff3_lines": len([line for line in filtered_gff3_text.splitlines() if line.strip()]),
            "by_call": dict(call_counts),
            "by_status": dict(status_counts),
            "by_qc_result": dict(qc_counts),
            "by_segment": dict(segment_counts),
            "passed_contigs": qc_counts.get("pass", 0),
            "failed_contigs": qc_counts.get("fail", 0),
            "thresholds": thresholds.as_dict(),
            "miniprot": _miniprot_settings_for_summary(hit_settings),
            "outputs": output_names,
            "run": {
                "target": target,
                "targets_scanned": [target],
                "isolate": isolate_value,
                "output_stem": stem,
                "counts": run_counts,
            },
            "contigs": contig_summaries,
            "downloads": _build_download_manifest(output_names),
            "run_log": {
                "thresholds": thresholds.as_dict(),
                "miniprot": _miniprot_settings_for_summary(hit_settings),
                "messages": log_text,
            },
        }
        summary_text = json.dumps(summary, indent=2, sort_keys=True) + chr(10)
        return json.dumps(
            {
                "summary": summary,
                "log": log_text,
                "outputs": {
                    f"{stem}.gff3": filtered_gff3_text,
                    f"{stem}.gbk": gbk_text,
                    f"{stem}.cds.fna": cds_text,
                    f"{stem}.faa": faa_text,
                    f"{stem}.hits.tsv": hit_tsv_text,
                    f"{stem}.summary.json": summary_text,
                    f"{stem}.log": log_text,
                },
            }
        )
    except Exception as exc:
        return json.dumps(
            {
                "error": {
                    "type": exc.__class__.__name__,
                    "message": str(exc) or "ganflu failed",
                    "traceback": traceback.format_exc(),
                    "log": _join_log(stdout_buf, stderr_buf),
                }
            }
        )
    finally:
        for handler, stream in original_streams:
            handler.setStream(stream)

def _load_web_reference_bundle(target, work_dir):
    target = str(target or "").upper()
    if target not in SUPPORTED_TARGETS:
        raise ValueError(f"Unsupported target: {target}")
    target_root = resources.files("ganflu").joinpath("db", target)
    toml_text = _read_resource_text(target_root.joinpath(f"{target}.toml"))
    config = _load_toml_text(toml_text)
    prot_text = _read_resource_text(_resource_join(target_root, config["metadata"]["prot_faa"]))
    protein_lengths = {
        record.id: len(record.seq)
        for record in SeqIO.parse(io.StringIO(prot_text), "fasta")
    }
    toml_path = os.path.join(work_dir, f"{target}.toml")
    with open(toml_path, "w", encoding="utf-8") as handle:
        handle.write(toml_text)
    return auto_mode.ReferenceBundle(
        target=target,
        ref_dir=work_dir,
        toml_path=toml_path,
        prot_faa="",
        config=config,
        segment_keys=list(config.get("segments", {}).keys()),
        gene_configs=config.get("genes", {}),
        protein_lengths=protein_lengths,
    )

def _read_text_file(path):
    with open(path, "r", encoding="utf-8") as handle:
        return handle.read()

def run_ganflu_auto_web(input_fasta, gff3_by_target_json, isolate, output_stem="ganflu", preserve_original_id=False, hit_settings_json=None):
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    original_streams = []
    try:
        if not str(input_fasta or "").strip():
            raise ValueError("Input FASTA is empty.")
        isolate_value = str(isolate or "").strip() or _safe_stem(output_stem)

        if isinstance(gff3_by_target_json, str):
            gff3_by_target = json.loads(gff3_by_target_json)
        else:
            gff3_by_target = dict(gff3_by_target_json or {})

        targets = [
            target
            for target in AUTO_TARGETS
            if str(gff3_by_target.get(target, "")).strip()
        ]
        if not targets:
            raise ValueError("No auto GFF3 payloads were provided.")

        work_dir = "/tmp/ganflu-web"
        os.makedirs(work_dir, exist_ok=True)
        for old_path in glob.glob(os.path.join(work_dir, "*")):
            try:
                os.remove(old_path)
            except OSError:
                pass

        stem = _safe_stem(output_stem)
        input_path = os.path.join(work_dir, f"{stem}.fasta")
        with open(input_path, "w", encoding="utf-8") as handle:
            handle.write(str(input_fasta))

        contigs = list(SeqIO.parse(io.StringIO(str(input_fasta)), "fasta"))
        if not contigs:
            raise ValueError("Input FASTA contains no records.")
        contigs_by_id = {record.id: record for record in contigs}
        if len(contigs_by_id) != len(contigs):
            raise ValueError("Input FASTA contains duplicate record IDs.")

        root_logger = logging.getLogger()
        for handler in root_logger.handlers:
            if isinstance(handler, logging.StreamHandler):
                original_streams.append((handler, handler.stream))
                handler.setStream(stdout_buf)

        hit_settings = _normalize_hit_settings(hit_settings_json)
        thresholds = _thresholds_from_hit_settings(hit_settings)
        references = {}
        scan_gff3_by_target = {}
        candidates_by_contig = defaultdict(list)

        for target in targets:
            references[target] = _load_web_reference_bundle(target, work_dir)
            scan_gff3 = os.path.join(work_dir, f"{stem}.{target}.scan.gff3")
            scan_gff3_by_target[target] = scan_gff3
            with open(scan_gff3, "w", encoding="utf-8") as handle:
                handle.write(str(gff3_by_target[target]))
            candidates = auto_mode.parse_miniprot_gff3(
                scan_gff3,
                references[target],
                contigs_by_id,
                thresholds,
            )
            for candidate in candidates:
                candidates_by_contig[candidate.contig_id].append(candidate)

        calls = auto_mode.classify_contigs(contigs, candidates_by_contig, thresholds)
        accepted_by_target = auto_mode.make_accepted_segments(calls)
        outputs = {}
        feature_data_by_input = {}

        for target in sorted(accepted_by_target):
            accepted_segments = accepted_by_target[target]
            if not accepted_segments:
                continue
            target_fasta = os.path.join(work_dir, f"{stem}.{target}.accepted.fasta")
            target_gff3 = os.path.join(work_dir, f"{stem}.{target}.gff3")
            target_gbk = os.path.join(work_dir, f"{stem}.{target}.gbk")
            target_cds = os.path.join(work_dir, f"{stem}.{target}.cds.fna")
            target_faa = os.path.join(work_dir, f"{stem}.{target}.faa")

            auto_mode.write_target_fasta(contigs, accepted_segments, target_fasta)
            auto_mode.filter_gff3_for_target(
                scan_gff3_by_target[target],
                target_gff3,
                accepted_segments,
                references[target],
            )

            args = [
                "-i", target_fasta,
                "-g", target_gff3,
                "--toml", references[target].toml_path,
                "-o", target_gbk,
                "--isolate", isolate_value,
                "--cds-fna", target_cds,
                "--faa", target_faa,
            ]
            if bool(preserve_original_id):
                args.append("--preserve_original_id")

            with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
                gff3togbk.main(args)

            target_gff3_text = _read_text_file(target_gff3)
            outputs[f"{stem}.{target}.gff3"] = target_gff3_text
            outputs[f"{stem}.{target}.gbk"] = _read_text_file(target_gbk)
            outputs[f"{stem}.{target}.cds.fna"] = _read_text_file(target_cds)
            outputs[f"{stem}.{target}.faa"] = _read_text_file(target_faa)
            feature_data_by_input.update(
                _build_genbank_feature_data(
                    target_gbk,
                    calls,
                    target_gff3_text,
                    target,
                    accepted_input_ids=list(accepted_segments),
                )
            )

        auto_tsv_path = os.path.join(work_dir, f"{stem}.auto.tsv")
        auto_mode.write_auto_tsv(calls, auto_tsv_path)
        outputs[f"{stem}.auto.tsv"] = _read_text_file(auto_tsv_path)
        log_text = _join_log(stdout_buf, stderr_buf)
        outputs[f"{stem}.auto.log"] = log_text

        summary_outputs = {name: name for name in outputs}
        summary_outputs[f"{stem}.auto.summary.json"] = f"{stem}.auto.summary.json"
        summary = auto_mode.build_summary(
            input_fasta=input_path,
            output_stem=stem,
            targets=targets,
            thresholds=thresholds,
            calls=calls,
            outputs=summary_outputs,
        )
        summary["auto"] = True
        summary["schema_version"] = 1
        summary["mode"] = "auto"
        summary["miniprot"] = _miniprot_settings_for_summary(hit_settings)
        summary["run"] = {
            "target": "auto",
            "targets_scanned": targets,
            "isolate": isolate_value,
            "output_stem": stem,
            "counts": _run_counts_for_summary(calls, summary["counts"].get("accepted", 0)),
        }
        summary["contigs"] = _build_contig_summaries(calls, feature_data_by_input)
        summary["downloads"] = _build_download_manifest(summary_outputs)
        summary["run_log"] = {
            "thresholds": thresholds.as_dict(),
            "miniprot": _miniprot_settings_for_summary(hit_settings),
            "messages": log_text,
        }
        summary_text = json.dumps(summary, indent=2, sort_keys=True) + chr(10)
        outputs[f"{stem}.auto.summary.json"] = summary_text

        return json.dumps(
            {
                "summary": summary,
                "log": log_text,
                "outputs": outputs,
            }
        )
    except Exception as exc:
        return json.dumps(
            {
                "error": {
                    "type": exc.__class__.__name__,
                    "message": str(exc) or "ganflu auto mode failed",
                    "traceback": traceback.format_exc(),
                    "log": _join_log(stdout_buf, stderr_buf),
                }
            }
        )
    finally:
        for handler, stream in original_streams:
            handler.setStream(stream)
`;
