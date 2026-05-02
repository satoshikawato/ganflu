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
        if best and best.segment and call.contig_id not in accepted_segments:
            accepted_segments[call.contig_id] = best.segment
    return accepted_segments

def run_ganflu_web(input_fasta, gff3_text, target, isolate, output_stem="ganflu", preserve_original_id=False):
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
        thresholds = auto_mode.AutoThresholds()
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

        output_names = {
            f"{stem}.gff3": f"{stem}.gff3",
            f"{stem}.gbk": f"{stem}.gbk",
            f"{stem}.cds.fna": f"{stem}.cds.fna",
            f"{stem}.faa": f"{stem}.faa",
            f"{stem}.hits.tsv": f"{stem}.hits.tsv",
            f"{stem}.summary.json": f"{stem}.summary.json",
        }
        call_counts = Counter(call.call for call in calls)
        status_counts = Counter(call.status for call in calls)
        segment_counts = Counter(accepted_segments.values())
        summary = {
            **_count_genbank_outputs(gbk_path),
            "target": target,
            "isolate": isolate_value,
            "output_stem": stem,
            "input_contigs": len(contigs),
            "annotated_contigs": len(accepted_segments),
            "skipped_contigs": len(contigs) - len(accepted_segments),
            "gff3_lines": len([line for line in filtered_gff3_text.splitlines() if line.strip()]),
            "by_call": dict(call_counts),
            "by_status": dict(status_counts),
            "by_segment": dict(segment_counts),
            "outputs": output_names,
        }
        summary_text = json.dumps(summary, indent=2, sort_keys=True) + "\\n"
        return json.dumps(
            {
                "summary": summary,
                "log": "\\n".join(part for part in [stdout_buf.getvalue().strip(), stderr_buf.getvalue().strip()] if part),
                "outputs": {
                    f"{stem}.gff3": filtered_gff3_text,
                    f"{stem}.gbk": gbk_text,
                    f"{stem}.cds.fna": cds_text,
                    f"{stem}.faa": faa_text,
                    f"{stem}.hits.tsv": hit_tsv_text,
                    f"{stem}.summary.json": summary_text,
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
                    "log": "\\n".join(part for part in [stdout_buf.getvalue().strip(), stderr_buf.getvalue().strip()] if part),
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

def run_ganflu_auto_web(input_fasta, gff3_by_target_json, isolate, output_stem="ganflu", preserve_original_id=False):
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

        thresholds = auto_mode.AutoThresholds()
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

            outputs[f"{stem}.{target}.gff3"] = _read_text_file(target_gff3)
            outputs[f"{stem}.{target}.gbk"] = _read_text_file(target_gbk)
            outputs[f"{stem}.{target}.cds.fna"] = _read_text_file(target_cds)
            outputs[f"{stem}.{target}.faa"] = _read_text_file(target_faa)

        auto_tsv_path = os.path.join(work_dir, f"{stem}.auto.tsv")
        auto_mode.write_auto_tsv(calls, auto_tsv_path)
        outputs[f"{stem}.auto.tsv"] = _read_text_file(auto_tsv_path)

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
        summary_text = json.dumps(summary, indent=2, sort_keys=True) + "\\n"
        outputs[f"{stem}.auto.summary.json"] = summary_text

        return json.dumps(
            {
                "summary": summary,
                "log": "\\n".join(part for part in [stdout_buf.getvalue().strip(), stderr_buf.getvalue().strip()] if part),
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
                    "log": "\\n".join(part for part in [stdout_buf.getvalue().strip(), stderr_buf.getvalue().strip()] if part),
                }
            }
        )
    finally:
        for handler, stream in original_streams:
            handler.setStream(stream)
`;
