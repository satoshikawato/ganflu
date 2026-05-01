export const PYTHON_HELPERS = `
import contextlib
import glob
import io
import json
import logging
import os
import re
import traceback
from importlib import resources

from Bio import SeqIO
from ganflu.scripts import gff3togbk

SUPPORTED_TARGETS = {"IAV", "IBV", "ICV", "IDV"}

def _safe_stem(value):
    value = str(value or "ganflu").strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    return value or "ganflu"

def _read_resource_text(path):
    return path.read_text(encoding="utf-8")

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
        if not str(isolate or "").strip():
            raise ValueError("Isolate is required.")

        work_dir = "/tmp/ganflu-web"
        os.makedirs(work_dir, exist_ok=True)
        for old_path in glob.glob(os.path.join(work_dir, "*")):
            try:
                os.remove(old_path)
            except OSError:
                pass

        stem = _safe_stem(output_stem)
        input_path = os.path.join(work_dir, f"{stem}.fasta")
        gff3_path = os.path.join(work_dir, f"{stem}.gff3")
        toml_path = os.path.join(work_dir, f"{target}.toml")
        gbk_path = os.path.join(work_dir, f"{stem}.gbk")
        cds_path = os.path.join(work_dir, f"{stem}.cds.fna")
        faa_path = os.path.join(work_dir, f"{stem}.faa")

        target_root = resources.files("ganflu").joinpath("db", target)
        with open(input_path, "w", encoding="utf-8") as handle:
            handle.write(str(input_fasta))
        with open(gff3_path, "w", encoding="utf-8") as handle:
            handle.write(str(gff3_text))
        with open(toml_path, "w", encoding="utf-8") as handle:
            handle.write(_read_resource_text(target_root.joinpath(f"{target}.toml")))

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
            "--isolate", str(isolate).strip(),
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

        summary = _count_genbank_outputs(gbk_path)
        return json.dumps(
            {
                "summary": summary,
                "log": "\\n".join(part for part in [stdout_buf.getvalue().strip(), stderr_buf.getvalue().strip()] if part),
                "outputs": {
                    f"{stem}.gbk": gbk_text,
                    f"{stem}.cds.fna": cds_text,
                    f"{stem}.faa": faa_text,
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
`;
