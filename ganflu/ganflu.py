#!/usr/bin/env python
# coding: utf-8

import os
import sys
import logging
import argparse
import time
import functools
import http.server
import webbrowser
import tempfile
from importlib import resources
from . import __version__
from .launchers.miniprot import MiniprotCommandLine
from .scripts import auto_mode, gff3_prune, gff3togbk, validate_reference_files

SUPPORTED_TARGETS = ["IAV", "IBV", "ICV", "IDV"]
CLI_TARGETS = SUPPORTED_TARGETS + ["auto"]
GUI_COMMAND = "gui"

def _version():
    """
    Returns the version of the pipeline
    """
    return __version__

def is_positive_integer(parameter_name, value):
    try:
        value = int(value)
        if value <= 0:
            raise argparse.ArgumentTypeError(f"{parameter_name} must be a positive integer")
    except ValueError:
        raise Exception(f"{parameter_name} must be a positive integer")
    return value

def setup_logging(log_file, verbose=False):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    log_format = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    console_handler.setFormatter(log_format)
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)

    return logger

def resolve_isolate(isolate, output_stem):
    return (isolate or "").strip() or gff3togbk.get_output_id_prefix(output_stem)

def get_webapp_dir():
    web_dir = resources.files("ganflu").joinpath("web")
    index_path = web_dir.joinpath("index.html")
    if not index_path.is_file():
        raise FileNotFoundError(
            "ganflu web app assets are not installed. Reinstall ganflu from a package "
            "that includes web assets, or run from a source checkout."
        )
    return web_dir

def _get_gui_args(raw_args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="ganflu gui",
        description=f"ganflu v{_version()}: launch the browser-based ganflu Web app",
    )
    parser.add_argument("--host", default="127.0.0.1", help="Host/interface to bind (default: 127.0.0.1)")
    parser.add_argument("--port", default=0, type=int, help="Port to bind (default: 0, auto-select a free port)")
    parser.add_argument("--no-port-fallback", dest="port_fallback", action="store_false", help="Fail instead of trying the next port when an explicit --port is busy")
    parser.add_argument("--open-browser", action="store_true", help="Open the web app URL in the default browser")
    parser.add_argument("-v", "--version", action="version", version=_version())
    args = parser.parse_args(raw_args)
    args.command = GUI_COMMAND
    return args

def _get_args() -> argparse.Namespace:
    if len(sys.argv) > 1 and sys.argv[1] == GUI_COMMAND:
        return _get_gui_args(sys.argv[2:])

    parser = argparse.ArgumentParser(
        description=f"ganflu v{_version()}: Influenza virus genome annotation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Subcommands:\n"
            "  ganflu gui              Launch the browser-based ganflu Web app\n"
            "  ganflu gui --help       Show ganflu Web app options\n"
        ),
    )
    
    # Input/Output options
    parser.add_argument("-i", "--input", required=True, type=str, help="Input FASTA file")
    parser.add_argument("-o", "--output", dest="output", default=None, type=str, help="basename for Output GenBank file name (default: <input>)")
    parser.add_argument("-t", "--target", dest="target", required=True, help="Target virus", choices=CLI_TARGETS)
    parser.add_argument("-d", "--db_dir", dest="db_dir", help="Data path (optional; default: ganflu/db)", default=None)
    parser.add_argument("--isolate", dest="isolate", default=None, help='isolate name (default: output stem prefix; e.g. "A/Narita/1/2009", "A/goose/Guangdong/1/1996", "B/Lee/1940", "C/Ann_Arbor/1/1950", "D/swine/Oklahoma/1334/2011")')
    parser.add_argument("--preserve_original_id", "--preserve-original-id", dest="preserve_original_id", action="store_true", help="Preserve original FASTA record IDs in GenBank output")
    parser.add_argument("--log-file", dest="log_file", default=None, help="Log file path (default: <output>.log; auto mode: <output>.auto.log)")
    parser.add_argument("--verbose", action="store_true", help="Show debug logs in the terminal")
    parser.add_argument("--auto-targets", dest="auto_targets", default="IAV,IBV,ICV,IDV", help="Comma-separated targets to scan in auto mode (default: IAV,IBV,ICV,IDV)")
    parser.add_argument("--auto-min-identity", dest="auto_min_identity", default=0.55, type=float, help="Minimum amino-acid identity for auto candidate hits")
    parser.add_argument("--auto-min-aa-coverage", dest="auto_min_aa_coverage", default=0.35, type=float, help="Minimum reference amino-acid coverage for auto candidate hits")
    parser.add_argument("--auto-min-score", dest="auto_min_score", default=0.25, type=float, help="Minimum normalized auto score (identity * coverage)")
    parser.add_argument("--auto-min-margin", dest="auto_min_margin", default=0.10, type=float, help="Minimum score margin between best and second-best target in auto mode")
    parser.add_argument("--auto-complete-aa-coverage", dest="auto_complete_aa_coverage", default=0.90, type=float, help="Reference amino-acid coverage required to call an auto hit complete")
    parser.add_argument("--auto-write-rejected", dest="auto_write_rejected", action="store_true", help="Write rejected/review contigs to <output>.auto.rejected.fasta")
    parser.add_argument("--auto-report-prefix", dest="auto_report_prefix", default=None, help="Output prefix for auto TSV/summary reports (default: <output>)")
    parser.add_argument("-v", "--version", action="version", version=_version())
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()
    args.command = "annotate"

    return args 

def bind_gui_server(host, port, web_dir, port_fallback=True):
    handler = functools.partial(
        http.server.SimpleHTTPRequestHandler,
        directory=str(web_dir),
    )
    current_port = int(port)
    last_error = None
    for _ in range(100 if port_fallback else 1):
        try:
            return http.server.ThreadingHTTPServer((host, current_port), handler)
        except OSError as exc:
            last_error = exc
            if exc.errno not in {48, 98, 10048}:
                break
            current_port += 1
    raise OSError(f"Could not bind ganflu gui server starting at {host}:{port}: {last_error}")

def get_server_url(server):
    host, port = server.server_address[:2]
    display_host = "127.0.0.1" if host in {"0.0.0.0", ""} else host
    if ":" in display_host and not display_host.startswith("["):
        display_host = f"[{display_host}]"
    return f"http://{display_host}:{port}/"

def run_gui(args) -> int:
    web_dir = get_webapp_dir()
    server = bind_gui_server(
        args.host,
        args.port,
        web_dir,
        port_fallback=args.port_fallback,
    )
    url = get_server_url(server)
    print(f"Serving ganflu Web from {web_dir}")
    print(f"Open {url}")
    print("Press Ctrl+C to stop.")
    if args.open_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping ganflu Web.")
    finally:
        server.server_close()
    return 0

def main():
    start_time = time.time()
    args = _get_args()

    if getattr(args, "command", None) == GUI_COMMAND:
        return run_gui(args)

    input_fasta = args.input
    if args.output:
        out_stem = os.path.abspath(args.output)
    else:
        input_dir = os.path.dirname(os.path.abspath(input_fasta))
        input_stem = os.path.splitext(os.path.basename(input_fasta))[0]
        out_stem = os.path.join(input_dir, input_stem)
    # workdir is the directory where the output files will be saved
    work_dir = os.path.dirname(out_stem)
    os.makedirs(work_dir, exist_ok=True)
    args.isolate = resolve_isolate(args.isolate, out_stem)

    default_log_file = f"{out_stem}.auto.log" if args.target == "auto" else f"{out_stem}.log"
    log_file = os.path.abspath(args.log_file) if args.log_file else default_log_file
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    logger = setup_logging(log_file, args.verbose)
    logger.info(f"ganflu v{_version()} started")
    logger.info(f"Log file: {log_file}")
    logger.debug(f"Command line: {' '.join(sys.argv)}")
    logger.info(f"Input FASTA: {os.path.abspath(input_fasta)}")
    logger.info(f"Output stem: {out_stem}")
    logger.info(f"Work directory: {work_dir}")
    logger.info(f"Target: {args.target}")

    target = args.target
    if target == "auto":
        try:
            auto_mode.run_auto(args, out_stem, work_dir, logger)
            logger.info(f"ganflu auto mode completed in {time.time() - start_time:.2f} seconds")
        except Exception:
            logger.exception(f"ganflu auto mode failed after {time.time() - start_time:.2f} seconds")
            raise
        return 0

    if args.db_dir:
        ref_dir = os.path.abspath(args.db_dir)
    else:
        db_path_traversable = resources.files('ganflu').joinpath(f'db/{args.target}')
        ref_dir = db_path_traversable.resolve()

    logger.info(f"Reference directory: {ref_dir}")
    ref_toml = validate_reference_files.validate_reference_files(target=args.target, db_dir=args.db_dir, logger=logger)
    ref_faa = os.path.join(ref_dir, ref_toml["metadata"]["prot_faa"])
    logger.info(f"Reference protein FASTA: {ref_faa}")

    gff3_file = f"{out_stem}.gff3"
    gbk_file = f"{out_stem}.gbk"
    cds_fna_file = f"{out_stem}.cds.fna"
    faa_file = f"{out_stem}.faa"
    raw_gff3_file = None
    try:
        logger.info("Running miniprot")
        with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".raw.gff3",
            prefix=f"{os.path.basename(out_stem)}.",
            dir=work_dir,
            delete=False,
        ) as raw_gff3:
            raw_gff3_file = raw_gff3.name
        miniprot = MiniprotCommandLine(
            input=input_fasta, work_dir=work_dir, output=raw_gff3_file,
            prot_faa=ref_faa, miniprot_bin="miniprot", stderr_filename="miniprot.stderr", kmer_size=15,
            max_secondary_alignments=gff3_prune.RELAXED_MAX_SECONDARY_ALIGNMENTS,
            secondary_to_primary_ratio=gff3_prune.RELAXED_SECONDARY_TO_PRIMARY_RATIO,
            output_score_ratio=gff3_prune.RELAXED_OUTPUT_SCORE_RATIO,
            )
        miniprot.run_piped_commands()
        logger.info("Pruning miniprot GFF3")
        prune_result = gff3_prune.prune_gff3(
            raw_gff3_file,
            gff3_file,
            prot_faa=ref_faa,
            antigen_names=ref_toml.get("serotype", {}).keys(),
        )
        logger.info(
            f"Pruned GFF3 output: {gff3_file} "
            f"({prune_result.selected_parent_count} parent alignment(s))"
        )
        try:
            os.remove(raw_gff3_file)
            raw_gff3_file = None
        except OSError:
            logger.debug(f"Could not remove temporary raw GFF3: {raw_gff3_file}", exc_info=True)

        logger.info("Converting GFF3 to GenBank")
        gff3togbk_args = [
            "-g", gff3_file,
            "-o", gbk_file,
            "-i", input_fasta,
            "--toml", os.path.join(ref_dir, f'{target}.toml'),
            "--isolate", args.isolate,
            "--cds-fna", cds_fna_file,
            "--faa", faa_file,
        ]
        if args.preserve_original_id:
            gff3togbk_args.append("--preserve_original_id")
        gff3togbk.main(gff3togbk_args)
        logger.info(f"GenBank output: {gbk_file}")
        logger.info(f"ganflu completed in {time.time() - start_time:.2f} seconds")
    except Exception:
        logger.exception(f"ganflu failed after {time.time() - start_time:.2f} seconds")
        raise

    return 0

if __name__ == "__main__":
    main()
