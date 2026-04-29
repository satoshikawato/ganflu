#!/usr/bin/env python
# coding: utf-8

import os
import sys
import logging
import argparse
import time
from importlib import resources
from .launchers.miniprot import MiniprotCommandLine
from .scripts import gff3togbk, validate_reference_files

def _version():
    """
    Returns the version of the pipeline
    """
    return "0.1.0"

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

def _get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"ganflu v{_version()}: Influenza virus genome annotation")
    
    # Input/Output options
    parser.add_argument("-i", "--input", required=True, type=str, help="Input FASTA file")
    parser.add_argument("-o", "--output", dest="output", default=None, type=str, help="basename for Output GenBank file name (default: <input>)")
    parser.add_argument("-t", "--target", dest="target", required=True, help="Target (IAV or IBV)", choices=["IAV", "IBV"])
    parser.add_argument("-d", "--db_dir", dest="db_dir", help="Data path (optional; default: ganflu/db)", default=None)
    parser.add_argument("--isolate", dest="isolate", default=None, help='isolate name (e.g. "A/Narita/1/2009", "A/goose/Guangdong/1/1996", "B/Lee/1940")')
    parser.add_argument("--preserve_original_id", "--preserve-original-id", dest="preserve_original_id", action="store_true", help="Preserve original FASTA record IDs in GenBank output")
    parser.add_argument("--log-file", dest="log_file", default=None, help="Log file path (default: <output>.log)")
    parser.add_argument("--verbose", action="store_true", help="Show debug logs in the terminal")
    parser.add_argument("-v", "--version", action="version", version=_version())
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()

    return args 

def main():
    start_time = time.time()
    args = _get_args()

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

    log_file = os.path.abspath(args.log_file) if args.log_file else f"{out_stem}.log"
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
    try:
        logger.info("Running miniprot")
        miniprot = MiniprotCommandLine(
            input=input_fasta, work_dir=work_dir, output=gff3_file,
            prot_faa=ref_faa, miniprot_bin="miniprot", stderr_filename="miniprot.stderr", kmer_size=15
            )
        miniprot.run_piped_commands()
        logger.info(f"miniprot GFF3 output: {gff3_file}")

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
