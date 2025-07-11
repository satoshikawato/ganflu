#!/usr/bin/env python
# coding: utf-8

import os
import sys
import logging
import argparse
from importlib import resources
from .launchers.miniprot import MiniprotCommandLine
from .scripts import gff3togbk, validate_reference_files

rootlogger = logging.getLogger() 
rootlogger.handlers.pop() 

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

def _get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"ganflu v{_version()}: Influenza virus genome annotation")
    
    # Input/Output options
    parser.add_argument("-i", "--input", required=True, type=str, help="Input FASTA file")
    parser.add_argument("-o", "--output", dest="output", default=None, type=str, help="basename for Output GenBank file name (default: <input>)")
    parser.add_argument("-t", "--target", dest="target", required=True, help="Target (IAV or IBV)", choices=["IAV", "IBV"])
    parser.add_argument("-d", "--db_dir", dest="db_dir", help="Data path (optional; default: ganflu/db)", default=None)
    parser.add_argument("--isolate", dest="isolate", default=None, help='isolate name (e.g. "A/Narita/1/2009", "A/goose/Guangdong/1/1996", "B/Lee/1940")')
    parser.add_argument("-v", "--version", action="version", version=_version())
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    args = parser.parse_args()

    return args 

def main():
    logger = logging.getLogger()
    args = _get_args()

    input_fasta = args.input
    if args.output:
        out_stem = args.output
    else:
        out_stem = os.path.splitext(os.path.basename(input_fasta))[0]
    # workdir is the directory where the output files will be saved
    if args.output:
        work_dir = os.path.dirname(os.path.abspath(args.output))
    else:
        work_dir = os.path.dirname(os.path.abspath(input_fasta))

    target = args.target
    if args.db_dir:
        ref_dir = os.path.abspath(args.db_dir)
    else:
        db_path_traversable = resources.files('ganflu').joinpath(f'db/{args.target}')
        ref_dir = db_path_traversable.resolve()

    ref_toml = os.path.join(ref_dir, f"{target}.toml")
    ref_toml = validate_reference_files.validate_reference_files(target=args.target, db_dir=args.db_dir, logger=logger)
    ref_faa = os.path.join(ref_dir, ref_toml["metadata"]["prot_faa"])

    MiniprotCommandLine(
        input=input_fasta, work_dir=work_dir, output=f"{out_stem}.gff3",
        prot_faa=ref_faa, miniprot_bin="miniprot", stderr_filename="miniprot.stderr", kmer_size=15
        ).run_piped_commands()
    gff3togbk.main(["-g", f"{out_stem}.gff3", "-o", f"{out_stem}.gbk", "-i", input_fasta, "--toml", os.path.join(ref_dir, f'{target}.toml'), "--isolate", args.isolate]) 




    return 0

if __name__ == "__main__":
    main()
