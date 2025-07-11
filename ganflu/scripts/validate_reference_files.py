#!/usr/bin/env python
# coding: utf-8

import os
import sys
import toml
import logging
from importlib import resources
from importlib.abc import Traversable

logger = logging.getLogger("ganflu")

def find_reference_toml(db_dir, logger, extension="toml"):
    # Get a list of all files in the directory
    files = os.listdir(db_dir)

    # Filter files with the specified extension
    toml_files = [file for file in files if file.endswith(extension)]

    if len(toml_files) == 1:
        # If there's only one file with the specified extension
        ref_toml = os.path.join(db_dir, toml_files[0])
        logger.info(f"Reference TOML file found: {ref_toml}")
        with open(ref_toml, 'r') as file:
            ref_toml = toml.load(file)
        return ref_toml
    else:
        # If there are no files or more than one file with the specified extension
        if len(toml_files) == 0:
            logger.error(f"No TOML file found in the directory: {db_dir}")
        else:
            logger.error(f"Multiple TOML files found in the directory: {db_dir}")
        sys.exit(1)

def create_faidx_if_needed(config_dict, db_dir, target, logger):
    # Load the reference TOML configuration
    ref_file = config_dict["ref_toml"]
    with open(ref_file, 'r') as file:
        ref_config = toml.load(file)
    db_nucl_dir = os.path.join(db_dir, "nucl")
    # Iterate over each virus in the reference configuration
    #for virus, virus_config in ref_config.items():
    virus_config = ref_config
    # Check if the virus has a "segments" section
    try:
        virus_config["segments"]
        # Iterate over each segment in the "segments" section
        for key in virus_config["segments"].keys():
            # Get the FASTA file path for the segment
            segment = virus_config["segments"][key]
            fasta_file = segment["file"]
            absolute_file_path = os.path.join(db_nucl_dir, fasta_file)
            # Check if the FASTA file exists
            try:
                os.path.isfile(absolute_file_path)
                # Get the FASTA index file path
                faidx_file = f"{absolute_file_path}.fai"
                # Check if the FASTA index file exists
                if os.path.isfile(faidx_file):
                    # Get the last edit time of the FASTA file and FASTA index file
                    fasta_mtime = os.path.getmtime(absolute_file_path)
                    faidx_mtime = os.path.getmtime(faidx_file)
                    # Check if the FASTA index file is older than the FASTA file
                    if faidx_mtime < fasta_mtime:
                        # Recreate the FASTA index file using samtools faidx
                        command = f"samtools faidx {absolute_file_path}"
                        os.system(command)
                        logger.info(f"Recreated FASTA index file: {faidx_file}")
                else:
                    # Create the FASTA index file using samtools faidx
                    command = f"samtools faidx {absolute_file_path}"
                    os.system(command)
                    logger.info(f"Created FASTA index file: {faidx_file}")
            except FileNotFoundError:
                logger.error(f"FASTA file not found: {absolute_file_path}")
    except KeyError:
        logger.warning(f"No 'segments' section found for target: {target}")
        sys.exit(1)

def create_fasta_list_if_needed(config_dict, db_dir, target, logger):
    # Load the reference TOML configuration
    ref_file = config_dict["ref_toml"]
    with open(ref_file, 'r') as file:
        ref_config = toml.load(file)
    db_nucl_dir = os.path.join(db_dir, "nucl")
    # Iterate over each virus in the reference configuration
    #for virus, virus_config in ref_config.items():
    virus_config = ref_config
    # Check if the virus has a "segments" section
    try:
        virus_config["segments"]
        # Get the FASTA list file path
        fasta_list_file = os.path.join(db_nucl_dir, "fasta_list.txt")
        # Check if the FASTA list file exists
        if os.path.isfile(fasta_list_file):

            # Get the last edit time of the FASTA list file
            fasta_list_mtime = os.path.getmtime(fasta_list_file)
            # Iterate over each segment in the "segments" section
            for key in virus_config["segments"].keys():
                # Get the FASTA file path for the segment
                segment = virus_config["segments"][key]
                fasta_file = segment["file"]
                absolute_file_path = os.path.join(db_nucl_dir, fasta_file)
                # Check if the FASTA file exists
                try:
                    os.path.isfile(absolute_file_path)
                    # Get the last edit time of the FASTA file
                    fasta_mtime = os.path.getmtime(absolute_file_path)
                    # Check if the FASTA list file is older than the FASTA file
                    if fasta_list_mtime < fasta_mtime:
                        # Create the FASTA list file
                        with open(fasta_list_file, 'w') as file:
                            # Iterate over each segment in the "segments" section
                            for key in virus_config["segments"].keys():
                                # Get the FASTA file path for the segment
                                segment = virus_config["segments"][key]
                                fasta_file = segment["file"]
                                file.write(f"{absolute_file_path}\n")
                        logger.info(f"Created FASTA list file: {fasta_list_file}")
                except FileNotFoundError:
                    logger.error(f"FASTA file not found: {absolute_file_path}")
        else:
            # Create the FASTA list file
            with open(fasta_list_file, 'w') as file:
                # Iterate over each segment in the "segments" section
                for key in virus_config["segments"].keys():
                    # Get the FASTA file path for the segment
                    segment = virus_config["segments"][key]
                    fasta_file = segment["file"]
                    absolute_file_path = os.path.join(db_nucl_dir, fasta_file)
                    file.write(f"{absolute_file_path}\n")
            logger.info(f"Created FASTA list file: {fasta_list_file}")
    except KeyError:
        logger.warning(f"No 'segments' section found for target: {target}")
        sys.exit(1)


def validate_reference_files(target="", db_dir="", logger=""):
    if db_dir:
        if os.path.exists(db_dir):
            logger.info(f"INFO: Found reference directory: {db_dir}")
        else:
            logger.error(f"Reference directory {db_dir} does not exist")
            sys.exit(1)
    else:
        try:
            db_path_traversable: Traversable = resources.files(
                'ganflu.db').joinpath(target)
            db_dir = db_path_traversable.resolve()
        except FileNotFoundError as e:
            logger.error(f"Reference directory {db_dir} does not exist: {e}")
            sys.exit(1)
    # Check the directory for the required files
    ref_toml = find_reference_toml(db_dir, logger)

    return ref_toml