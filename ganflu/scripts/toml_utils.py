
#!/usr/bin/env python
# coding: utf-8

import os
import sys
import platform
import toml
from collections import defaultdict
from importlib import resources
from importlib.abc import Traversable

def load_toml(logger, file_name=None, absolute_file_path=None) -> dict:
    toml_dict = {}  # Initialize to empty dict to handle cases where loading fails.
    if absolute_file_path and file_name:
        raise ValueError("Cannot specify both file name and absolute file path for load_toml()")

    elif absolute_file_path:  # Initialize outside of try for scope in exception block
        try:
            # To check if the file exists, do this:
            os.path.isfile(absolute_file_path)
            logger.info(f"INFO: Found TOML file: {absolute_file_path}")
            file_path = absolute_file_path
        except FileNotFoundError as e:
            logger.error(f"TOML file {absolute_file_path} does not exist: {e}")
            sys.exit(1)
    else:
        try:
            # Generate the path object for the 'config.toml' file
            file_path_traversable: Traversable = resources.files(
                'ganflu.data').joinpath(file_name)
            # Convert the path to an absolute path
            file_path = file_path_traversable.resolve()  # type: ignore
        except FileNotFoundError as e:
            logger.error(f"Config file {file_name} does not exist: {e}")
            sys.exit(1)
    try:
        with open(file_path, 'r') as toml_file:
            toml_dict: dict = toml.load(toml_file)
            logger.info(f"INFO: Loaded TOML file: {file_path}")
    except FileNotFoundError as e:
        logger.error(f"Failed to load TOML file {file_path}: {e}")
        sys.exit(1)
    return toml_dict    

def format_toml_dict(toml_dict):
    # Sort sections and keys
    sorted_dict = {}
    for section, content in toml_dict.items():
        if isinstance(content, dict):
            sorted_dict[section] = dict(sorted(content.items()))
        else:
            sorted_dict[section] = content

    return sorted_dict

def update_toml_file(toml_dict, save_file):
    formatted_dict = format_toml_dict(toml_dict)
    with open(save_file, "w") as save_file_handle:
        toml.dump(formatted_dict, save_file_handle)

def create_toml_dict(args, ref_toml, cfg):
    toml_dict = defaultdict(dict)
    toml_dict['system'] = {}
    toml_dict['system']['ganflu'] = cfg["version"]
    toml_dict['system']['os'] = platform.platform()
    toml_dict['system']['python'] = ".".join(map(str, sys.version_info[0:3]))
    toml_dict['params'] = {}
    toml_dict['time'] = {}
    toml_dict['stages'] = {}
    toml_dict['errors'] = []
    toml_dict['results'] = defaultdict(dict)
    
    toml_dict['warnings'] = defaultdict(dict)
    toml_dict['params']['cmd'] = " ".join(sys.argv)
    toml_dict['params']['input'] = args.input
    toml_dict['params']['out_dir'] = args.out_dir
    toml_dict['params']['target'] = args.target
    toml_dict['params']['isolate'] = args.isolate
    toml_dict['params']['config'] = cfg['config_file']
    toml_dict['params']['reference'] = ref_toml
    toml_dict['params']['db_dir'] = args.db_dir
    toml_dict['params']['primers'] = args.primers
    toml_dict['params']['num_threads'] = args.num_threads
    toml_dict['params']['debug'] = args.debug
    toml_dict['params']['unambiguous'] = args.unambiguous
    toml_dict['params']['skip_subtyping'] = args.skip_subtyping
    toml_dict['params']['delete_tmp_files'] = args.delete_tmp_files
    toml_dict['params']['keep_primers'] = args.keep_primers
    toml_dict['params']['skip_amplicon_filter'] = args.skip_amplicon_filter
    toml_dict['params']['min_depth'] = args.min_depth
    toml_dict['params']['min_read_len'] = args.min_read_len
    toml_dict['params']['max_depth_per_segment'] = args.max_depth_per_segment
    toml_dict['params']['cluster_identity_threshold'] = args.cluster_identity_threshold
    toml_dict['params']['cluster_size_threshold'] = args.cluster_size_threshold
    toml_dict['params']['mashmap_min_identity'] = args.mashmap_min_identity
    toml_dict['params']['mashmap_kmer_size'] = args.mashmap_kmer_size
    toml_dict["results"]["variants"] = []
    toml_dict["results"]["duplicated_sequences"] = []
    
    return toml_dict