#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import logging
import re

from .base import CommandLineTool

logger = logging.getLogger()
ANSI_ESCAPE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")

class MiniprotCommandLine(CommandLineTool):
    def __init__(self, input=None, work_dir=None, output=None, prot_faa=None, miniprot_bin="miniprot", stderr_filename="miniprot.stderr", kmer_size=15, prefix="MP"):
        super().__init__(miniprot_bin, stderr_filename, work_dir)
        self.miniprot_bin = miniprot_bin
        self.input = input
        self.output = output
        self.kmer_size = str(kmer_size)
        self.prot_faa = prot_faa
        self.prefix = prefix
        #miniprot -J 15 --gff A_duck_Japan_AQ-HE29-22_2017_H7N9.fa IAV_proteome_consensus.faa >A_duck_Japan_AQ-HE29-22_2017_H7N9.gff3
    def run_piped_commands(self):
        stderr_path = os.path.join(self.work_dir, self.stderr_file_name)
        cmdline = [self.bin_path, '-P', self.prefix, '--gff', '-J', self.kmer_size, self.input, self.prot_faa]
        logger.info(f"Running command: {' '.join(cmdline)}")
        logger.debug(f"miniprot stdout file: {self.output}")
        logger.debug(f"miniprot stderr file: {stderr_path}")
        with open(self.output, 'w') as output_file:
            with open(stderr_path, 'w') as stderr_file:
                miniprot_proc = subprocess.Popen(cmdline, stdout=output_file, stderr=stderr_file)
                miniprot_proc.communicate()
        if miniprot_proc.returncode != 0:
            logger.error(f"miniprot failed with exit code {miniprot_proc.returncode}")
            logger.error(f"miniprot stderr: {stderr_path}")
            try:
                with open(stderr_path, 'r') as stderr_file:
                    stderr_tail = stderr_file.readlines()[-20:]
                if stderr_tail:
                    stderr_message = ANSI_ESCAPE.sub("", "".join(stderr_tail).rstrip())
                    logger.error("Last lines from miniprot stderr:\n" + stderr_message)
            except OSError:
                logger.debug("Could not read miniprot stderr file", exc_info=True)
            raise RuntimeError(f"miniprot failed with exit code {miniprot_proc.returncode}")
        logger.info("miniprot completed successfully")
        return 0
