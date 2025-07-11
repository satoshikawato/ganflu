#!/usr/bin/env python
# coding: utf-8

import os
import subprocess

from .base import CommandLineTool

class MiniprotCommandLine(CommandLineTool):
    def __init__(self, input=None, work_dir=None, output=None, prot_faa=None, miniprot_bin="miniprot", stderr_filename="miniprot.stderr", kmer_size=15, prefix="MP"):
        super().__init__(miniprot_bin, stderr_filename, work_dir)
        self.miniprot_bin = miniprot_bin
        self.input = input
        self.output = output
        self.kmer_size = str(kmer_size)
        self.prot_faa = prot_faa
        self.prefix = prefix
        self.run_piped_commands()
        #miniprot -J 15 --gff A_duck_Japan_AQ-HE29-22_2017_H7N9.fa IAV_proteome_consensus.faa >A_duck_Japan_AQ-HE29-22_2017_H7N9.gff3
    def run_piped_commands(self):
        stderr_file = os.path.join(self.work_dir, self.stderr_file_name)
        with open(self.output, 'w') as output_file:
            with open(stderr_file, 'w') as stderr_file:
                miniprot_proc = subprocess.Popen([self.miniprot_bin, '-P', self.prefix, '--gff', '-J', self.kmer_size, self.input, self.prot_faa], stdout=output_file, stderr=stderr_file)
                miniprot_proc.communicate()
        return 0