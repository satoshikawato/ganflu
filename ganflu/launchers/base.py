#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import logging
import subprocess
import shutil
import logging

logger = logging.getLogger()

class CommandLineTool:
    def __init__(self, binary_name, stderr_file_name, work_dir):
        self.bin_path = self.check_binaries(binary_name)
        self.stderr_file_name = stderr_file_name
        self.cmdline = [self.bin_path]
        self.work_dir = work_dir
    def check_binaries(self, binary_name):
        if not shutil.which(binary_name):
            raise Exception(f"Binary {binary_name} was not found. ")
        try:
            # Check if the binary is executable. If found, return the path to the binary. Also displays the path in debug message
            logger.debug(f"Found {binary_name} at {shutil.which(binary_name)}")
            return shutil.which(binary_name)
        except subprocess.CalledProcessError as e:
            if e.returncode == -9:
                logger.error("Looks like the system ran out of memory")
            raise Exception(str(e))
        except OSError as e:
            raise Exception(str(e))

    def construct_query(self):
        pass
    def run_piped_commands(self):
        pass
    def run(self):
        stderr_file = os.path.join(self.work_dir, self.stderr_file_name)
        try:
            # Run the command and capture the output. Cmd is displayed in debug message
            message=f"Running command: {' '.join(self.cmdline)}"
            logger.info(message)
            cmdline = subprocess.Popen(self.cmdline, stderr=open(stderr_file, "w"), shell=False)
            cmdline.communicate()
        except (subprocess.CalledProcessError, OSError) as e:
            logger.error(f"Error running {self.__class__.__name__}, terminating. See the alignment error log for details: " + stderr_file)
            logger.error("Cmd: " + " ".join(self.cmdline))
            raise Exception(str(e))

