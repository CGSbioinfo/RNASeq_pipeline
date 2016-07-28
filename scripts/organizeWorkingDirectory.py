#!/usr/bin/env python

import os
import sys
import re
import errno
import glob
import time
import pickle
import logging
from joblib import Parallel, delayed
import multiprocessing
import subprocess
sys.path.insert(0,'/usr/local/bin/')
import functions
import argparse

__version__='v02'

if __name__ == '__main__':

    # Parser
    parser = argparse.ArgumentParser(prog='organizeWorkingDirectory.py',description = 'Organize working directory of the analysis')
    parser.add_argument('-v','--version', action='version',version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names_info.txt', default='sample_names.txt')
    args=parser.parse_args()

    # Set path of working directory
    params_file=args.analysis_info_file
    path=functions.read_analysis_info_file(params_file)['Working directory']
    reads_dir=functions.read_analysis_info_file(params_file)['reads_dir']
    os.chdir(path)

    # Read sample names
    sample_names_file=args.sample_names_file
    sampleNames = functions.read_sample_names(sample_names_file)
     
    # Create rawReads folder
    functions.create_rawReads_folder(sampleNames)
