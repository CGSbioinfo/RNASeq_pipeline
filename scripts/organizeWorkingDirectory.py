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
sys.path.insert(0,'/usr/local/bin')
import functions
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Organize working directory of the analysis')
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')

    args=parser.parse_args()

    # Set path of working directory
    params_file=args.analysis_info_file
    path=functions.read_parameters_file(params_file)['Working directory']
    os.chdir(path)

    sampleNames = functions.read_sample_names()

    functions.create_rawReads_folder(sampleNames)


