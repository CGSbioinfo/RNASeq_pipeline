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

def trimming(i):
    allFiles = os.listdir(in_dir + "/" + i )
    pairedReads_temp = [allFiles[y] for y, x in enumerate(allFiles) if re.findall("_R2", x)]
    if pairedReads_temp:
        os.system("trim_galore --paired --fastqc --fastqc_args '--nogroup --extract' --output_dir " + out_dir + " " + 
            in_dir + i + "/" + i + "*_R1*.fastq" + gz + " " + 
            in_dir + i + "/" + i + "*_R2*.fastq" + gz)
    else:
        os.system("trim_galore --fastqc --fastqc_args '--nogroup --extract' --output_dir " + out_dir + " " + in_dir + i + "/" + i + "*_R1*.fastq" + gz)


#########################


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'QC and adapter trimming using Trim Galore')
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=rawReads/', default='rawReads/')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=trimmedReads/', default='trimmedReads/')
    args=parser.parse_args()

    # Read analysis info file
    params_file=args.analysis_info_file
    path=functions.read_parameters_file(params_file)['Working directory']
    os.chdir(path)

    # Read sample names text file
    sampleNames = functions.read_sample_names()

    # Set input and output directories if not 'rawReads/'
    in_dir=args.in_dir
    out_dir=args.out_dir

    # Detect if files are gz
    gz = functions.check_gz(in_dir)

    functions.make_sure_path_exists(out_dir)
    Parallel(n_jobs=7)(delayed(trimming)(i) for i in sampleNames)


