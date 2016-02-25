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
sys.path.insert(0,'/usr/local/bin')
import functions
import argparse

def tables(i):
    outdir = re.sub('fastqc_data.txt', '', i)
    os.system("fastqc_plot_data.py " + i + " all " + outdir)


##############################################################

#os.system("preprocessing_numbers.R .")
#os.system("trimgalore_summary.R .")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Generates tables used for FastQC plots')
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=rawReads/', default='rawReads')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=rawReads/', default='rawReads')
    args=parser.parse_args()

    params_file=args.analysis_info_file
    path=functions.read_parameters_file(params_file)['Working directory']
    os.chdir(path)

    # Read sample names text file
    sampleNames = functions.read_sample_names()

    # Set input and output directories if not 'rawReads/'
    in_dir=args.in_dir
    out_dir=args.out_dir

    files=functions.get_filepaths(in_dir)
    files = [files[y] for y, x in enumerate(files) if re.findall("fastqc_data.txt", x)] 
    Parallel(n_jobs=8)(delayed(tables)(i) for i in files)


#os.system('ls rawReads/*/*fastqc  |  grep -v trimmed  | grep ":"  | sed \'s/://g\' > sample_names2.txt')
#os.system('fastqc_summary.py ./sample_names2.txt ./summary_fastqc.txt')


