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

def qc_check(i):
    allFiles = os.listdir(in_dir + "/" + i )
    pairedReads_temp = [allFiles[y] for y, x in enumerate(allFiles) if re.findall("_R2", x)]
    functions.make_sure_path_exists(out_dir+'/'+i)
    os.system("fastqc " + in_dir + "/" + i + "/" + i + "*_R1*.fastq" + gz + " --outdir=" + out_dir + "/" + i + " --nogroup --extract")
    if pairedReads_temp:
        os.system("fastqc " + path + "/rawReads/" + i + "/" + i + "*_R2*.fastq" + gz + " --outdir=" + out_dir + "/" + i + " --nogroup --extract")

####################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Quality control of data using FastQC')
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=rawReads/', default='rawReads/')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=rawReads/', default='rawReads/')
    args=parser.parse_args()

    # Set path of working directory
    params_file=args.analysis_info_file
    path=functions.read_parameters_file(params_file)['Working directory']
    os.chdir(path)
    
    # Read sample names text file
    sampleNames = functions.read_sample_names()

    # Set input and output directories if not 'rawReads/'
    in_dir=args.in_dir
    out_dir=args.out_dir
    #try:
    #    in_dir=sys.argv[2]
    #    organizeWorkingDirectory.make_sure_path_exists(in_dir)
    #except:
    #    in_dir='rawReads/'

    #try:
    #    out_dir=sys.argv[3]
    #    organizeWorkingDirectory.make_sure_path_exists(out_dir)
    #except:
    #    out_dir='rawReads/'

    # Detect if files are gz 
    gz = functions.check_gz(in_dir)

    Parallel(n_jobs=7)(delayed(qc_check)(i) for i in sampleNames)










