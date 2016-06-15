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
    os.system("fastqc "  + in_dir + "/" + i + "/" + i + "*_R1*.fastq" + gz + " --outdir=" + out_dir + "/" + i + " --nogroup --extract ")
    if pairedReads_temp:
        os.system("fastqc " + in_dir + "/" + i + "/" + i + "*_R2*.fastq" + gz + " --outdir=" + out_dir + "/" + i + " --nogroup --extract")

####################

if __name__ == '__main__':

    """ This script will run fastqc on all samples in parallel. 
     - You can specify the number of cores used
     - The input and output default is rawReads/samplesubfolder/
     - It also returns a table and a plot with the number of reads for each sample, the output by default is Report/figure/data """

    # Parser
    parser = argparse.ArgumentParser(description = 'Quality control of data using FastQC')
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=rawReads/', default='rawReads/')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=rawReads/', default='rawReads/')
    parser.add_argument('--out_dir_report', help='Path to out put folder. Default=Report/figure/data/', default='Report/figure/data/')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names_info.txt', default='sample_names.txt')
    parser.add_argument('--ncores', help='Number of cores to use. Default=8', default='8')
    args=parser.parse_args()

    # Set path of working directory
    params_file=args.analysis_info_file
    path=functions.read_analysis_info_file(params_file)['Working directory']
    os.chdir(path)
    
    #Ncores
    ncores=int(args.ncores)

    # Read sample names text file
    sample_names_file=args.sample_names_file
    sampleNames = functions.read_sample_names(sample_names_file)

    # Set input and output directories if not 'rawReads/'
    in_dir=path + '/' + args.in_dir
    out_dir=path + '/' +args.out_dir
    out_dir_report=path + '/' + args.out_dir_report

    # Create out_dir_report
    functions.make_sure_path_exists(out_dir_report)

    # Detect if files are gz 
    gz = functions.check_gz(in_dir)

    # Run fastqc
    Parallel(n_jobs=ncores)(delayed(qc_check)(i) for i in sampleNames)

    # Number of reads per sample
    os.system("Rscript /usr/local/bin/preprocessing_numbers.R " + in_dir + " " + out_dir_report) 








