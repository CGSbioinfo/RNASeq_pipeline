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
        os.system("trim_galore --gzip --paired --fastqc --fastqc_args '--nogroup --extract' --output_dir " + out_dir + " " + 
            in_dir + i + "/" + i + "*_R1*.fastq" + gz + " " + 
            in_dir + i + "/" + i + "*_R2*.fastq" + gz)
    else:
        os.system("trim_galore --gzip --fastqc --fastqc_args '--nogroup --extract' --output_dir " + out_dir + " " + in_dir + i + "/" + i + "*_R1*.fastq" + gz)


#########################


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'QC and adapter trimming using Trim Galore')
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=rawReads/', default='rawReads/')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=trimmedReads/', default='trimmedReads/')
    parser.add_argument('--out_dir_report', help='Path to out put folder. Default=Report/figure/data/', default='Report/figure/data/')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names_info.txt', default='sample_names.txt')
    parser.add_argument('--ncores', help='Number of cores to use. Default=8', default='8')
    args=parser.parse_args()

    # Set path of working directory
    params_file=args.analysis_info_file
    #path=functions.read_analysis_info_file(params_file)['Working directory']
    path=os.getcwd()
    os.chdir(path)
    
    #Ncores
    ncores=int(args.ncores)

    # Read sample names text file
    sample_names_file=args.sample_names_file
    sampleNames = functions.read_sample_names(sample_names_file)

    # Set input and output directories if not 'rawReads/'
    in_dir=path + '/' + args.in_dir
    out_dir=path + '/' + args.out_dir
    out_dir_report=path + '/' + args.out_dir_report

    # Detect if files are gz
    gz = functions.check_gz(in_dir)

    # Run trimm galore
    functions.make_sure_path_exists(out_dir)
    Parallel(n_jobs=7)(delayed(trimming)(i) for i in sampleNames)
    functions.make_sure_path_exists(out_dir_report)
    
    # Nreads
    os.system("Rscript /usr/local/bin/trimming_summary.R " + in_dir + " " + out_dir + " " + out_dir_report )

