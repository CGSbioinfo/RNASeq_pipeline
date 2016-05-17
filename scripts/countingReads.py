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

def counting(i):
    os.system('samtools view ' + in_dir + '/' + i + 'Aligned.sortedByCoord.sortedByName.out.bam | htseq-count -a 10 -m union -s ' + strand_htseq + ' - ' + gtfFile + ' > ' + out_dir + '/' + i + '_count.txt')


#########################


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Counting Reads')
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=alignedReads/', default='alignedReads/')
    parser.add_argument('--out_dir', help='Path to output folder. Default=countedReads/', default='countedReads/')
    parser.add_argument('--mapping_summary_file', help='Mapping summary file. Default=mapping_summary.csv', default='mapping_summary.csv')

    args=parser.parse_args()

    params_file=args.analysis_info_file
    path=functions.read_parameters_file(params_file)['Working directory']
    refGenome=functions.read_parameters_file(params_file)['Reference Genome']
    strand=functions.read_parameters_file(params_file)['strand']
    strand_piccard, strand_htseq = functions.get_strand(strand)
    gtfFile=functions.read_parameters_file(params_file)['GTF File']

    os.chdir(path)

    # Read sample names text file
    sampleNames = functions.read_sample_names()

    # Set input and output directories if not '/'
    in_dir=args.in_dir
    out_dir=args.out_dir
    functions.make_sure_path_exists(out_dir)
    mapping_summary_file=args.mapping_summary_file


    # Detect if files are gz
    gz = functions.check_gz(in_dir)

    # Count command
    Parallel(n_jobs=7)(delayed(counting)(i) for i in sampleNames)
    
    # QC
    os.system("Rscript /usr/local/bin/countsLog_rnaseq.R " + out_dir + ' ' + mapping_summary_file)
    os.system("Rscript /usr/local/bin/library_proportion.R " + out_dir + ' ' + out_dir + ' ' + gtfFile)

