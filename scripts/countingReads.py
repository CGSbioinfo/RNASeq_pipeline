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
    os.system('samtools sort -n --output-fmt sam ' + in_dir + '/' + i + 'Aligned.sortedByCoord.out.bam | htseq-count -a 10 -m union -s ' + strand_htseq + ' - ' + gtfFile + ' > ' + out_dir + '/' + i + '_count.txt')


#########################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Counting Reads')
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=alignedReads/', default='alignedReads/')
    parser.add_argument('--out_dir', help='Path to output folder. Default=countedReads/', default='countedReads/')
    parser.add_argument('--mapping_summary_file', help='Mapping summary file. Default=alignedReads/QC/mapping_summary.csv', default='alignedReads/QC/mapping_summary.csv')
    parser.add_argument('--out_dir_report', help='Path to out put folder. Default=Report/figure/countingQC/', default='Report/figure/countingQC/')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names_info.txt', default='sample_names.txt')
    parser.add_argument('--ncores', help='Number of cores to use. Default=8', default='8')
    args=parser.parse_args()

    params_file=args.analysis_info_file
    #path=functions.read_analysis_info_file(params_file)['Working directory']
    path=os.getcwd()
    refGenome=functions.read_analysis_info_file(params_file)['Reference Genome']
    strand=functions.read_analysis_info_file(params_file)['strand']
    strand_piccard, strand_htseq = functions.get_strand(strand)
    gtfFile=functions.read_analysis_info_file(params_file)['GTF File']
    os.chdir(path)
    

    #Ncores
    ncores=int(args.ncores)

    # Read sample names text file
    sample_names_file=args.sample_names_file
    sampleNames = functions.read_sample_names(sample_names_file)

    # Set input and output directories if not '/'
    in_dir=path + '/' + args.in_dir
    out_dir=path + '/' + args.out_dir
    functions.make_sure_path_exists(out_dir)
    mapping_summary_file=path + '/' + args.mapping_summary_file
    out_dir_report=path + '/' + args.out_dir_report
    functions.make_sure_path_exists(out_dir_report)

    # Detect if files are gz
    gz = functions.check_gz(in_dir)

    # Count command
    Parallel(n_jobs=ncores)(delayed(counting)(i) for i in sampleNames)
    
    # QC
    os.system("Rscript /usr/local/bin/countsLog_rnaseq.R " + out_dir + ' ' +  out_dir_report + ' ' + mapping_summary_file)
    os.system("Rscript /usr/local/bin/library_proportion.R " + out_dir + ' ' + out_dir_report + ' ' + gtfFile)

