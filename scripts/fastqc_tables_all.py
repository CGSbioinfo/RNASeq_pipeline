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

def plots(i):
    outdir = re.sub('fastqc_data.txt', '', i)
    os.system('Rscript /usr/local/bin/fastqc_plots_all.R ' + in_dir + ' ' + i + ' ' + readType + ' ' + 
        out_dir_plots + ' ' + suffix_name  )


##############################################################

#os.system("preprocessing_numbers.R .")
#os.system("trimgalore_summary.R .")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Generates tables used for FastQC plots')
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=rawReads/', default='rawReads')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=rawReads/', default='rawReads')
    parser.add_argument('--readType', help='Read Type: pairedEnd or singleEnd. Default=pairedEnd', default='pairedEnd')
    parser.add_argument('--out_dir_plots', help='Path to out put folder. Default=Report/figure', default='Report/figure')
    parser.add_argument('--suffix_name', help='Suffix to optionally put to the output name. Default=', default='')
    args=parser.parse_args()

    params_file=args.analysis_info_file
    path=functions.read_parameters_file(params_file)['Working directory']
    os.chdir(path)

    # Read sample names text file
    sampleNames = functions.read_sample_names()

    # Set input and output directories if not 'rawReads/'
    in_dir=args.in_dir
    out_dir=args.out_dir
    out_dir_plots=args.out_dir_plots
    readType=args.readType
    suffix_name=args.suffix_name
    
    files=functions.get_filepaths(in_dir)
    files = [files[y] for y, x in enumerate(files) if re.findall("fastqc_data.txt", x)] 
    Parallel(n_jobs=8)(delayed(tables)(i) for i in files)
    
    functions.make_sure_path_exists(out_dir_plots)
    Parallel(n_jobs=8)(delayed(plots)(i) for i in sampleNames)

    


#os.system('ls rawReads/*/*fastqc  |  grep -v trimmed  | grep ":"  | sed \'s/://g\' > sample_names2.txt')
#os.system('fastqc_summary.py ./sample_names2.txt ./summary_fastqc.txt')


