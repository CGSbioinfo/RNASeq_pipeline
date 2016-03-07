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

def mapping(i):
    trimmedReads = os.listdir(in_dir)
    trimmedReads = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall(i, x)]
    r2 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall("_R2.*val_2.fq", x)]
    if r2:
        r1 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall( "_R1.*val_1.fq", x)]
        r2 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall("_R2.*val_2.fq", x)]
        os.system("STAR --genomeDir " + refGenome + 
            " --readFilesIn " + in_dir + '/' + r1[0] + " " + in_dir + '/' + r2[0] + 
            " --runThreadN 4 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat " + 
            " --outFileNamePrefix " + out_dir + '/' + i)
    else:
        r1 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall( "_R1.*trimmed.fq", x)]
        os.system("STAR --genomeDir " + refGenome + 
            " --readFilesIn " + in_dir + "/" + r1[0] + 
            " --runThreadN 8 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat " + 
            " --outFileNamePrefix " + out_dir + '/' + i)


def indexing(i):
     os.system("samtools index " + out_dir + '/' + i + "Aligned.sortedByCoord.out.bam")


def sortByName(i):
     os.system('samtools sort -n ' + out_dir + '/' + i + 'Aligned.sortedByCoord.out.bam ' + out_dir + '/' + i + 'Aligned.sortedByCoord.sortedByName.out')


#########################


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Mapping reads with STAR')
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=trimmedReads/', default='trimmedReads/')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=alignedReads/', default='alignedReads/')
    args=parser.parse_args()

    params_file=args.analysis_info_file
    path=functions.read_parameters_file(params_file)['Working directory']
    refGenome=functions.read_parameters_file(params_file)['Reference Genome']
    os.chdir(path)

    # Read sample names text file
    sampleNames = functions.read_sample_names()

    # Set input and output directories if not '/'
    in_dir=args.in_dir
    out_dir=args.out_dir
    functions.make_sure_path_exists(out_dir)

    # Detect if files are gz
    gz = functions.check_gz(in_dir)

    Parallel(n_jobs=2)(delayed(mapping)(i) for i in sampleNames)
    Parallel(n_jobs=7)(delayed(indexing)(i) for i in sampleNames)
    Parallel(n_jobs=8)(delayed(sortByName)(i) for i in sampleNames)


