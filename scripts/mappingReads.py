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
    parser.add_argument('--ncores', help='Number of cores to use. Default=8', default='2')
    parser.add_argument('--temp_dir', help='Path to a temp dir folder. Default=none', default='none')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names_info.txt', default='sample_names.txt')
    args=parser.parse_args()

    params_file=args.analysis_info_file
    #path=functions.read_analysis_info_file(params_file)['Working directory']
    path=os.getcwd()
    refGenome=functions.read_analysis_info_file(params_file)['Reference Genome']
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

    # Detect if files are gz
    gz = functions.check_gz(in_dir)

    # Create temp dir in local node if specifyied
    temp_dir=args.temp_dir
    if temp_dir != 'none':
        functions.make_sure_path_exists(temp_dir)
        os.system("cp " + sample_names_file + " " + temp_dir )
        os.system("cp " + params_file + " " + temp_dir )
        os.chdir(temp_dir)
        os.system("sed '1d' " + params_file + " -i")
        os.system("echo 'Working directory = ' "+ temp_dir + " | cat - " + params_file  + " > " + "temp_file && mv " + "temp_file " + params_file)
        out_dir_original = out_dir
        out_dir= temp_dir + '/' +  args.out_dir
        functions.make_sure_path_exists(out_dir)

    # Run STAR
    Parallel(n_jobs=ncores)(delayed(mapping)(i) for i in sampleNames)
    Parallel(n_jobs=ncores)(delayed(indexing)(i) for i in sampleNames)
    #Parallel(n_jobs=ncores)(delayed(sortByName)(i) for i in sampleNames)

    # Move files if temp locar folder was created
    if temp_dir != 'none':
        os.system("mv " + out_dir + " " + path)



