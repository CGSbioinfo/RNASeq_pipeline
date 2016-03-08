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

def junctions(i):
    os.system("python ~/bin/junction_annotation.py -i " + 
        in_dir + "/" + i +"Aligned.sortedByCoord.out.bam -o " + 
        out_dir + '/' + i + " -r " + bedFile + " & ")
    os.system("python ~/bin/junction_saturation.py -i " + 
        in_dir + '/' + i +"Aligned.sortedByCoord.out.bam  -o " + 
        out_dir + '/' + i + " -r " + bedFile)


# Collect Metrics
def piccard_collect_metrics(i):
    bamfiles = os.listdir(in_dir)
    bamfiles = [bamfiles[y] for y, x in enumerate(bamfiles) if re.findall(i, x)]
    bamfiles = [bamfiles[y] for y, x in enumerate(bamfiles) if re.findall(r'sortedByCoord\.out\.bam$', x)][0]
    os.system("java -jar ~/tools/picard-tools-1.127/picard.jar CollectRnaSeqMetrics " + 
        " REF_FLAT=" + refFlat + 
        " RIBOSOMAL_INTERVALS=" + rRNA_interval_list + 
        " STRAND_SPECIFICITY=" + strand_piccard + 
        " INPUT=" + in_dir + "/" + bamfiles +  
        " OUTPUT=" + out_dir + "/" + i + "_metrics.txt")


def pct(i):
    os.system('mapping_distribution.R ' + out_dir + '/ ' + i)


#########################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Quality control of mapped data')
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=alignedReads/', default='alignedReads/')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=alignedReads/QC/', default='alignedReads/QC/')
    parser.add_argument('--out_dir_plots', help='Path to out put folder to copy plots. Default=Report/figure/', default='Report/figure/')
    args=parser.parse_args()


    params_file=args.analysis_info_file
    path=functions.read_parameters_file(params_file)['Working directory']
    refGenome=functions.read_parameters_file(params_file)['Reference Genome']
    bedFile_10k=functions.read_parameters_file(params_file)['BedFile10K']
    bedFile=functions.read_parameters_file(params_file)['BedFile']
    refFlat=functions.read_parameters_file(params_file)['refFlat']
    rRNA_interval_list=functions.read_parameters_file(params_file)['rRNA_interval_list']
    strand=functions.read_parameters_file(params_file)['strand']
    strand_piccard, strand_htseq = functions.get_strand(strand)
    #if strand == 'reverse':
    #    strand = 'SECOND_READ_TRANSCRIPTION_STRAND'
    #os.chdir(path)

    # Read sample names text file
    sampleNames = functions.read_sample_names()

    # Set input and output directories if not '/'
    in_dir=args.in_dir
    out_dir=args.out_dir
    functions.make_sure_path_exists(out_dir)
    out_dir_plots=args.out_dir_plots
    functions.make_sure_path_exists(out_dir_plots)

    # Detect if files are gz
    gz = functions.check_gz(in_dir)

    os.system("mapping_summary.R " + in_dir + '/')

    # Gene body coverage
    os.system("ls " + in_dir + "/*Aligned.sortedByCoord.out.bam > tempbamfiles.txt")
    os.system("python ~/bin/geneBody_coverage.py -r " + bedFile_10k + " -i tempbamfiles.txt -o " + out_dir + "/10KGenes")
    os.system("rm tempbamfiles.txt")
    os.system('cp ' + out_dir + '/10KGenes.geneBodyCoverage.curves.pdf ' + out_dir_plots + '/10KGenes_geneBodyCoverage_curves.pdf')

    # Junction QC
    Parallel(n_jobs=8)(delayed(junctions)(i) for i in sampleNames)
    os.system('grep "y=c(" ' + out_dir + '/*junctionSaturation*  | sed \'s/:y=c(/,/g\' | sed \'s/.junctionSaturation_plot.r//g\' | sed \'s/)//g\' | sed \"s/.*\///g\"  > ' + out_dir + '/junctionSat_all.csv')
    os.system('Rscript ~/bin/junctionPlotAll.R ' + out_dir + ' ' + out_dir)
    os.system('cp ' + out_dir + '/junctionSaturationAll.pdf ' + out_dir_plots)

    # Piccard tools
    Parallel(n_jobs=8)(delayed(piccard_collect_metrics)(i) for i in sampleNames)
    Parallel(n_jobs=7)(delayed(pct)(i) for i in sampleNames)

