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


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = 'Creating differential expression files')
	parser.add_argument('--out_dir', help='Path to output folder. Default=./', default='./')
	parser.add_argument('--deseq2_arguments', help='Output text file with deseq2 arguments. Example=deseq2_arguments.txt')
	parser.add_argument('--edger_arguments', help='Output text file with edgeR arguments. Example=edger_arguments.txt')
	parser.add_argument('--sample_info', help='Output csv file to fill with sample names and groups. Example=sample_info.csv')
	parser.add_argument('--comparisons', help='Output csv file to fill with comparisons to be performed. Example=comparisons.csv')
	args=parser.parse_args()

	if args.deseq2_arguments:
		outfile_name=args.deseq2_arguments
		lines=['indir = ', 'outdir = ', 'sample_info = ', 'comparisons = ', 'design = ', 'gtfFile = ']
		outfile=open(outfile_name,'w')
		for l in lines:
			outfile.write(l + '\n')
		outfile.close()

	if args.edger_arguments:
		outfile_name=args.edger_arguments
		lines=['indir = ', 'outdir = ', 'sample_info = ', 'comparisons = ', 'min.count = ', 'min.nsamples = ' ,'design = ', 'gtfFile = ']
		outfile=open(outfile_name,'w')
		for l in lines:
			outfile.write(l + '\n')
		outfile.close()

	if args.sample_info:
		outfile_name=args.sample_info
		lines=['SampleID,Group']
		outfile=open(outfile_name,'w')
		for l in lines:
			outfile.write(l + '\n')
		outfile.close()

	if args.comparisons:
		outfile_name=args.comparisons
		lines=['baselineGroup,comparisonGroup']
		outfile=open(outfile_name,'w')
		for l in lines:
			outfile.write(l + '\n')
		outfile.close()



    
