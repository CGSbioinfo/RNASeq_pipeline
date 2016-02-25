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


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raisec

def get_filepaths(directory):
    """
    This function will generate the file names in a directory
    tree by walking the tree either top-down or bottom-up. For each
    directory in the tree rooted at directory top (including top itself),
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths

def read_parameters_file(params_file):
    params_file = params_file
    args = {}
    with open(params_file, 'r') as f:
        for line in f:
            entry = line.strip().split("=")
            if entry[0]:
                args[entry[0].strip(" ")] = entry[1].strip()
    global path
    path = args['Working directory'].replace("\\", "")
    global gtfFile
    gtfFile = args['GTF File']
    global refGenome
    refGenome = args['Reference Genome']
    global bedFile
    bedFile = args['BedFile']
    global bedFile_10k
    bedFile_10k = args['BedFile10K']
    global refFlat
    refFlat = args['refFlat']
    global rRNA_interval_list
    rRNA_interval_list = args['rRNA_interval_list']
    global strand
    strand = args['strand']
    global reads_dir
    reads_dir = args['reads_dir']
    return(args)

def create_rawReads_folder(sampleNames):
    folders = os.listdir('.')
    readsFiles = [folders[i] for i, x in enumerate(folders) if re.findall('rawReads',x)]
    allFiles = get_filepaths(reads_dir)
    allFiles= [allFiles[i] for i, x in enumerate(allFiles) if re.findall("_R\d.*.fastq", x)]
    if not readsFiles:
        make_sure_path_exists('rawReads')
        sampleDir = []
        for sample in sampleNames:
            reads = [allFiles[i] for i,x in enumerate(allFiles) if re.findall(sample,x)]
            if sample not in sampleDir:
                make_sure_path_exists('rawReads/'+sample)
            for r in reads:
               os.system('mv ' + '"' + r + '"' + ' rawReads/' + sample)
            sampleDir.append(sample)

def read_sample_names():
    sampleNames = []
    sample_names_file = open('sample_names.txt','r')
    for line in sample_names_file:
        sampleNames.append(line.strip())
    return(sampleNames)

def check_gz(dir):
    readFiles = []
    for root, dir, files in os.walk(dir):
        readFiles.extend(files)
    indicesgzFiles = [i for i, x in enumerate(readFiles) if re.findall(".fastq.gz", x)]
    if indicesgzFiles:
        gz =".gz"
    else:
        gz =""
    return(gz)


