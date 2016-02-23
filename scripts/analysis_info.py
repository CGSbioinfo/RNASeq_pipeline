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


if __name__ == '__main__':
    try: 
        outfile_name=sys.argv[1]
    except:
        outfile_name='analysis_info.txt'
    lines=['Working directory = ', 'reads_dir = ', 'Reference Genome = ', 'GTF File = ', 'BedFile = ', 'BedFile10K = ', 'refFlat = ', 'rRNA_interval_list = ', 'strand = ']
    outfile=open(outfile_name,'w')
    for l in lines:
        outfile.write(l + '\n')
    outfile.close()
    
