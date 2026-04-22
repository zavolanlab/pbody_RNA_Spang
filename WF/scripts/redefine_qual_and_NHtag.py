#!/usr/bin/env python3

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import warnings
warnings.simplefilter('ignore')

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import subprocess
import pandas as pd
import numpy as np
import HTSeq

def write_to_new_bam(cur_alignment_list,cur_NH,bam_writer):
    if cur_NH==1:
        MAPQ=255
    else:
        # MAPQ = int(-10*np.log10(1-1/cur_NH)) # formula from STAR manual
        MAPQ=0 # useful to visualize in IGV, MM number is anyway accessible from NH tag
    for alignment in cur_alignment_list:
        alignment.aQual = MAPQ
        alignment.optional_fields = [(('NH',cur_NH) if elem[0]=='NH' else elem) for elem in alignment.optional_fields]
        bam_writer.write( alignment )

def main():
    """ 
        Parse dirty bam file and redefine correctly the quality string (STAR approach) and NH tag value
    """

    __doc__ = "Parse dirty bam file and redefine correctly the quality string (STAR approach) and NH tag value"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--input_bam_file",
                        dest="input_bam_file",
                        help="Path to the bam file",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--out_bam_file",
                        dest="out_bam_file",
                        help="path to the output bam file",
                        required=True,
                        metavar="DIR")    
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    input_bam_file = options.input_bam_file
    out_bam_file = options.out_bam_file

    almnt_file = HTSeq.SAM_Reader( input_bam_file )
    bam_writer = HTSeq.BAM_Writer.from_BAM_Reader( out_bam_file, almnt_file )
    
    # all .iv elements follow bed format
    j=0
    cur_name,cur_NH,cur_alignment_list = '',0,[]
    for almnt in almnt_file:
        read_name = almnt.read.name
        if cur_name!=read_name:
            # define NH and qual for the previous read group if not the start
            if cur_name!='':
                write_to_new_bam(cur_alignment_list,cur_NH,bam_writer)
                cur_alignment_list = []
            # next read group
            cur_NH = 0
            cur_alignment_list = []    
        cur_NH = cur_NH+1
        cur_name = read_name
        cur_alignment_list.append(almnt)
        j=j+1

    # if there are still unwritten alignments
    if len(cur_alignment_list)>0:
        write_to_new_bam(cur_alignment_list,cur_NH,bam_writer)
    bam_writer.close()
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)