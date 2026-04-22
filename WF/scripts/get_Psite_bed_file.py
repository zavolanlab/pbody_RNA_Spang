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
import csv
import re

def main():
    """ Parse bed file obtained with bedtools bamtobed -cigar
    """

    __doc__ = "get Psite positions for each read"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--input_bed_file",
                        dest="input_bed_file",
                        help="Path to the input bed file, expect to have 'cigar' at 7th column",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--output_bed_file_plus",
                        dest="output_bed_file_plus",
                        help="Path to write the bed file with P-site positions, on + strand",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--output_bed_file_minus",
                        dest="output_bed_file_minus",
                        help="Path to write the bed file with P-site positions, on - strand",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--offset",
                        dest="offset",
                        help="estimated offset from read 5p end to the P-site",
                        required=False,
                        default=-1,
                        metavar="FILE",)

    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    input_bed_file_path = options.input_bed_file
    output_bed_file_plus = options.output_bed_file_plus
    output_bed_file_minus = options.output_bed_file_minus
    offset = int(options.offset)

    cigar_pattern = re.compile(r'(\d+)([MIDNSHPX=])')
    # we need to look only at cigar elements that "consume" reference
    reference_consuming_elements = ['M','D','N','=','X']
    offset_consuming_elements = ['M','=','X']

    if offset>=0 and os.path.isfile(input_bed_file_path) and os.path.getsize(input_bed_file_path)>0:
        bed = pd.read_csv(input_bed_file_path,delimiter="\t",index_col=None,header=None)
        bed = bed.rename(columns={5:'strand',6:'cigar'})
        cigars = bed[['strand','cigar']].drop_duplicates().reset_index(drop=True)
        
        def get_shift(x,offset,reference_consuming_elements,offset_consuming_elements):
            shift = 0
            consumed_nts = 0
            reference_shift = 0
            cigar_elements = cigar_pattern.findall(x['cigar'])
            
            if x['strand']=='-':
                cigar_elements = cigar_elements[::-1]
            for length_str, op_type in cigar_elements:
                length_str = int(length_str)
                if op_type in reference_consuming_elements:
                    if op_type in offset_consuming_elements:
                        if length_str+consumed_nts<offset:
                            # proceed
                            consumed_nts = consumed_nts+length_str
                        else:
                            # need to stop at this element
                            shift = reference_shift+(offset-consumed_nts)
                            break
                    reference_shift = reference_shift+length_str
            return shift
            
        cigars['shift'] = cigars.apply(lambda x:get_shift(x,offset,reference_consuming_elements,offset_consuming_elements),1)
        bed = pd.merge(bed,cigars,how='left',on=['strand','cigar'])
        tmp1 = bed.loc[bed['strand']=='+']
        tmp1['start'] = tmp1[1]+tmp1['shift']
        tmp1['end'] = tmp1['start']+1
        tmp1[[0,'start','end']].sort_values([0,'start','end']).to_csv(output_bed_file_plus, sep=str('\t'),header=False,index=None,quoting=csv.QUOTE_NONE)
        
        tmp2 = bed.loc[bed['strand']=='-']
        tmp2['end'] = tmp2[2]-tmp2['shift']
        tmp2['start'] = tmp2['end']-1
        tmp2[[0,'start','end']].sort_values([0,'start','end']).to_csv(output_bed_file_minus, sep=str('\t'),header=False,index=None,quoting=csv.QUOTE_NONE)
    else:
        print('[INFO] input file is empty or offset is not determined\n')
        # create empty files
        tmp = pd.DataFrame(columns=[0,'start','end'])
        tmp.to_csv(output_bed_file_plus, sep=str('\t'),header=False,index=None,quoting=csv.QUOTE_NONE)
        tmp.to_csv(output_bed_file_minus, sep=str('\t'),header=False,index=None,quoting=csv.QUOTE_NONE)
            
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)