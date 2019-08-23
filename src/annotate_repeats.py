#!/usr/bin/env python


import argparse
import sys, os
import subprocess
from collections import defaultdict
import logging
import re

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../PyLib"]))

import ctat_util

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    
    #add options to inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Adds repeat feature annotations to vcf file.\n")
    
    parser.add_argument('--input_vcf', required=True, help="input vcf file")
    parser.add_argument('--repeats_bed', required=True, help='repeat features bed file')
    parser.add_argument('--output_vcf', required=True, help="output vcf fiel including repeat feature annotations")

    parser.add_argument('--debug', default=False, action='store_true', help='debug mode, retains temporary intermediate files')


    args = parser.parse_args()
    


    input_vcf_file = args.input_vcf
    repeats_bed_file = args.repeats_bed
    output_vcf_file = args.output_vcf
    DEBUG_MODE = args.debug
    

    ## identify repeat feature overlaps:
    rpt_intersect_file = "{}.rpt_intersect".format(output_vcf_file)
    cmd = "bedtools intersect -wb -a {} -b {} | cut -f1,2,14 > {}".format(input_vcf_file,
                                                                          repeats_bed_file,
                                                                          rpt_intersect_file)

    logger.info("CMD: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)

    rpt_dict = defaultdict(set)
    with open(rpt_intersect_file) as fh:
        for line in fh:
            line = line.rstrip()
            (chr, pos, rpt_type) = line.split("\t")
            chrpos = "{}:{}".format(chr, pos)
            rpt_dict[chrpos].add(rpt_type)
    
    
    logger.info("Adding repeat annotations.")
    ## make output a vcf formatted file:
    with open(output_vcf_file, 'w') as ofh:
        
        with ctat_util.open_file_for_reading(input_vcf_file) as fh:
            for line in fh:
                if line[0] == "#":
                    
                    if re.match("#CHROM\t", line):
                        # add header info line for the repeat annotation type
                        ofh.write("##INFO=<ID=RPT,Type=String,Description=\"Repeat family from UCSC Genome Browser Repeatmasker Annotations\">\n")

                    ofh.write(line)
                else:
                    line = line.rstrip()
                    vals = line.split("\t")
                    chrpos = "{}:{}".format(vals[0], vals[1])
                    if chrpos in rpt_dict:
                        rpt_type = ",".join(sorted(rpt_dict[chrpos]))
                        vals[7] += ";RPT={}".format(rpt_type)
                    ofh.write("\t".join(vals) + "\n")


    # cleanup
    if not DEBUG_MODE:
        os.remove(rpt_intersect_file)
    
    sys.exit(0)
    

if __name__ == "__main__":

    main()
