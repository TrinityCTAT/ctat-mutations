#!/usr/bin/env python3

import sys, os, re
import csv
import gzip
import argparse
import logging


logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description = "extracts boosted variants subset from vcf")

    parser.add_argument("--vcf_in", required = True, help="Input vcf.")
    parser.add_argument("--boosted_variants_matrix", required=True,
                        help="matrix indicating boosted variants. Requires columns: 'chr:pos' and 'boosted'",
                        nargs='+')
    parser.add_argument("--vcf_out", required=True, help="Output vcf.")
    
    args = parser.parse_args()

    vcf_input_file = args.vcf_in
    boosted_matrices = args.boosted_variants_matrix # can provide one for the indels and one for the snps, and we combine them.
    vcf_output_file = args.vcf_out
    
    # get list of boosted variants
    logger.info("-parsing list of boosted variants from: {}".format(vcf_input_file))
    boosted_variants = set()
    for boosted_matrix in boosted_matrices:
        with open(boosted_matrix) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row['boosted'] == 'True':
                    boosted_variants.add(row['chr:pos'])

    logger.info("-identified {} boosted variants".format(len(boosted_variants)))
    
    # output vcf

    if re.search(".gz$", vcf_input_file):
        fh = gzip.open(vcf_input_file, 'rt')
    else:
        fh = open(vcf_file, "rt")

    num_reported_variants = 0
        
    with open(vcf_output_file, 'wt') as ofh:
        for line in fh:
            if line[0] == "#":
                ofh.write(line)
                continue

            vals = line.split("\t")
            chrpos = ":".join(vals[0:2])

            if chrpos in boosted_variants:
                ofh.write(line)
                num_reported_variants += 1
                continue

    if num_reported_variants < len(boosted_variants):
        raise RuntimeError("Error, only {} of {} variants were extracted from the vcf".format(num_reported_variants, len(boosted_variants)))

    
    sys.exit(0)

    
if __name__=='__main__':
    main()


    

         
