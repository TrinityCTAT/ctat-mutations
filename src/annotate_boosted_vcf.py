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
                                     description = "annotate boosted variants subset from vcf")

    parser.add_argument("--vcf_in", required = True, help="Input vcf.")
    parser.add_argument("--boosted_variants_matrix", required=True,
                        help="matrix indicating boosted variants. Requires columns: 'chr:pos' and 'boosted'",
                        nargs='+')
    parser.add_argument("--boost_type", required=True, help="boosting model type")
    parser.add_argument("--vcf_out", required=True, help="Output vcf.")
    
    args = parser.parse_args()

    vcf_input_file = args.vcf_in
    boosted_matrices = args.boosted_variants_matrix # can provide one for the indels and one for the snps, and we combine them.
    vcf_output_file = args.vcf_out
    boost_type = args.boost_type
    
    # get list of boosted variants
    logger.info("-parsing list of boosted variants from: {}".format(vcf_input_file))
    boosted_variants = set()
    for boosted_matrix in boosted_matrices:
        with open(boosted_matrix) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row['boosted'].lower() == 'true':
                    boosted_variants.add(row['chr:pos'])

    logger.info("-identified {} boosted variants".format(len(boosted_variants)))
    
    # output vcf

    if re.search(".gz$", vcf_input_file):
        fh = gzip.open(vcf_input_file, 'rt', encoding='utf-8')
    else:
        fh = open(vcf_file, "rt", encoding='utf-8')

    num_reported_variants = 0

    
    with open(vcf_output_file, 'wt', encoding='utf-8') as ofh:
        for line in fh:
            if line[0] == "#":
                if re.match("#CHROM\t", line):
                    print('##INFO=<ID=BOOSTselect,Number=A,Type=String,Description="selected according to boosting method">', file=ofh)
                ofh.write(line)
                continue

            vals = line.split("\t")
            chrpos = ":".join(vals[0:2])

            if chrpos in boosted_variants:
                # include variant annotation for boosting type
                vals = line.split("\t")
                vals[7] += f";BOOSTselect={boost_type}"
                line = "\t".join(vals)
                num_reported_variants += 1

            ofh.write(line)

    if num_reported_variants < len(boosted_variants):
        raise RuntimeError("Error, only {} of {} variants were extracted from the vcf".format(num_reported_variants, len(boosted_variants)))

    
    sys.exit(0)

    
if __name__=='__main__':
    main()


    

         
