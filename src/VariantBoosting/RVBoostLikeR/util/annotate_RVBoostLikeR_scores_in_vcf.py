#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os, sys, re
import logging
import argparse
import subprocess

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    parser = argparse.ArgumentParser(description="add rvboost score annotations",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--input_vcf", type=str,
                        required=True, help="input vcf file")

    parser.add_argument("--rvboost_outdir", type=str,
                        required=True, help="rvboost output directory")

    parser.add_argument("--output_vcf", type=str,
                        required=True, help="output vcf file including score annotations")
    
    args = parser.parse_args()

    
    ## get orig score annotations.

    orig_score_file = os.path.join(args.rvboost_outdir,
                                   "original_score.txt")

    if not os.path.exists(orig_score_file):
        logger.critical("Error, cannot find file: {}".format(orig_score_file))
        sys.exit(1)
    
    orig_scores_vec = subprocess.check_output("cat {}".format(orig_score_file), shell=True, encoding='utf8').split("\n")

    ## get rv qscore annotations

    rv_qscore_file = os.path.join(args.rvboost_outdir,
                                  "RV.Qscore.txt")

    if not os.path.exists(rv_qscore_file):
        logger.critical("Error, cannot find file: {}".format(rv_qscore_file))
        sys.exit(1)


    rv_qscores_vec = subprocess.check_output("cat {}".format(rv_qscore_file), shell=True, encoding='utf8').split("\n")
    

    ## write annotated vcf
    counter = -1

    ofh = open(args.output_vcf, 'w')
    with open(args.input_vcf) as fh:

        for line in fh:
            line = line.rstrip()
            if line[0] == "#":
                # header line
                
                if re.match("#CHROM\tPOS", line):
                    print("##INFO=<ID=OrgScore,Number=1,Type=Float,Description=\"RV Boosting algorithm Original score\">\n" +
                          "##INFO=<ID=QScore,Number=1,Type=Float,Description=\"RV Boosting algorithm Q-Score\">",
                          file=ofh)
                
                
                print(line, file=ofh)
                continue

            counter += 1
            line = line.rstrip()
            vals = line.split("\t")
            
            vals[7] += ";OrgScore={:0.3f};QScore={:0.3f}".format(float(orig_scores_vec[counter]),
                                                                 float(rv_qscores_vec[counter]))

            print("\t".join(vals), file=ofh)


    print("Done.", file=sys.stderr)

    sys.exit(0)
    
    
    
    
        
 
####################
 
if __name__ == "__main__":
    main()
