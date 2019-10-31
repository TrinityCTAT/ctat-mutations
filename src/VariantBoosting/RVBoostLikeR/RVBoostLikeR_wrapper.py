#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os, sys, re

if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)

import datetime
import subprocess
import argparse
import logging
FORMAT = "%(asctime)-15s: %(levelname)s %(module)s.%(name)s.%(funcName)s %(message)s"
logger = logging.getLogger('ctat_mutations')
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)


sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"]))
from Pipeliner import Pipeliner, Command, run_cmd, ParallelCommandList

UTILDIR = os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "util"])

def main():

    parser = argparse.ArgumentParser(description="wrapper for running rvboost", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--input_vcf", type=str, required=True, help="input vcf file")

    parser.add_argument("--output_dir", type=str, required=True, help="output directory name")


    parser.add_argument("--attributes", type=str, required=False, help="vcf info attributes to use for scoring purposes",
                        default="DJ,PctExtPos,ReadPosRankSum,QD,FS,ED")
    

    parser.add_argument("--score_threshold", type=float, required=False, default=0.05, help="score threshold for filtering rvboost results")

    args = parser.parse_args()
    
    
    ## prep for run
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    checkpts_dir = os.path.join(output_dir, "__chckpts")

    pipeliner = Pipeliner(checkpts_dir)
    
    ## build pipeline

    # run rvboost
    rvboost_score_dir = os.path.join(args.output_dir, "rvboost_dir")
    cmd = " ".join([ os.path.join(UTILDIR, "my.RVboost.R"),
                     args.input_vcf,
                     rvboost_score_dir,
                     args.attributes ])
    
    pipeliner.add_commands([Command(cmd, "rvboost_core.ok")])
    
    # incorporate score annotations into vcf file:
    vcf_w_rvb_score_filename = os.path.join(args.output_dir, os.path.splitext(os.path.basename(args.input_vcf))[0] + ".RVB.vcf")
    cmd = " ".join([ os.path.join(UTILDIR, "my.RVboost.add_score_annotations.py"),
                     "--input_vcf", args.input_vcf,
                     "--rvboost_outdir", rvboost_score_dir,
                     "--output_vcf", vcf_w_rvb_score_filename ])
    
    pipeliner.add_commands([Command(cmd, "annot_rvb_score.ok")])
    

    # apply a filter on the score threshold:

    score_filtered_vcf = os.path.join(args.output_dir, os.path.splitext(os.path.basename(vcf_w_rvb_score_filename))[0] + ".filt_Q{}.vcf".format(args.score_threshold))
    
    cmd = " ".join([ os.path.join(UTILDIR, "filter_rvboost_by_Qscore.pl"),
                     vcf_w_rvb_score_filename,
                     str(args.score_threshold),
                     " > {} ".format(score_filtered_vcf) ])

    pipeliner.add_commands([Command(cmd, "filt_qscore_{}.ok".format(args.score_threshold))])
    
    
    
    pipeliner.run()

    sys.exit(0)
    


if __name__ == '__main__':
    main()
    
