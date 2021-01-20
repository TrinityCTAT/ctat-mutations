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


sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../../PyLib"]))
from Pipeliner import Pipeliner, Command, run_cmd, ParallelCommandList

RVB_UTILDIR = os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "util"])
CTAT_UTILDIR = os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../.."])


def main():

    parser = argparse.ArgumentParser(description="wrapper for running rvboost-like-R", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--input_vcf", type=str, required=True, help="input vcf file")

    parser.add_argument("--work_dir", type=str, required=True, help="working directory name for intermediates")

    parser.add_argument("--attributes", type=str, required=False, help="vcf info attributes to use for scoring purposes",
                        default="DJ,ReadPosRankSum,QD,FS,ED,PctExtPos,RS")
    

    parser.add_argument("--score_threshold", type=float, required=False, default=0.05, help="score threshold for filtering rvboost results")


    parser.add_argument("--output_filename", required=True, help="name of output file containing final list of variants")
    
    args = parser.parse_args()
    
    
    ## prep for run
    output_dir = args.work_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    checkpts_dir = os.path.join(output_dir, "__chckpts")

    pipeliner = Pipeliner(checkpts_dir)
    
    ## build pipeline

    # extract feature matrix
    matrix_file = os.path.join(output_dir, "features.matrix")

    # ensure we capture the common variant annots
    if "RS" not in args.attributes.split(","):
        args.attributes += ",RS"
    
    
    cmd = " ".join([ os.path.join(RVB_UTILDIR, "extract_attribute_annotation_matrix.pl"),
                    args.input_vcf,
                    args.attributes,
                    ">",
                    matrix_file])
    pipeliner.add_commands([Command(cmd, "feature_extraction_to_matrix.ok")])
    
    
    # run rvboost-like-R
    boost_scores_file = os.path.join(output_dir, "RVB_like_R.var_scores")
    cmd = " ".join([ os.path.join(RVB_UTILDIR, "RVBoostLike.R"),
                    matrix_file,
                    boost_scores_file,
                    args.attributes ])
    
    pipeliner.add_commands([Command(cmd, "rvboost_core.ok")])
    
    # incorporate score annotations into vcf file:
    vcf_w_rvb_score_filename = os.path.join(output_dir, os.path.splitext(os.path.basename(args.input_vcf))[0] + ".RVBLR.vcf")
    cmd = " ".join([ os.path.join(RVB_UTILDIR, "annotate_RVBoostLikeR_scores_in_vcf.py"),
                     "--input_vcf", args.input_vcf,
                     "--scores", boost_scores_file,
                     "--output_vcf", vcf_w_rvb_score_filename ])
    
    pipeliner.add_commands([Command(cmd, "annot_rvb_score.ok")])
    

    # apply a filter on the score threshold:

    score_filtered_vcf = args.output_filename
    
    cmd = " ".join([ os.path.join(RVB_UTILDIR, "filter_by_RVBLRQ.pl"),
                     vcf_w_rvb_score_filename,
                     str(args.score_threshold),
                     " > {} ".format(score_filtered_vcf) ])

    pipeliner.add_commands([Command(cmd, "filt_qscore_{}.ok".format(args.score_threshold))])
    
    
    
    pipeliner.run()

    sys.exit(0)
    


if __name__ == '__main__':
    main()
    
