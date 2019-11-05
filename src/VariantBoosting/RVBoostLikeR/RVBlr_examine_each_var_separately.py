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

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

def main():

    parser = argparse.ArgumentParser(description="wrapper for running rvboost", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--input_vcf", type=str, required=True, help="input vcf file")

    parser.add_argument("--output_dir", type=str, required=True, help="output directory name")

    parser.add_argument("--attributes", type=str, required=False, help="vcf info attributes to use for scoring purposes",
                        default="DJ,PctExtPos,ReadPosRankSum,QD,FS,ED")
    
    parser.add_argument("--score_threshold", type=float, required=False, default=0.05, help="score threshold for filtering rvboost results")

    parser.add_argument("--num_threads", type=int, required=False, default=4, help="number of concurrent processes")

    args = parser.parse_args()
    
    
    ## prep for run
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    checkpts_dir = os.path.join(output_dir, "__chckpts")

    pipeliner = Pipeliner(checkpts_dir)
    
    ## build pipeline

    atts_list = args.attributes.split(",")

    cmd_list = []
    for att in atts_list:
        cmd = " ".join([ os.path.join(SCRIPTDIR, "RVBoostLikeR_wrapper.py"),
                         "--input_vcf", args.input_vcf,
                         "--attributes", att,
                         "--score_threshold {} ".format(args.score_threshold),
                         "--work_dir {}/rvb-att-{}".format(args.output_dir, att),
                         "--output_filename {}/rvb-att-{}.thresh{:.3f}.vcf".format(args.output_dir, att, args.score_threshold) ])
        cmd_list.append(cmd)
    
    
    pipeliner.add_commands([ParallelCommandList(cmd_list, "all_atts-rvb.ok", args.num_threads)])
    
    
    pipeliner.run()

    sys.exit(0)
    


if __name__ == '__main__':
    main()
    
