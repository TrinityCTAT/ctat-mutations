#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


#import inspect
import os,sys
import csv 
import argparse
import subprocess
import gzip
import glob
import logging


FORMAT = "%(asctime)-15s: %(message)s"
logger = logging.getLogger()
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)


parser = argparse.ArgumentParser()
parser.add_argument("--genome_lib_dir" ,dest="genome_lib_dir", type=str,
                                        default=os.environ.get('CTAT_GENOME_LIB'),
                                        help="genome lib directory - see http://FusionFilter.github.io for details. Uses env var CTAT_GENOME_LIB as default")
args=parser.parse_args()

csv.field_size_limit(sys.maxsize)

#Check for Picard env
picard_path=os.environ.get("PICARD_HOME")

if not picard_path:
    sys.exit("Set $PICARD_HOME to your picard software")


#Generating ref_genome.dict in ctat_genome_lib
#if not present
ref_dict=os.path.join(genome_lib_dir,"ref_genome.dict")
if not os.path.exists(ref_dict): 
    logger.info("Generating "+ ref_dict)
    picard_jar = os.path.join(picard_path,"picard.jar")
    ref_fa=os.path.join(genome_lib_dir,"ref_genome.fa")
    cmd = ["java", "-jar", picard_jar,
           "CreateSequenceDictionary",
           "R=",ref_fa,
           "O=",ref_dict,
           "VALIDATION_STRINGENCY=LENIENT"]
    subprocess.check_call(cmd)

if os.path.exists(ref_dict):
    logger.info("ref dict created for gatk use in pipe")
else:
    raise RuntimeError("Error, ref dict could not be generated with picard")


logger.info("Done prepping ctat mutation lib.")

sys.exit(0)
