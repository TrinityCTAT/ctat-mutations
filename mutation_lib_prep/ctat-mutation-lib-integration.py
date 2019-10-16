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

UTILDIR = os.path.dirname(__file__)


parser = argparse.ArgumentParser()
parser.add_argument("--genome_lib_dir", dest="genome_lib_dir", type=str, required=True,
                    help="CTAT genome lib build dir")
args=parser.parse_args()


genome_lib_dir = args.genome_lib_dir

if not genome_lib_dir:
    sys.exit("Error, genome lib dir not set...  see usage")

csv.field_size_limit(sys.maxsize)

#Check for Picard env
picard_path=os.environ.get("PICARD_HOME")

if not picard_path:
    sys.exit("Set $PICARD_HOME to your picard software")



ctat_mutation_lib_dir = os.path.join(genome_lib_dir, "ctat_mutation_lib")
if not os.path.exists(ctat_mutation_lib_dir):
    sys.exit("Error, please be sure to unpack the ctat_mutation_lib.tar.gz in the ctat genome lib dir first")


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



# Generate splice adjacent targets

splice_adjacent_bed_file = os.path.join(ctat_mutation_lib_dir, "ref_annot.splice_adj.bed")

if not os.path.exists(splice_adjacent_bed_file + ".gz.tbi"):
    logger.info("extracting splice-adjacent regions")
    cmd = str(os.path.join(UTILDIR, "gencode_gtf_to_splice_adjacent_regions.pl") +
              " " + os.path.join(genome_lib_dir, "ref_annot.gtf") +
              " > {}".format(splice_adjacent_bed_file) )
    subprocess.check_call(cmd, shell=True)

    cmd = "bgzip {}".format(splice_adjacent_bed_file)
    subprocess.check_call(cmd, shell=True)

    cmd = "tabix -p bed {}.gz".format(splice_adjacent_bed_file)
    subprocess.check_call(cmd, shell=True)



refGene_sorted_bed_file = os.path.join(ctat_mutation_lib_dir, "refGene.sort.bed")

if not os.path.exists(refGene_sorted_bed_file + ".gz.tbi"):

    logger.info("prepping gene regions bed")
    
    refGene_bed_file = os.path.join(ctat_mutation_lib_dir, "refGene.bed")

    cmd = str(os.path.join(UTILDIR, "gencode_gtf_to_bed.pl") +
              " " + os.path.join(genome_lib_dir, "ref_annot.gtf") +
              " > {}".format(refGene_bed_file) )
    subprocess.check_call(cmd, shell=True)

    cmd = "sort -k 1,1 -k2,2g -k3,3g < {} > {}".format(refGene_bed_file, refGene_sorted_bed_file)
    subprocess.check_call(cmd, shell=True)

    cmd = "bgzip {}".format(refGene_sorted_bed_file)
    subprocess.check_call(cmd, shell=True)

    cmd = "tabix -p bed {}.gz".format(refGene_sorted_bed_file)
    subprocess.check_call(cmd, shell=True)

    

logger.info("Done prepping ctat mutation lib.")

sys.exit(0)
