#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import argparse
import logging
import os
import subprocess

import sys


FORMAT = "%(asctime)-15s: %(message)s"
logger = logging.getLogger()
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)

UTILDIR = os.path.dirname(__file__)
parser = argparse.ArgumentParser()
parser.add_argument("--genome_lib_dir", dest="genome_lib_dir", type=str,
    default=os.environ.get("CTAT_GENOME_LIB"),
    help="CTAT genome lib dir")
parser.add_argument("--gatk", dest="gatk", type=str, help="Optional path to GATK", default='gatk')
args = parser.parse_args()

gatk = args.gatk
ctat_genome_lib_path = args.genome_lib_dir

if ctat_genome_lib_path is None or not os.path.exists(ctat_genome_lib_path):
    exit("Missing path to CTAT_GENOME_LIB in $CTAT_GENOME_LIB.")


def is_in_path(path):
    try:
        if (
                subprocess.run(
                    ["which", path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT
                ).returncode
                != 0
        ):
            return False
    except:
        return False
    return True


if not is_in_path('gatk'):
    sys.exit("GATK not found")

ctat_mutation_lib_dir = os.path.join(ctat_genome_lib_path, "ctat_mutation_lib")
if not os.path.exists(ctat_mutation_lib_dir):
    sys.exit("Error, please be sure to unpack the ctat_mutation_lib.tar.gz in the ctat genome lib dir first")

# Generating ref_genome.dict in ctat_genome_lib
# if not present
ref_dict = os.path.join(ctat_genome_lib_path, "ref_genome.dict")
if not os.path.exists(ref_dict):
    logger.info("Generating " + ref_dict)
    ref_fa = os.path.join(ctat_genome_lib_path, "ref_genome.fa")
    cmd = [gatk,
           "CreateSequenceDictionary",
           "R=", ref_fa,
           "O=", ref_dict,
           "VALIDATION_STRINGENCY=LENIENT"]
    subprocess.check_call(cmd)

if os.path.exists(ref_dict):
    logger.info("ref dict created for gatk use in pipe")
else:
    raise RuntimeError("Error, ref dict could not be generated")

# Generate splice adjacent targets

splice_adjacent_bed_file = os.path.join(ctat_mutation_lib_dir, "ref_annot.splice_adj.bed")

if not os.path.exists(splice_adjacent_bed_file + ".gz.tbi"):
    logger.info("extracting splice-adjacent regions")
    cmd = str(os.path.join(UTILDIR, "gencode_gtf_to_splice_adjacent_regions.pl") +
              " " + os.path.join(ctat_genome_lib_path, "ref_annot.gtf") +
              " > {}".format(splice_adjacent_bed_file))
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
              " " + os.path.join(ctat_genome_lib_path, "ref_annot.gtf") +
              " > {}".format(refGene_bed_file))
    subprocess.check_call(cmd, shell=True)

    cmd = "sort -k 1,1 -k2,2g -k3,3g < {} > {}".format(refGene_bed_file, refGene_sorted_bed_file)
    subprocess.check_call(cmd, shell=True)

    cmd = "bgzip {}".format(refGene_sorted_bed_file)
    subprocess.check_call(cmd, shell=True)

    cmd = "tabix -p bed {}.gz".format(refGene_sorted_bed_file)
    subprocess.check_call(cmd, shell=True)

logger.info("Done prepping ctat mutation lib.")

sys.exit(0)
