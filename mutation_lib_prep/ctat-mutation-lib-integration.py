#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


__author__ = "Asma Bankapur"
__copyright__ = "Copyright 2014"
__credits__ = [ "Vrushali Fangal", "Brian Haas" ]
__license__ = "MIT"
__maintainer__ = "Asma Bankapur"
__email__ = "bankapur@broadinstitute.org"
__status__ = "Development"

#import inspect
import os,sys
import csv 
import argparse
import subprocess
import gzip
import glob
import pandas as pd
import logging


FORMAT = "%(asctime)-15s: %(message)s"
logger = logging.getLogger()
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)


parser = argparse.ArgumentParser()
parser.add_argument("--CosmicCodingMuts", required = True ,help="CosmicCodingMut VCF file")
parser.add_argument("--CosmicMutantExport" ,required = True ,help="CosmicMutantExport TSV file")
parser.add_argument("--genome_lib_dir" ,dest="genome_lib_dir", type=str,
                                        default=os.environ.get('CTAT_GENOME_LIB'),
                                        help="genome lib directory - see http://FusionFilter.github.io for details. Uses env var CTAT_GENOME_LIB as default")
args=parser.parse_args()

csv.field_size_limit(sys.maxsize)

#Check for Picard env
picard_path=os.environ.get("PICARD_HOME")

if not picard_path:
    sys.exit("Set $PICARD_HOME to your picard software")

#Check for Genome lib
genome_lib_dir = args.genome_lib_dir
if genome_lib_dir is None:
    raise RuntimeError("Error, must specify --genome_lib_dir or set env var CTAT_GENOME_LIB")      
genome_lib_dir = os.path.abspath(genome_lib_dir)

#Search for compressed mutation lib
compressed_mutation_lib=glob.glob(os.path.join(genome_lib_dir,
                                     "mutation_lib.*.tar.gz"))
if len(compressed_mutation_lib) != 1:
    raise RuntimeError("Error, didn't locate single mutation_lib.*.tar.gz in {}".format(genome_lib_dir))
compressed_mutation_lib = compressed_mutation_lib[0]

#Uncompressing ctat_mutation_lib
ctat_mutation_lib=os.path.join(genome_lib_dir,"ctat_mutation_lib")
if not os.path.exists(ctat_mutation_lib):
    if compressed_mutation_lib:
        logger.info("Uncompressing mutation lib in genome_lib_dir")
        subprocess.check_call(["tar","-xzvf",compressed_mutation_lib,"-C",genome_lib_dir])
    else:
        raise RuntimeError("Cannot find mutation lib in CTAT_GENOME_LIB")

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

##Add lines to header
add_header_lines = [
'##INFO=<ID=COSMIC_ID,Type=String,Description="COSMIC mutation id (unique).">\n',
'##INFO=<ID=TISSUE,Type=String,Description="The primary tissue/cancer and subtype from which the sample originated.">\n',
'##INFO=<ID=TUMOR,Type=String,Description="The histological classification of the sample.">\n',
'##INFO=<ID=FATHMM,Type=String,Description="FATHMM (Functional Analysis through Hidden Markov Models). \'Pathogenic\'=Cancer or damaging, \'Neutral\'=Passanger or Tolerated.">\n',
'##INFO=<ID=SOMATIC,Type=String,Description="Information on whether the sample was reported to be Confirmed Somatic. \'Confirmed somatic\'=if the mutation has been confimed to be somatic in the experiment by sequencing both the tumour and a matched normal from the same patient, \'Previously Observed\'=when the mutation has been reported as somatic previously but not in current paper, \'variant of unknown origin\'=when the mutation is known to be somatic but the tumour was sequenced without a matched normal">\n',
'##INFO=<ID=PUBMED_COSMIC,Type=String,Description="The PUBMED ID for the paper that the sample was noted in COSMIC.">\n'
]

#GENE,STRAND,CDS,AA,CNT
#COSMIC_ID,TISSUE,TUMOR,FATHMM,SOMATIC,PUBMED_COSMIC,GENE,STRAND,GENE,STRAND,CDS,AA,CNT
mutant_dict_necessary_info={}
logger.info("Capturing info from: {}".format(args.CosmicMutantExport))
with gzip.open(args.CosmicMutantExport,"r") as mt:
    mutant_reader=csv.DictReader(mt, delimiter=str("\t"), quoting=csv.QUOTE_NONE)
    for row in mutant_reader:
        info_items=["COSMIC_ID="+row.get("Mutation ID",""),
                    "TISSUE="+row.get("Primary site",""),
                    "TUMOR="+row.get("Primary histology","")+" -- "+row.get("Histology subtype 1",""),
                    "FATHMM="+row.get("FATHMM prediction",""),
                    "SOMATIC="+row.get("Mutation somatic status",""),
                    "PUBMED_COSMIC="+row.get("Pubmed_PMID",""),
                    "GENE="+row.get("Gene name",""),
                    "STRAND="+row.get("Mutation strand",""),
                    "CDS="+row.get("Mutation CDS",""),
                    "AA="+row.get("Mutation AA","")]
        info=";".join(info_items)
        mutant_dict_necessary_info[row["Mutation ID"]]=info



logger.info("Capturing info from {}".format(args.CosmicCodingMuts))
comments=[]
entries=[]
with gzip.open(args.CosmicCodingMuts,"r") as cv:
    for row in cv:
        if row.startswith("##"):
            comments.append(row)   
        elif row.startswith("#"):
            header=row
        else:
            entries.append(row)

only_entries=os.path.join(ctat_mutation_lib,"only_entries_tmp.tsv.gz")
with gzip.open(only_entries,"w") as oe:
    lines=[header]+entries
    oe.writelines(lines)

logger.info("made temp file with header and entries only")

# write cosmic summary file
cosmic_vcf=os.path.join(ctat_mutation_lib,"cosmic.vcf")
logger.info("writing summary file: {}".format(cosmic_vcf))
with gzip.open(only_entries,"r") as oer:
    only_entries=csv.DictReader(oer,delimiter=str("\t"),quoting=csv.QUOTE_NONE)
    logger.debug("entries read into dict")
    with open(cosmic_vcf,"w") as iw:
         final_entries=[]
         for entry in only_entries:
             info=mutant_dict_necessary_info[entry["ID"]]
             CNT=entry["INFO"].split(";")[-1]
             full_info=";".join([info,CNT])
             if not entry["#CHROM"].startswith("chr"):
                 entry["#CHROM"]="chr"+entry["#CHROM"] 
             final_entries.append("\t".join([entry["#CHROM"],
                                                entry["POS"],
                                                entry["ID"],
                                                entry["REF"],
                                                entry["ALT"],
                                                entry["QUAL"],
                                                entry["FILTER"],
                                                full_info])+"\n")
         final_lines=comments+add_header_lines+[header]+final_entries
         iw.writelines(final_lines)  


logger.info("bgzip compressing {}".format(cosmic_vcf))
subprocess.check_call(["bgzip", cosmic_vcf])

logger.info("indexing {}".format(cosmic_vcf))
subprocess.check_call(["bcftools", "index", "{}.gz".format(cosmic_vcf)])

logger.info("Done prepping ctat mutation lib.")

sys.exit(0)
