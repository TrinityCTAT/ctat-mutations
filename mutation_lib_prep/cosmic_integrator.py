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



##
## This script decorates the Cosmic coding variants with cancer census annotations.

FORMAT = "%(asctime)-15s: %(message)s"
logger = logging.getLogger()
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)


parser = argparse.ArgumentParser()
parser.add_argument("--CosmicCodingMuts", required = True ,help="CosmicCodingMut VCF file")
parser.add_argument("--CosmicMutantExport", required = True ,help="CosmicMutantExport TSV file")
parser.add_argument("--output_vcf", required=True, help="output vcf file")

args=parser.parse_args()

csv.field_size_limit(sys.maxsize)

##Add lines to header
add_header_lines = [
'##INFO=<ID=COSMIC_ID,Type=String,Description="COSMIC mutation id (unique).">\n',
'##INFO=<ID=TISSUE,Type=String,Description="The primary tissue/cancer and subtype from which the sample originated.">\n',
'##INFO=<ID=TUMOR,Type=String,Description="The histological classification of the sample.">\n',
'##INFO=<ID=FATHMM,Type=String,Description="FATHMM (Functional Analysis through Hidden Markov Models). \'Pathogenic\'=Cancer or damaging, \'Neutral\'=Passanger or Tolerated.">\n',
'##INFO=<ID=SOMATIC,Type=String,Description="Information on whether the sample was reported to be Confirmed Somatic. \'Confirmed somatic\'=if the mutation has been confimed to be somatic in the experiment by sequencing both the tumour and a matched normal from the same patient, \'Previously Observed\'=when the mutation has been reported as somatic previously but not in current paper, \'variant of unknown origin\'=when the mutation is known to be somatic but the tumour was sequenced without a matched normal">\n',
'##INFO=<ID=PUBMED_COSMIC,Type=String,Description="The PUBMED ID for the paper that the sample was noted in COSMIC.">\n'
]



####################################
# parsing the cancer gene census: CosmicMutantExport

#GENE,STRAND,CDS,AA,CNT
#COSMIC_ID,TISSUE,TUMOR,FATHMM,SOMATIC,PUBMED_COSMIC,GENE,STRAND,GENE,STRAND,CDS,AA,CNT
mutant_dict_necessary_info={}
logger.info("Capturing info from: {}".format(args.CosmicMutantExport))
with gzip.open(args.CosmicMutantExport,"rt") as mt:
    mutant_reader=csv.DictReader(mt, delimiter=str("\t"), quoting=csv.QUOTE_NONE)
    for row in mutant_reader:
        info_items=["COSMIC_ID="+row.get("GENOMIC_MUTATION_ID",""),
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
        mutant_dict_necessary_info[row["GENOMIC_MUTATION_ID"]]=info




logger.info("Now annotating {}".format(args.CosmicCodingMuts))

coding_muts_gzip_fh = gzip.open(args.CosmicCodingMuts,"rt")

cosmic_vcf=os.path.join(args.output_vcf)
logger.info("writing summary file: {}".format(cosmic_vcf))
ofh = open(cosmic_vcf, 'wt')

annotated_set = set()
not_annotated = set()

with gzip.open(args.CosmicCodingMuts,"rt") as fh:
    for line in fh:
        if line.startswith("##"):
            ofh.write(line)
            continue
        if line.startswith("#CHROM"):
            ofh.write("".join(add_header_lines))
            ofh.write(line)
            continue
        
        line = line.rstrip()
        vals = line.split("\t")
        vals[0] = "chr" + vals[0]
        cosmic_id = vals[2]
        if cosmic_id in mutant_dict_necessary_info:
            current_info = vals[7]
            vals[7] += ";" + mutant_dict_necessary_info[cosmic_id]
            annotated_set.add(cosmic_id)
        else:
            not_annotated.add(cosmic_id)

        ofh.write("\t".join(vals) + "\n")

ofh.close()

logger.info("-number of variants with annotations added: {}".format(len(annotated_set)))
logger.info("-number of variants w/o added annotations: {}".format(len(not_annotated)))


logger.info("bgzip compressing {}".format(cosmic_vcf))
subprocess.check_call("bgzip -f {}".format(cosmic_vcf), shell=True)

logger.info("indexing {}".format(cosmic_vcf))
subprocess.check_call(["bcftools", "index", "{}.gz".format(cosmic_vcf)])

logger.info("Done prepping cosmic vcf: {}".format(cosmic_vcf))

sys.exit(0)
