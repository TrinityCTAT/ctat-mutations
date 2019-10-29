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
parser.add_argument("--CosmicCodingMuts", required = True ,help="CosmicCodingMut VCF file")
parser.add_argument("--CosmicMutantExport" ,required = True ,help="CosmicMutantExport TSV file")
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

# should close oe?

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
