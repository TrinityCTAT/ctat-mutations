#!/usr/bin/env python


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

parser = argparse.ArgumentParser()
parser.add_argument("--CosmicCodingMuts", required = True ,help="CosmicCodingMut VCF file")
parser.add_argument("--CosmicMutantExport" ,required = True ,help="CosmicMutantExport TSV file")
parser.add_argument("--outfile" ,required = True ,help="final results file")
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
with gzip.open(args.CosmicMutantExport,"r") as mt:
    mutant_reader=csv.DictReader(mt,delimiter="\t", quoting=csv.QUOTE_NONE)
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

print "read mutant"
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
print "read in vcf"
        
with gzip.open("only_entries_tmp.tsv.gz","w") as oe:
    lines=[header]+entries
    oe.writelines(lines)

print "made temp file with header and entries only"

with gzip.open("only_entries_tmp.tsv.gz","r") as oer:
    only_entries=csv.DictReader(oer,delimiter="\t",quoting=csv.QUOTE_NONE)
    print "entries read into dict"
    if not (args.outfile).endswith(".gz"):
        outfile=args.outfile+".gz"
    else:
        outfile=args.outfile
    with gzip.open(outfile,"w") as iw:
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
         print "Done."
