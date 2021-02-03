#!/usr/bin/env python
# coding: utf-8

#------------------------------
# Import the needed libraries
#------------------------------

import argparse
import csv
import gzip
import logging
import os
import re
import subprocess

import pandas as pd


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger('ctat_boosting')


def main():

    #add options to inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Adds ED.\n")

    parser.add_argument('--input_vcf', required=True, help="input vcf file")
    parser.add_argument('--output_vcf', required=True, help="output vcf file including annotation for distance to splice neighbor")
    parser.add_argument('--reference', required=True, help="Reference Genome.")
    parser.add_argument("--temp_dir", default="/tmp", help="tmp directory")


    args = parser.parse_args()



    input_vcf = args.input_vcf
    out_file = args.output_vcf
    reference = args.reference
    temp_dir = args.temp_dir



    logger.info("\n################################\n Annotating VCF: Calculating ED \n################################\n")


    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Process the VCF input
    # ~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("Processing VCF Positions")
    input_vcf_df = pd.read_csv(input_vcf,
                         sep='\t', low_memory=False, comment='#', header =None,
                         names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "ENCODING"])

    temp = input_vcf_df.loc[:,["CHROM", "POS"]]
    temp["START"] = temp["POS"] - 50
    temp["END"] = temp["POS"] + 50
    positions = temp["CHROM"] + ":" + temp["START"].astype(str) + "-" + temp["END"].astype(str)

    # create the FA file (samtools faidx cant take in stdin)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    if os.path.exists(temp_dir):
        positions_file = os.path.join(temp_dir, "positions.fa")
        faidx_output = os.path.join(temp_dir, "faidx_output.fa")
    # Write the positions
    outfile = open(positions_file, "w")
    outfile.write("\n".join(positions))
    outfile.close()


    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Run SAMTOOLS - FAIDX
    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Run samtools FAIDX to get sequences
    logger.info("Running samtools faidx")
    cmd = "samtools faidx {} --region-file {} > {}".format(reference, positions_file, faidx_output)
    print(cmd)
    fsa = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', shell=True).communicate()[0]

    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Run BLAT
    # ~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("Running Blat")

    # check if temp exists, if not make it
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    # Path to the psl file
    psl_output = os.path.join(temp_dir, "blat_output.psl")

    cmd = "pblat {} {} -threads=4 -noHead -minScore=70 -minIdentity=90 {}".format(reference, faidx_output, psl_output)
    print(cmd)
    subprocess.Popen(cmd, stdin=subprocess.PIPE, encoding='utf8', shell=True).communicate(fsa)

    # Process the psl output
    logger.info("Processing Output")
    header = ["match", "mis-match", "rep. match", "N's", "Q gap count", "Q gap bases", "T gap count", "T gap bases", "strand", "Q name", "Q size", "Q start", "Q end", "T name", "T size", "T start", "T end", "block count", "blockSizes", "qStarts", "tStarts"]
    psl_df = pd.read_csv(psl_output,
                         sep='\t', low_memory=False, comment='#', names = header)

    # ~~~~~~~~~~~~~~
    # Find ED Value
    # ~~~~~~~~~~~~~~
    # Create a dictionary to calculate the ED values
    logger.info("Creating ED features")
    dic = {}
    names = [re.split(":|-", i) for i in psl_df["Q name"]]
    for i in range(len(names)):
        val1 = abs(int(names[i][1]) - psl_df["T start"][i])
        val2 = abs(int(names[i][2]) - psl_df["T end"][i])
        if val1 >= 1 or val2 >= 1:
            if psl_df["Q name"][i] not in dic:
                dic[psl_df["Q name"][i]] = 1
            else:
                dic[psl_df["Q name"][i]] += 1
    ## Add the ED to each variant
    for i in range(len(positions)):
        if positions[i] in dic:
            input_vcf_df.loc[i,["INFO"]] = input_vcf_df["INFO"][i] + ";ED={}".format(dic[positions[i]])
        else:
            input_vcf_df.loc[i,["INFO"]] = input_vcf_df["INFO"][i] + ";ED=-1"


    # Remove the output file if it exist
    if os.path.exists(out_file):
        os.remove(out_file)

    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Process the VCF Header
    # ~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("Outputing the annotated VCF.")
    header = []
    with gzip.open(input_vcf, 'rt') if input_vcf.endswith('.gz') else open(input_vcf, "rt") as vcf:
        for line in vcf.readlines():
            # check to see if header, is so append the line to header
            if line[0] == "#":
                header.append(line)
                continue
    ## add the new DJ annotation info
    header.insert(-1,"##INFO=<ID=ED,Number=1,Type=Integer,Description=\"Number of blat hits to reference genome, not counting self-hit\">\n")
    ## Join each line and write to a file
    vcf_header = "".join(header)
    add_header = open(out_file,"w")
    add_header.write(vcf_header)
    add_header.close()

    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Write variants
    # ~~~~~~~~~~~~~~~~~~~~~~~
    # write the variants to the new VCF
    input_vcf_df.to_csv(out_file, mode='a', header=False, index = False, sep = "\t", quoting=csv.QUOTE_NONE)


if __name__ == "__main__":

    main()



