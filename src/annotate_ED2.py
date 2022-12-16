#!/usr/bin/env python
# coding: utf-8

# ------------------------------
# Import the needed libraries
# ------------------------------

import argparse
import csv
import gzip
import logging
import sys, os, re
import subprocess
import csv
import pandas as pd
from collections import defaultdict

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("ctat_boosting")


def main():

    # add options to inputs
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter, description="Adds ED.\n"
    )

    parser.add_argument("--input_vcf", required=True, help="input vcf file")
    parser.add_argument(
        "--output_vcf",
        required=True,
        help="output vcf file including annotation for distance to splice neighbor",
    )
    parser.add_argument("--reference", required=True, help="Reference Genome.")
    parser.add_argument("--temp_dir", default="/tmp", help="tmp directory")
    parser.add_argument(
        "--threads",
        default=1,
        type=int,
        required=False,
        help="number of threads to use by pblat",
    )

    args = parser.parse_args()

    input_vcf = args.input_vcf
    out_file = args.output_vcf
    reference = args.reference
    temp_dir = args.temp_dir
    threads = args.threads

    logger.info(
        "\n################################\n Annotating VCF: Calculating ED \n################################\n"
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Process the VCF input
    # ~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("Processing VCF Positions")
    input_vcf_df = pd.read_csv(
        input_vcf,
        sep="\t",
        low_memory=False,
        comment="#",
        header=None,
        names=[
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "ENCODING",
        ],
    )

    temp = input_vcf_df.loc[:, ["CHROM", "POS"]]
    temp["START"] = temp["POS"].apply(lambda x: max(1, x - 50))
    temp["END"] = temp["POS"] + 50
    positions = (
        temp["CHROM"] + ":" + temp["START"].astype(str) + "-" + temp["END"].astype(str)
    )

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

    del temp # no longer needed, free mem
    
    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Run SAMTOOLS - FAIDX
    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Run samtools FAIDX to get sequences
    logger.info("Running samtools faidx")
    cmd = "samtools faidx {} --region-file {} > {}".format(
        reference, positions_file, faidx_output
    )
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    # fsa = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', shell=True).communicate()[0]

    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Run BLAT
    # ~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("Running Blat")

    # check if temp exists, if not make it
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    # Path to the psl file
    psl_output = os.path.join(temp_dir, "blat_output.psl")

    cmd = "pblat {} {} -threads={} -noHead -minScore=70 -minIdentity=90 {}".format(
        reference, faidx_output, threads, psl_output
    )
    print(cmd)
    # subprocess.Popen(cmd, stdin=subprocess.PIPE, encoding='utf8', shell=True).communicate(fsa)
    if not os.path.exists("blat_output.psl.ok"):
        subprocess.check_call(cmd, shell=True)

    subprocess.check_call("touch blat_output.psl.ok", shell=True)
    
    # Process the psl output
    logger.info("Processing Output")
    header = [
        "match",
        "mis-match",
        "rep. match",
        "N's",
        "Q gap count",
        "Q gap bases",
        "T gap count",
        "T gap bases",
        "strand",
        "Q name",
        "Q size",
        "Q start",
        "Q end",
        "T name",
        "T size",
        "T start",
        "T end",
        "block count",
        "blockSizes",
        "qStarts",
        "tStarts",
    ]
    psl_df = pd.read_csv(
        psl_output, sep="\t", low_memory=False, comment="#", names=header
    )

    psl_fh = open(psl_output, "rt")
    reader = csv.DictReader(psl_fh, delimiter="\t", fieldnames=header)
    

    # ~~~~~~~~~~~~~~
    # Find ED Value
    # ~~~~~~~~~~~~~~
    # Create a dictionary to calculate the ED values
    logger.info("Creating ED features")
    hit_counts = defaultdict(int)
    for row in reader:
        names = re.split(":|-", row["Q name"])
        val1 = abs(int(names[1]) - int(row["T start"]) )
        val2 = abs(int(names[2]) - int(row["T end"]) )
        hit_counts[row["Q name"]] += 1

        
    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Process the VCF Header
    # ~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("Outputing the annotated VCF.")
    header = []
    with gzip.open(input_vcf, "rt") if input_vcf.endswith(".gz") else open(
        input_vcf, "rt"
    ) as vcf:
        for line in vcf:
            # check to see if header, is so append the line to header
            if re.match("##", line):
                header.append(line)
            else:
                break
                
    ## add the new ED annotation info
    header.append(
        '##INFO=<ID=ED,Number=1,Type=Integer,Description="Number of blat hits to reference genome, not counting self-hit">\n'
    )
    ## Join each line and write to a file
    vcf_header = "".join(header)
    add_header = open(out_file, "w")
    add_header.write(vcf_header)
    add_header.close()

    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Write variants and add annotations
    # ~~~~~~~~~~~~~~~~~~~~~~~
    # write the variants to the new VCF

    vcf_fh = gzip.open(input_vcf, "rt") if input_vcf.endswith(".gz") else open(input_vcf, "rt")
    vcf_reader = csv.DictReader(filter(lambda row: row[0]!='#' or re.match("#CHROM", row), vcf_fh), delimiter="\t")
    fieldnames = list(vcf_reader.fieldnames)
    vcf_ofh = open(out_file, "at")
    vcf_writer = csv.DictWriter(vcf_ofh, fieldnames=fieldnames, delimiter="\t", lineterminator='\n')
    vcf_writer.writeheader()
    pos_index = -1
    for row in vcf_reader:
        pos_index += 1
        pos_name = positions[pos_index]
        if pos_name in hit_counts:
            row["INFO"] = row["INFO"] + ";ED={}".format(
                hit_counts[pos_name]
            )
        else:
            row["INFO"] = row["INFO"] + ";ED=-1"

        vcf_writer.writerow(row)


    sys.exit(0)
    

if __name__ == "__main__":

    main()
