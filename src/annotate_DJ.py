#!/usr/bin/env python
# coding: utf-8

#------------------------------
# Import the needed libraries
#------------------------------

import argparse
import csv
import logging
import os
import subprocess

import pandas as pd


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger('ctat_boosting')


def main():

    #add options to inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Adds exon splice distance annotations to vcf file (report up to len 10 distance away from splice).\n")

    parser.add_argument('--input_vcf', required=True, help="input vcf file")
    parser.add_argument('--splice_bed', required=True, help='Path to CTAT Mutations library splice adj bed file.')
    parser.add_argument('--output_vcf', required=True, help="output vcf file including annotation for distance to splice neighbor")
    parser.add_argument("--temp_dir", default="/tmp", help="tmp directory")

    args = parser.parse_args()

    input_vcf = args.input_vcf
    if input_vcf.endswith('.gz'):
        uncompressed_vcf = input_vcf[:len(input_vcf) - 3]
        if not os.path.exists(uncompressed_vcf):
            cmd = "gunzip -c {} > {}".format(input_vcf, uncompressed_vcf)
            logger.info("CMD: {}".format(cmd))
            subprocess.check_call(cmd, shell=True)
        input_vcf = uncompressed_vcf
        
    splice_bed = args.splice_bed
    out_file = args.output_vcf
    temp_dir = args.temp_dir

    # path to the temp sorted file
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    temp_sorted_vcf = os.path.join(temp_dir, "temp.sorted.vcf")


    logger.info("\n################################\n Annotating VCF: Calculating DJ \n################################\n")

    #~~~~~~~~~~~~~~
    # Sort VCF
    #~~~~~~~~~~~~~~
    # Sorting the VCF file by lexicographically
    ## have to do this for bedtools closest
    logger.info("Sorting VCF")
    cmd = "grep '^#' {} > {} && grep -v '^#' {} | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> {}".format(input_vcf, temp_sorted_vcf, input_vcf, temp_sorted_vcf)
    logger.info("CMD: {}".format(cmd))
    subprocess.run(cmd, shell=True, executable='/bin/bash')

    #~~~~~~~~~~~~~~
    # Load VCF
    #~~~~~~~~~~~~~~
    # Read in the input vcf as a data frame
    logger.info("Loading input VCF")
    input_vcf_df = pd.read_csv(temp_sorted_vcf,
                             sep='\t', low_memory=False, comment='#', header =None,
                             names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "ENCODING"])


    # Run BEDTools closeestBed
    logger.info("Running closestBed")
    cmd = "bedtools closest -header -t first -a {} -b {}".format(temp_sorted_vcf, splice_bed)
    logger.info("CMD: {}".format(cmd))
    distance_output = subprocess.check_output(cmd, shell=True).decode()


    # ~~~~~~~~~~~~~~~
    # Process Distances
    # ~~~~~~~~~~~~~~~
    # Convert the VCF from a string to a pandas dataframe
    temp = distance_output.split('\n')
    temp.remove('')
    variants = [x for x in temp if x[0] != "#"]
    test = pd.DataFrame(variants)
    if len(test) >0:
        ## split the one column into many
        vcf = test[0].str.split("\t",expand = True)

        # find the distances
        logger.info("Generating Distances")
        
        def compute_distance(row):
            var_pos = row[1]
            splice_adj_lend = row[11]
            splice_adj_rend = row[12]
            splice_adj_token = row[13]

            if splice_adj_token == ".":
                return ""
            
            splice_direction = splice_adj_token.split(":")[1]
            if splice_direction == "L":
                DJ = abs(int(var_pos) - (int(splice_adj_lend)-1))
                distance_string = ";DJ=" + str(DJ)
                return distance_string
            elif splice_direction == "R":
                DJ = abs(int(var_pos) - (int(splice_adj_rend)+1))
                distance_string = ";DJ=" + str(DJ)
                return distance_string
            else:
                raise RuntimeError("not recognizing {} as splice adj L or R from {}".format(splice_direction, splice_adj_token))

            
        distances = vcf.apply(compute_distance, axis=1)
        input_vcf_df["INFO"] = vcf[7] + distances

    # Remove the output file if it exist
    if os.path.exists(out_file):
      os.remove(out_file)


    #~~~~~~~~~~~~~~~~~~~
    # output results
    #~~~~~~~~~~~~~~~~~~~

    # Configure the vcfheader
    ## Parse out the vcf header lines
    header = [x for x in temp if x[0] == "#"]
    ## add the new DJ annotation info
    header.insert(-1,"##INFO=<ID=DJ,Number=1,Type=Integer,Description=\"Distance to closest junction\">" )
    ## Join each line and write to a file
    vcf_header = "\n".join(header) + "\n"
    add_header = open(out_file,"w")
    add_header.write(vcf_header)
    add_header.close()

    # write the variants to the new VCF
    input_vcf_df.to_csv(out_file, mode='a', header=False, index = False, sep = "\t", quoting=csv.QUOTE_NONE)


    #~~~~~~~~~~~~~~~~~~~
    # Sort the output
    #~~~~~~~~~~~~~~~~~~~
    cmd = "bcftools sort {} -o {}".format(out_file, out_file)
    logger.info("CMD: {}".format(cmd))
    subprocess.run(cmd, shell=True)

if __name__ == "__main__":

    main()



