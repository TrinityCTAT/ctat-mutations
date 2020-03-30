#!/usr/bin/env python

# Import libraries
import sys
import argparse
import csv
import gzip
import os
from collections import defaultdict

# Constants
CHR_COMMENT = "#"
STR_VCF_DELIMITER = "\t"
I_INFO_INDEX = 7
CHR_INFO_DELIMITER = ";"
STR_COSMIC_ID = "COSMIC_ID"
STR_COMMON_VARIANT = "RS"
STR_DEPTH = "DP"
#STR_ORIGIN = "SAO"
#STR_ORIGIN_GERMLINE = "1"
STR_FATHMM = "FATHMM"
#STR_FATHMM_NEUTRAL = "NEUTRAL"
STR_RNAEDIT = "RNAEDIT"


# Arguments
prog_desc = "".join(["Extracts common variants which do",
                     "not have COSMIC ids entries from a vcf file"])
prsr_arguments = argparse.ArgumentParser(prog="filter_vcf_for_cancer.py",
    description=prog_desc,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
prsr_arguments.add_argument("str_input_file", help="Input vcf.gz file.")
prsr_arguments.add_argument("str_output_file", help="Output filtered vcf file.")
args = prsr_arguments.parse_args()

# Stores the vcf info
lstr_vcf = []
# Buffer size in lines
i_write_amount = 1000


def func_split_token(str_token):
    """
    Splits a VCF info token while guarenteeing 2 derived pieces.
    * str_token : String token to split at '='
        : String

        * return : List of 2 strings
    """

    if str_token:
        lstr_pieces = str_token.split("=")
        i_pieces = len(lstr_pieces)
    if i_pieces == 1:
        return(lstr_pieces + [""])
    if i_pieces == 2:
        return lstr_pieces
    elif i_pieces > 2:
        return [lstr_pieces[0], "=".join(lstr_pieces[1:])]
    return ["", ""]

# Get handle
str_input_ext = os.path.splitext(args.str_input_file)[1]
hndl_vcf = gzip.open(args.str_input_file, "rt") if str_input_ext == ".gz" else open(args.str_input_file, "r")


filter_counter = defaultdict(int)
pass_counter = 0

# Read in vcf file

with open(args.str_output_file, "w") as hndl_out:
    for lstr_line in csv.reader(hndl_vcf, delimiter = STR_VCF_DELIMITER):

        # Keep comments
        if lstr_line[0][0] == CHR_COMMENT:
            lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))
            continue

        # Parse INFO tokens
        dict_info_tokens = dict([func_split_token(str_token) for str_token in lstr_line[I_INFO_INDEX].split(CHR_INFO_DELIMITER)])


        ## automatic pass if in cosmic:
        if STR_COSMIC_ID in dict_info_tokens:
            # Store passing variant
            lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))
            pass_counter += 1
            continue

        ##
        # FATHMM=PATHOGENIC
        if STR_FATHMM in dict_info_tokens and dict_info_tokens[STR_FATHMM] in ("PATHOGENIC","CANCER"):
            # Store passing variant
            lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))
            pass_counter += 1
            continue


        #########################
        ## Negative filters below

        filter_flag = False

        if STR_RNAEDIT in dict_info_tokens:
            filter_counter['rnaedit'] += 1
            filter_flag = True

        ## Filter out common variant that do not have cosmic ids
        if STR_COMMON_VARIANT in dict_info_tokens:
            if (dict_info_tokens[STR_COMMON_VARIANT]):
                filter_counter['common'] += 1
                filter_flag = True

        ## Filter DP < 10
        if STR_DEPTH in dict_info_tokens:
            if (int(dict_info_tokens[STR_DEPTH]) < 10):
                filter_counter['insufficient_depth'] += 1
                filter_flag = True


        ## No gene match
        if not dict_info_tokens["GENE"]:
            filter_counter['no_gene'] += 1
            filter_flag = True


        ############################
        ############################

        # Store passing variant
        if not filter_flag:
            lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))
            pass_counter += 1

        # If buffer is large enough, write to file
        if len(lstr_vcf) >= i_write_amount:
            str_write = "\n".join(lstr_vcf)
            if str_write:
                str_write = str_write + "\n"
                hndl_out.write(str_write)
                lstr_vcf = []

    # Last write of buffer
    hndl_out.write("\n".join(lstr_vcf))


# Close input handle
hndl_vcf.close()

print("Audit for pre-cravat filter:" + str(filter_counter), file=sys.stderr)
print("Count of variants prepped for cravat: {}".format(pass_counter), file=sys.stderr)

sys.exit(0)

