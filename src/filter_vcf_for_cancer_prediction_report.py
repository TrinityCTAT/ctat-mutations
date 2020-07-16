#!/usr/bin/env python

# Import libraries
import argparse
import csv
import gzip
import os
import codecs
import re

# Constants
# VCF related
CHR_COMMENT = "#"
STR_VCF_DELIMITER = "\t"
I_INFO_INDEX = 7
CHR_INFO_DELIMITER = ";"

# INFO features


STR_FATHMM = "FATHMM"
STR_CANCER = ["CANCER", "PATHOGENIC"]
STR_TISSUE = "TISSUE"
STR_CRAVAT = "CRV"

# Thresholds
MAX_PVALUE = 0.05
MAX_GNOMAD_AF = 0.01

# Arguments
prog_desc = "".join(["Filters VCF based on mutation priority predictions."])
prsr_arguments = argparse.ArgumentParser(
    prog="filter_vcf_for_predictions.py",
    description=prog_desc,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
prsr_arguments.add_argument("str_input_file", help="Input vcf.gz file.")
prsr_arguments.add_argument("str_output_file", help="Output filtered vcf file.")
# prsr_arguments.add_argument("tissue_type", help="Tissue type.")
args = prsr_arguments.parse_args()

# Stores the vcf info
lstr_vcf = []


def func_split_token(str_token):
    """
    Splits a VCF info token while guarenteeing 2 tokens
    * str_token : String token to split at '='
                : String
    * return : List of 2 strings
    """
    if str_token:
        lstr_pieces = str_token.split("=")
        i_pieces = len(lstr_pieces)
        if i_pieces == 1:
            return lstr_pieces + [""]
        if i_pieces == 2:
            return lstr_pieces
        elif i_pieces > 2:
            return [lstr_pieces[0], "=".join(lstr_pieces[1:])]
    return ["", ""]


# Get handle
str_input_ext = os.path.splitext(args.str_input_file)[1]
hndl_vcf = (
    gzip.open(args.str_input_file, "rb")
    if str_input_ext == ".gz"
    else open(args.str_input_file, "rb")
)


crv_field_positions_dict = dict()


def parse_crv_field_positions(crv_info_line):
    m = re.search('Format: (.*)">', crv_info_line)
    if not m:
        raise RuntimeError(
            "Error, not able to extract the Format info from crv info line: {}".format(
                crv_info_line
            )
        )

    crv_info_list_str = m.group(1)
    fields = crv_info_list_str.split("|")
    for idx in range(len(fields)):
        field = fields[idx]
        field_split_list = field.split("=")
        field_store = field_split_list[0]
        crv_field_positions_dict[idx] = field_store


# Read in vcf file
if args.str_input_file:
    with open(args.str_output_file, "w") as hndl_out:
        for lstr_line in csv.reader(
            codecs.iterdecode(hndl_vcf, "utf-8"), delimiter=STR_VCF_DELIMITER
        ):  # csv.reader(hndl_vcf, delimiter = STR_VCF_DELIMITER):

            if re.match("^##INFO=<ID=CRV,", lstr_line[0]):
                parse_crv_field_positions(lstr_line[0])

            # Do not keep any line that does not meet a criteria.
            f_keep = False

            # Keep comments
            if lstr_line[0][0] == CHR_COMMENT:
                lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))
                continue

            # Parse INFO tokens
            dict_info_tokens = dict(
                [
                    func_split_token(str_token)
                    for str_token in lstr_line[I_INFO_INDEX].split(CHR_INFO_DELIMITER)
                ]
            )
            # print dict_info_tokens

            ## ---------------------
            ## Prune common variants

            STR_COMMON_VARIANT = "RS"
            if STR_COMMON_VARIANT in dict_info_tokens:
                # automatic skip if variant is common
                continue

            STR_GNOMAD_AF = "gnomad_AF"
            if STR_GNOMAD_AF in dict_info_tokens:
                for af in dict_info_tokens[STR_GNOMAD_AF].split(","):
                    af = float(af)
                    # might want to change this logic in the future
                    if af >= MAX_GNOMAD_AF:
                        # consider it a common variant
                        continue

            ## ------------------------
            ## Chasm and Vest selection

            if STR_CRAVAT in dict_info_tokens:
                cravat_str = dict_info_tokens[STR_CRAVAT]
                cravat_fields = cravat_str.split("|")
                for cravat_idx in range(len(cravat_fields)):
                    print(
                        "{}:\t{}".format(
                            crv_field_positions_dict[cravat_idx],
                            cravat_fields[cravat_idx],
                        )
                    )

            print("//\n\n\n")
            continue

            STR_CHASM_PVALUE = "chasmplus_pval"
            # Otherwise require CRAVAT or VEST to have an annotation.
            if (
                STR_CHASM_PVALUE in dict_info_tokens
                and float(dict_info_tokens.get(STR_CHASM_PVALUE, "1")) <= MAX_PVALUE
            ):
                f_keep = True

            STR_VEST_PVALUE = "vest_pval"
            if (
                STR_VEST_PVALUE in dict_info_tokens
                and float(dict_info_tokens.get(STR_VEST_PVALUE, "1")) <= MAX_PVALUE
            ):
                f_keep = True

            ## ----------------------
            ## Retain COSMIC entries

            # Keep everything that is a COSMIC ID at this point.
            STR_COSMIC_ID = "COSMIC_ID"
            # Keep FATHMM = Cancer and FATHMM=PATHOGENIC
            if (
                STR_COSMIC_ID in dict_info_tokens
                or dict_info_tokens.get(STR_FATHMM, "") in STR_CANCER
            ):
                f_keep = True

            # Store passing variant
            if f_keep:
                lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))

        # Last write of buffer
        hndl_out.write("\n".join(lstr_vcf))

# Close input handle
hndl_vcf.close()
