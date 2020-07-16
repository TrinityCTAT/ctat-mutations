#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

# Constants
## Comments
CHR_COMMENT = "#"
STR_INFO = "##INFO="
STR_VCF_DELIMITER = str("\t")
STR_FILTER = "##FILTER="
STR_DESCRIPTION = "Description"
STR_ID = "ID"
STR_NUMBER = "Number"
STR_TYPE = "Type"
STR_TYPE_STRING = "String"
LSTR_REQUIRED_INFO_FIELDS = [STR_ID, STR_NUMBER, STR_TYPE, STR_DESCRIPTION]
LSTR_REQUIRED_FILTER_FIELDS = [STR_ID, STR_DESCRIPTION]
## VCF entries related
I_INFO_COL = 7

import sys, os, re
import argparse
import csv

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../PyLib"]))
import ctat_util



prsr_arguments = argparse.ArgumentParser(prog="groom_vcf.py",
                                         description="Cleans VCF files and fixes formatting errors.",
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
prsr_arguments.add_argument("str_input_file",
                            help="Input vcf file.")
prsr_arguments.add_argument("str_output_file",
                            help="Output groomed vcf file.")
args = prsr_arguments.parse_args()

# Stores the vcf info
lstr_vcf = []

# Read in vcf file
if args.str_input_file:
    with ctat_util.open_file_for_reading(args.str_input_file) as hndl_vcf:
        for lstr_line in csv.reader(hndl_vcf, delimiter = STR_VCF_DELIMITER):
            
            # Work with comments
            if lstr_line[0][0] == CHR_COMMENT:
                # Make sure there are not spaces (especially at the end).
                lstr_line = [str_token.strip(" ") for str_token
                             in lstr_line if str_token]
                # Check if the comment is a type of interest
                str_header_type = None
                if (lstr_line[0][: len(STR_INFO)] == STR_INFO):
                    str_header_type = STR_INFO
                if (lstr_line[0][: len(STR_FILTER)] == STR_FILTER):
                    str_header_type = STR_FILTER

                if str_header_type:
                    # GATK does not like certain 4.2 features so, remove them here.
                    # COSMIC does not include number in the id so we have to add one
                    # the most liberal one for String is ".", which we use here.
                    if not(len(lstr_line) == 1):
                        print("Error, expected one token in the comment received more.")
                        print(str(lstr_line))
                        exit(101)

                    # Pull off header type info
                    str_header_token = lstr_line[0][len(str_header_type):]
                    # Remove brackets
                    str_header_token = str_header_token.strip("<>")
                    # Remove , and = from idescription free text.
                    i_description_text_start = str_header_token.index("Description=\"")
                    i_description_text_start += len("Description=\"")
                    i_description_next_quote = str_header_token[i_description_text_start:].index("\"")
                    i_description_next_quote += i_description_text_start
                    str_freetext = str_header_token[i_description_text_start: i_description_next_quote]
                    str_freetext = str_freetext.replace("=",":")
                    str_freetext = str_freetext.replace(","," ")
                    str_header_token = str_header_token[:i_description_text_start] + str_freetext + str_header_token[i_description_next_quote:]
                    # Split tokens and store as dict
                    dict_header_token = dict([str_header_token_piece.split("=")
                                              for str_header_token_piece
                                              in str_header_token.split(",")])
                    # Work with info field formats
                    if str_header_type == STR_INFO:
                        if (not STR_NUMBER in dict_header_token) and (STR_TYPE in dict_header_token):
                            if(dict_header_token[STR_TYPE].lower() == STR_TYPE_STRING.lower()):
                                dict_header_token[STR_NUMBER] = "."
                        for str_required_info in LSTR_REQUIRED_INFO_FIELDS:
                            if not str_required_info in dict_header_token:
                                print("INFO fields require the following fields " + ",".join(LSTR_REQUIRED_INFO_FIELDS))
                                print("This INFO line is missing fields")
                                print(STR_VCF_DELIMITER.join(lstr_line))
                                exit(102)

                        if re.match("extra_vcf_info", dict_header_token[STR_ID]):
                            # cravat added a field we dont want
                            continue
                                
                        # Order and write tokens
                        str_header_token = ",".join(["=".join([str_key,dict_header_token[str_key]]) for str_key in LSTR_REQUIRED_INFO_FIELDS])
                        str_header_token = str_header_type +"<" + str_header_token + ">"
                        lstr_line = [str_header_token]

                    # Work with the filter fields
                    if str_header_type == STR_FILTER:
                        for str_required_filter in LSTR_REQUIRED_FILTER_FIELDS:
                            if not str_required_filter in dict_header_token:
                                print("FILTER fields require the following fields " + ",".join(LSTR_REQUIRED_FILTER_FIELDS))
                                print("This FILTER line is missing fields")
                                print(STR_VCF_DELIMITER.join(lstr_line))
                                exit(103)
                        # Order and write tokens
                        lstr_token_order = LSTR_REQUIRED_FILTER_FIELDS
                        if len(LSTR_REQUIRED_FILTER_FIELDS) < len(dict_header_token):
                            lstr_token_order = lstr_token_order + [str_key for str_key in dict_header_token.keys() if not str_key in LSTR_REQUIRED_FILTER_FIELDS]
                        str_header_token = ",".join(["=".join([str_key,dict_header_token[str_key]]) for str_key in lstr_token_order])
                        str_header_token = str_header_type +"<" + str_header_token + ">"
                        lstr_line = [str_header_token]

                # Store
                lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))
                continue

            # Not a commment
            else:

                ## only yielding 'chr*' contigs
                if not re.search("^chr", lstr_line[0]):
                    continue
                
                # Remove spaces in the INFO data
                lstr_line[I_INFO_COL] = lstr_line[I_INFO_COL].replace(" ","_")

                # strip out cravat-added xtra fields we don't want.
                new_I_INFO_pts = list()
                for I_info_pt in lstr_line[I_INFO_COL].split(";"):
                    if not re.match("extra_vcf_info", I_info_pt):
                        new_I_INFO_pts.append(I_info_pt)
                lstr_line[I_INFO_COL] = ";".join(new_I_INFO_pts)
                            

            # Store
            lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))

    # Write cleaned lines to files.
    with open(args.str_output_file, "w") as hndl_out:
        for str_out_line in lstr_vcf:
            # some last bit of polishing
            str_out_line = re.sub(" +", " ", str_out_line)
            str_out_line = re.sub("_+", "_", str_out_line)
            str_out_line = re.sub("--+", "--", str_out_line)
            hndl_out.write(str_out_line + "\n")
