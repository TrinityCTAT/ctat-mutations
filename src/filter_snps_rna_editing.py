#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import csv

# Constants
CHR_COMMENT="#"
STR_VCF_DELIMITER="\t"
I_CHR_INDEX=0
I_POS_INDEX=1
I_REF_INDEX=3
I_ALT_INDEX=4


#RADAR columns
STR_RADAR_DELIMITER="\t"
I_CHR_RADAR=0
I_POS_RADAR=1
I_STRAND_RADAR=3


#REDIPortal columns
STR_REDIPORTAL_DELIMITER="\t"
I_CHR_REDIPORTAL=0
I_POS_REDIPORTAL=1
I_STRAND_REDIPORTAL=4
I_INCHR_REDIPORTAL=2
I_INRNA_REDIPORTAL=3


str_prog_des = " ".join(["Filters SNP entries from a VCF file which are",
                         "in provided RNA Editing databases."])
fmtr_cur = argparse.ArgumentDefaultsHelpFormatter
prsr_arguments = argparse.ArgumentParser(prog="filter_snps_rna_editing.py",
                                         description=str_prog_des,
                                         formatter_class=fmtr_cur)
prsr_arguments.add_argument("--radar",
                            action="store",
                            dest="str_radar_db",
                            help=" ".join(["Radar database for RNA editing",
                                           "matched to the correct species",
                                           "reference."]))

prsr_arguments.add_argument("--rediportal",
                            action="store",
                            dest="str_rediportal_db",
                            help=" ".join(["REDIPortal database for RNA editing",
                                           "matched to the correct species",
                                           "reference."]))
prsr_arguments.add_argument("str_input_file",
                            help="Input vcf file.")
prsr_arguments.add_argument("str_output_file",
                            help="Output SNP vcf file.")
args = prsr_arguments.parse_args()

# RNA editing
dict_radar = {}
dict_rediportal = {}
dict_seq_comp = {"A":"T","T":"A","G":"C","C":"G"}
dict_seq_interp = {"A":"A","C":"C","G":"G","T":"T","I":"G","U":"T"}

# Read in the RADAR data if given
if args.str_radar_db:
  with open(args.str_radar_db, "r") as hndl_radar:
    for lstr_line in csv.reader(hndl_radar, delimiter=STR_RADAR_DELIMITER):

      if not lstr_line:
        continue

      if lstr_line[I_STRAND_RADAR] == "+":
        str_rpos = lstr_line[I_CHR_RADAR].lower()+"-"+lstr_line[I_POS_RADAR]
        dict_radar[str_rpos] = ["A", dict_seq_interp["I"]]

      elif lstr_line[I_STRAND_RADAR] == "-":
        str_rneg = lstr_line[I_CHR_RADAR].lower()+"-"+lstr_line[I_POS_RADAR]
        dict_radar[str_rneg] = [dict_seq_comp["A"],
                                dict_seq_comp[dict_seq_interp["I"]]]



# Read in the REDIPORTAL data if given
if args.str_rediportal_db:
  with open(args.str_rediportal_db, "r") as hndl_rediportal:
    for lstr_line in csv.reader(hndl_rediportal, delimiter=STR_REDIPORTAL_DELIMITER):

      if not lstr_line:
        continue

      inchr = lstr_line[I_INCHR_REDIPORTAL]
      inrna = lstr_line[I_INRNA_REDIPORTAL]


      if lstr_line[I_STRAND_REDIPORTAL] == "+" or lstr_line[I_STRAND_REDIPORTAL] == "-":
        str_rpos = lstr_line[I_CHR_REDIPORTAL].lower()+"-"+lstr_line[I_POS_REDIPORTAL]
        dict_rediportal[str_rpos] = [dict_seq_interp[inchr],
                                 dict_seq_interp[inrna]]

      # elif lstr_line[I_STRAND_REDIPORTAL] == "-":
      #   str_rneg = lstr_line[I_CHR_REDIPORTAL].lower()+"-"+lstr_line[I_POS_REDIPORTAL]
      #   dict_rediportal[str_rneg] = [dict_seq_comp[dict_seq_interp[inchr]],
      #                            dict_seq_comp[dict_seq_interp[inrna]]]


# Stores the vcf info
lstr_vcf = []

# Read in vcf file
if args.str_input_file:
  with open(args.str_output_file, "w") as hndl_out:
    with open(args.str_input_file, "r") as hndl_vcf:
      for lstr_line in csv.reader(hndl_vcf, delimiter=STR_VCF_DELIMITER):

        # Keep comments
        if lstr_line[0][0] == CHR_COMMENT:
          lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))
          continue

        # Get ALT / REF
        str_alt = lstr_line[I_ALT_INDEX]
        str_ref = lstr_line[I_REF_INDEX]

        # Create ID (chr-pos)
        str_id = lstr_line[I_CHR_INDEX].lower()+"-"+lstr_line[I_POS_INDEX]

        #filter REDIportal
        if dict_rediportal:
          if str_id in dict_rediportal:
            lstr_rna_edit_entry = dict_rediportal[str_id]
            if(lstr_rna_edit_entry[0] == str_ref
               and lstr_rna_edit_entry[1] == str_alt):
              continue

        # Filter Radar
        if dict_radar:
          if str_id in dict_radar:
            lstr_rna_edit_entry = dict_radar[str_id]
            if(lstr_rna_edit_entry[0] == str_ref
               and lstr_rna_edit_entry[1] == str_alt):
              continue

        # Store SNP
        lstr_vcf.append(STR_VCF_DELIMITER.join(lstr_line))

      for str_out_line in lstr_vcf:
        hndl_out.write(str_out_line + "\n")
