#!/usr/bin/env python

import sys, os, re
import gzip
from collections import defaultdict

field_counter = defaultdict(int)

fields_auditing = { "AC",
                    "ALT",
                    "BaseQRankSum",
                    "DJ",
                    "DP",
                    "ED",
                    "Entropy",
                    "ExcessHet",
                    "FS",
                    "Homopolymer",
                    "LEN",
                    "MLEAF",
                    "MMF",
                    "QUAL",
                    "REF",
                    "RPT",
                    "RS",
                    "ReadPosRankSum",
                    "SAO",
                    "SOR",
                    "SPLICEADJ",
                    "TCR",
                    "TDM",
                    "VAF",
                    "VMMF",
                    "GT_1/2" }
                    

def main():

    if len(sys.argv) < 2:
        exit("\n\n\tusage: {} input.vcf\n\n".format(sys.argv[0]))

    input_vcf = sys.argv[1]

    if re.search(".gz$", input_vcf):
        fh = gzip.open(input_vcf, "rt", encoding='utf-8')
    else:
        fh = open(input_vcf, "rt", encoding='utf-8')

    for line in fh:
        if line[0] == "#":
            continue

        line = line.rstrip()
        vals = line.split("\t")

        info_line = vals[7]
        count_fields(info_line)

    fh.close()

    dump_counts()

    sys.exit(0)



def count_fields(info_line):
    fields = info_line.split(";")
    for field in fields:
        key_vals = field.split("=")
        key = key_vals[0]
        if key in fields_auditing:
            field_counter[key] += 1
                     



def dump_counts():

    keys = list(field_counter.keys())
    sorted(keys)

    for key in sorted(list(fields_auditing)):
        print("\t".join([key, str(field_counter.get(key, 0))]))


if __name__=='__main__':
    main()
