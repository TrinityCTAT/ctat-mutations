#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

## compile_RNAediting_library_vcf.py --rediportal rediportal.txt --radar radar.txt  > RNAediting.library.vcf


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



class variant:

    def __init__(self, chr, pos, refallele, altallele, datatype):

        self.chr = chr
        self.pos = int(pos)
        self.refallele = refallele
        self.altallele = altallele
        self.datatype = set()

        self.add_datatype(datatype)
        
    def add_datatype(self, datatype):
        self.datatype.add(datatype)



def add_radar_annotations(editing_sites_dict, rediportal_filename):


    with open(rediportal_filename) as fh:
        header = next(fh)
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            chr_val = vals[0]
            pos_val = vals[1]
            strand = vals[3]
            (refallele, altallele) = ('A', 'G') if (strand == '+') else ('T', 'C')

            chrpos_token = ":".join([chr_val, pos_val, refallele, altallele])
            if chrpos_token in editing_sites_dict:
                variant_obj = editing_sites_dict[ chrpos_token ]
                variant_obj.add_datatype("radar")
            else:
                variant_obj = variant(chr_val, pos_val, refallele, altallele, 'radar')
                editing_sites_dict[ chrpos_token ] = variant_obj

    return editing_sites_dict


def add_rediportal_annotations(editing_sites_dict, rediportal_filename):

    with open(rediportal_filename) as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            chr_val = vals[0]
            chr_pos = vals[1]
            refallele = vals[2]
            altallele = vals[3]

            chrpos_token = ":".join([chr_val, chr_pos, refallele, altallele])

            if chrpos_token in editing_sites_dict:
                variant_obj = editing_sites_dict[ chrpos_token ]
                variant_obj.add_datatype("rediportal")
            else:
                variant_obj = variant(chr_val, chr_pos, refallele, altallele, 'rediportal')
                editing_sites_dict[ chrpos_token ] = variant_obj

    return editing_sites_dict


def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--radar", required=False, default=None, help="radar data file")
    parser.add_argument("--rediportal", required=False, default=None, help="rediportal data file")
    
    args=parser.parse_args()

    
    if not (args.rediportal or args.radar):
        print("Error, must set --rediportal and/or --radar", file=sys.stderr)
        sys.exit(1)

    editing_sites_dict = dict()

    if args.radar:
        editing_sites_dict = add_radar_annotations(editing_sites_dict, args.radar)

    if args.rediportal:
        editing_sites_dict = add_rediportal_annotations(editing_sites_dict, args.rediportal)


    ## output vcf file:
    print("##fileformat=VCFv4.2")
    print("##INFO=<ID=RNAEDIT,Type=String,Description=\"A known or predicted RNA-editing site\">")

    variant_objs = editing_sites_dict.values()

    variant_objs = sorted(variant_objs, key=lambda variant: (variant.chr, variant.pos))
    
    for variant_obj in variant_objs:
        
        print("\t".join([variant_obj.chr,
                         str(variant_obj.pos),
                         ".",
                         variant_obj.refallele,
                         variant_obj.altallele,
                         ".",
                         ".",
                         "RNAEDIT=" + ",".join(sorted(variant_obj.datatype)) ]))
        
    

    sys.exit(0)
    

if __name__ == '__main__':
    main()


