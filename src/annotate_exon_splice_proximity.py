#!/usr/bin/env python


import argparse
import sys, os
import subprocess
from collections import defaultdict
import logging
import re

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../PyLib"]))
import ctat_util

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


## uses ref exon bed file:
## m ../ref_annot.gtf | perl -lane 'if ($F[2] eq "exon") { print;}' | print.pl 0 3 4 | sort -k1,1 -k2,2 -k3,3 -u > ref_exons.bed

MAX_DIST_TO_SPLICE = 10  # note this is based on the bed tools intersect target in the mutation lib

def main():
    
    #add options to inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Adds exon splice distance annotations to vcf file (report up to len 10 distance away from splice).\n")
    
    parser.add_argument('--input_vcf', required=True, help="input vcf file")
    parser.add_argument('--ctat_mutation_lib_dir', required=True, help='path to ctat mutation lib dir')
    parser.add_argument('--output_vcf', required=True, help="output vcf file including annotation for distance to splice neighbor")

    parser.add_argument('--debug', default=False, action='store_true', help='debug mode, retains temporary intermediate files')
    parser.add_argument("--tmpdir", default="/tmp", help="tmp directory")

    args = parser.parse_args()
    


    input_vcf_file = args.input_vcf
    ctat_mutation_lib_dir = args.ctat_mutation_lib_dir
    output_vcf_file = args.output_vcf
    DEBUG_MODE = args.debug

    tmpdir = args.tmpdir


    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)


    ref_splice_adj_regions_bed = os.path.join(ctat_mutation_lib_dir, "ref_annot.splice_adj.bed.gz")
    if not os.path.exists(ref_splice_adj_regions_bed):
        logger.critical("cannot locate required file: {}".format(ref_splice_adj_regions_bed))
        sys.exit(1)
            
    

    splice_adjacent_snps_vcf_file = os.path.join(tmpdir, os.path.basename(input_vcf_file) + ".splice_adjacent.vcf")
    cmd = "bedtools intersect -wb -a {} -b {} > {}".format(input_vcf_file, ref_splice_adj_regions_bed, splice_adjacent_snps_vcf_file)
    logger.info("CMD: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)

    # capture those variants that were identified as exon-adjacent
    exon_splice_adj_dict = dict()
    with open(splice_adjacent_snps_vcf_file) as fh:
        for line in fh:
            if re.match("#", line): continue
            line = line.rstrip()
            vals = line.split("\t")
            chr_val = vals[0]
            pos_val = vals[1]
            chrpos = "{}:{}".format(chr_val, pos_val)

            # determine distance from splice boundary
            splice_adj_name = vals[-1]
            # ex. chr1:L:6586064:6586073 or chr1:R:3662841:3662850 
            #
            # where L indicates left intron side and R indicates right intron side:
            #   ---->GT   intron     AG<-------
            #       |<-----|    |----->|
            #         L              R     segment type
            #

            (chromosome, segment_type, lend, rend) = splice_adj_name.split(":")
            
            pos_val = int(pos_val); lend = int(lend); rend = int(rend)
            # distance from splice boundary
            delta = pos_val - lend + 1 if segment_type == "L" else rend - pos_val + 1

            if not (delta > 0 and delta <= MAX_DIST_TO_SPLICE):
                continue
            
            # want to store the minimum distance in case there are competing introns.
            if chrpos in exon_splice_adj_dict:
                if delta < exon_splice_adj_dict[ chrpos ]:
                    exon_splice_adj_dict[ chrpos ] = delta
            else:
                # first occurence
                exon_splice_adj_dict[ chrpos ] = delta

                
    logger.info("Adding exon-neighboring annotations (SPLICEADJ = splice ajacent).")
    ## make output a vcf formatted file:
    with open(output_vcf_file, 'w') as ofh:
        
        with ctat_util.open_file_for_reading(input_vcf_file) as fh:
            for line in fh:
                if line[0] == "#":
                    
                    if re.match("#CHROM\t", line):
                        # add header info line for the repeat annotation type
                        ofh.write("##INFO=<ID=SPLICEADJ,Number=1,Type=Integer,Description=\"Variant is within distance of {} to a reference exon splice boundary\">\n".format(MAX_DIST_TO_SPLICE))
                        
                    ofh.write(line)
                else:
                    line = line.rstrip()
                    vals = line.split("\t")
                    chrpos = "{}:{}".format(vals[0], vals[1])
                    if chrpos in exon_splice_adj_dict:
                        vals[7] += ";SPLICEADJ={}".format(exon_splice_adj_dict[chrpos])

                    ofh.write("\t".join(vals) + "\n")


    # cleanup
    if not DEBUG_MODE:
        pass
    
    sys.exit(0)
    

if __name__ == "__main__":

    main()
