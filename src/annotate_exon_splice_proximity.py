#!/usr/bin/env python


import argparse
import sys, os
import subprocess
from collections import defaultdict
import logging
import re

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


## uses ref exon bed file:
## m ../ref_annot.gtf | perl -lane 'if ($F[2] eq "exon") { print;}' | print.pl 0 3 4 | sort -k1,1 -k2,2 -k3,3 -u > ref_exons.bed

def main():
    
    #add options to inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Adds repeat feature annotations to vcf file.\n")
    
    parser.add_argument('--input_vcf', required=True, help="input vcf file")
    parser.add_argument('--ref_exon_bed', required=True, help='reference exon bed file')
    parser.add_argument('--output_vcf', required=True, help="output vcf file including annotation for distance to splice neighbor")
    parser.add_argument("--max_dist_to_splice",  type=int, default=6, help="maximum distance from within intron to a splice boundary")

    parser.add_argument('--debug', default=False, action='store_true', help='debug mode, retains temporary intermediate files')
    parser.add_argument("--tmpdir", default="/tmp", help="tmp directory")

    args = parser.parse_args()
    


    input_vcf_file = args.input_vcf
    ref_exon_bed_file = args.ref_exon_bed
    max_dist_to_splice = args.max_dist_to_splice
    output_vcf_file = args.output_vcf
    DEBUG_MODE = args.debug

    tmpdir = args.tmpdir


    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)


    ## capture the variants that overlap exons
    exon_overlapping_variants_vcf = os.path.join(tmpdir, os.path.basename(input_vcf_file) + ".exon_overlapping.vcf")

    cmd = "bedtools intersect -header -a {} -b {} > {}".format(input_vcf_file, ref_exon_bed_file, exon_overlapping_variants_vcf)
    logger.info("CMD: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)

    ## get the remaining variants that do not overlap exons
    exon_nonoverlapping_vcf = os.path.join(tmpdir, os.path.basename(input_vcf_file) + ".exon_NONoverlapping.vcf")
    cmd = "bedtools subtract -header -a {} -b {} > {}".format(input_vcf_file, exon_overlapping_variants_vcf, exon_nonoverlapping_vcf)
    logger.info("CMD: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)
    
    ## see if the remaining exon-NONoverlapping variants are within distance
    # make a new snp vcf file that includes max_dist_to_splice padding on each side and see if that overlaps

    snp_padded_bed_file = os.path.join(tmpdir, os.path.basename(input_vcf_file) + ".padded_{}_variants.bed".format(max_dist_to_splice))
    with open(snp_padded_bed_file, "w") as ofh:
        with open(exon_nonoverlapping_vcf) as fh:
            for line in fh:
                if re.match("#", line):
                    continue
                vals = line.split("\t")
                chr = vals[0]
                pos = int(vals[1])
                lend = pos - max_dist_to_splice + 1
                rend = pos + max_dist_to_splice - 1
                chrpos = "{}:{}".format(chr, pos)
                
                ofh.write("\t".join([chr, str(lend), str(rend), chrpos]) + "\n")


    splice_adjacent_snps_bed = os.path.join(tmpdir, os.path.basename(input_vcf_file) + ".splice_adjacent.bed")
    cmd = "bedtools intersect -a {} -b {} > {}".format(snp_padded_bed_file, ref_exon_bed_file, splice_adjacent_snps_bed)
    logger.info("CMD: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)

    # capture those variants that were identified as exon-adjacent
    var_set = set()
    with open(splice_adjacent_snps_bed) as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            chrpos = vals[3]
            var_set.add(chrpos)
    
    
    logger.info("Adding exon-neighboring annotations (SPLICEADJ = splice ajacent).")
    ## make output a vcf formatted file:
    with open(output_vcf_file, 'w') as ofh:
        
        with open(input_vcf_file, 'r') as fh:
            for line in fh:
                if line[0] == "#":
                    
                    if re.match("#CHROM\t", line):
                        # add header info line for the repeat annotation type
                        ofh.write("##INFO=<ID=SPLICEADJ,Number=1,Type=Integer,Description=\"Variant is within {} distance to a reference exon splice boundary\">\n".format(max_dist_to_splice))

                    ofh.write(line)
                else:
                    line = line.rstrip()
                    vals = line.split("\t")
                    chrpos = "{}:{}".format(vals[0], vals[1])
                    if chrpos in var_set:
                        vals[7] += ";SPLICEADJ=1"

                    ofh.write("\t".join(vals) + "\n")


    # cleanup
    if not DEBUG_MODE:
        pass
    
    sys.exit(0)
    

if __name__ == "__main__":

    main()
