#!/usr/bin/env python


import argparse
import sys, os
import subprocess
from collections import defaultdict
import logging
import re
import math

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../PyLib"]))
import ctat_util


logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def contains_homopolymer(window_seq, edit_base_nuc):
    
    homopolymer = False

    window_seq = window_seq.upper()
    edit_base_nuc = edit_base_nuc.upper()
    
    # iterate of the sequence 
    #   check if the 4 consecutive nucleotides are the edited nucleotide 
    #   if true will set th variable homopolymer to TRUE
    for k in range(len(window_seq)-4):
        if (edit_base_nuc == window_seq[k] and
            edit_base_nuc == window_seq[k+1] and
            edit_base_nuc == window_seq[k+2] and
            edit_base_nuc == window_seq[k+3] ) :

            return True
    
    # no homopolymer found
    return False



def compute_entropy(window_seq):

    window_seq = window_seq.upper()
    
    window_length = len(window_seq)

    base_counter_dict = defaultdict(int)

    for nucleotide in window_seq:
        base_counter_dict[ nucleotide ] += 1


    entropy = 0
    for nucleotide, count in base_counter_dict.items():

        p_val = count / window_length

        entropy += -1 * p_val * math.log2(p_val)

    entropy = "{:0.3f}".format(entropy)
    
    return entropy


def main():
    
    #add options to inputs
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Adds repeat feature annotations to vcf file.\n")
    
    parser.add_argument('--input_vcf', required=True, help="input vcf file")
    parser.add_argument('--ref_genome_fa', required=True, help='reference  file')
    parser.add_argument('--output_vcf', required=True, help="output vcf file including annotation for distance to splice neighbor")
    parser.add_argument("--window_size",  type=int, default=10, help="window length centered at variant position for sequence analysis")
    
    parser.add_argument('--debug', default=False, action='store_true', help='debug mode, retains temporary intermediate files')
    parser.add_argument("--tmpdir", default="/tmp", help="tmp directory")

    args = parser.parse_args()
    


    input_vcf_file = args.input_vcf
    genome_fa_file = args.ref_genome_fa
    window_size = args.window_size
    output_vcf_file = args.output_vcf
    DEBUG_MODE = args.debug

    if DEBUG_MODE:
        logger.setLevel('DEBUG')

    tmpdir = args.tmpdir

    half_window = int(window_size/2.0 + 0.5)

    # 1) Create the bed file containing the windowed variants
    #------------------------
    windowed_variants_bed_file = os.path.join(tmpdir, os.path.basename(input_vcf_file) + ".windowed_{}_variants.bed".format(window_size))
    with open(windowed_variants_bed_file, 'w') as ofh:
        with ctat_util.open_file_for_reading(input_vcf_file) as fh:
            for line in fh:
                if re.match("#", line):
                    continue
                vals = line.split("\t")
                chr_val = vals[0]
                chr_coord = int(vals[1])
                edit_base_nuc = vals[4]
                
                chrpostoken = "{}:{}:{}".format(chr_val, chr_coord, edit_base_nuc)
                
                ofh.write("\t".join([chr_val, str(chr_coord-half_window), str(chr_coord+half_window), chrpostoken]) + "\n")



    ## get fasta coordinates for feature
    window_seqs_file = windowed_variants_bed_file + ".seqs"
    cmd = "fastaFromBed -name -tab -fi {} -bed {} -fo {}".format(genome_fa_file, windowed_variants_bed_file, window_seqs_file)
    logger.info("CMD: {}".format(cmd))
    subprocess.check_call(cmd, shell=True)
    
    chrpos_homopolymer_set = set()
    chrpos_entropy_dict = dict()

    with open(window_seqs_file) as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            chrpostoken = vals[0]
            window_seq = vals[1]
            print(line)
            print(chrpostoken)
            print(chrpostoken.split(":"))
            (chrom, position, edit_base_nuc) = chrpostoken.split(":")[0:3]
            
            chrpos = "{}:{}".format(chrom, position)
            
            # set the constant homopolymer to false for this line 
            homopolymer_flag = contains_homopolymer(window_seq, edit_base_nuc)
            if homopolymer_flag:
                chrpos_homopolymer_set.add(chrpos)

            entropy_val = compute_entropy(window_seq)
            chrpos_entropy_dict[chrpos] = entropy_val

            if DEBUG_MODE:
                logger.debug("\t".join([chrpostoken, window_seq, "entropy:{}".format(entropy_val), "homopolymer:{}".format(homopolymer_flag)]))

    ###############################################
    ## Add feature annotations to original vcf file
    
    logger.info("Adding entropy and homopolymer annotations")
    ## make output a vcf formatted file:
    with open(output_vcf_file, 'w') as ofh:
        
        with ctat_util.open_file_for_reading(input_vcf_file) as fh:
            for line in fh:
                if line[0] == "#":
                    
                    if re.match("#CHROM\t", line):
                        # add header info line for the repeat annotation type
                        ofh.write("##INFO=<ID=Homopolymer,Number=1,Type=Integer,Description=\"Variant is located in or near a homopolymer sequence\">\n")
                        ofh.write("##INFO=<ID=Entropy,Number=1,Type=Float,Description=\"Entropy for sequence in window of length {} centered at the variant position\">\n".format(window_size))
                        
                    ofh.write(line)
                else:
                    line = line.rstrip()
                    vals = line.split("\t")
                    chrpos = "{}:{}".format(vals[0], vals[1])
                    if chrpos in chrpos_homopolymer_set:
                        vals[7] += ";Homopolymer=1"

                    entropy_val = chrpos_entropy_dict[chrpos]
                    vals[7] += ";Entropy={}".format(entropy_val)
                    
                    ofh.write("\t".join(vals) + "\n")


    # cleanup
    if not DEBUG_MODE:
        pass


    
    sys.exit(0)
    

if __name__ == "__main__":

    main()
