#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import argparse
import datetime
import os,sys,re
import logging
import subprocess
import multiprocessing as mp
import gzip

logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)

'''
The following adds annotations to variants in a provided VCF.

Input:
    VCF: File containing the variants of interest  and the bam file
    BAM: THE BAM file that was used in the creation of the VCF input.
Output:
    VCF: The input VCF with annotations added to the INFO section

'''
def processVCFHead(line, outfile):
    '''
    Add the annotation flag description to the head of the VCF 
        VPR : Variant Passed Reads, reads that PASS filtering that contain the variation 
        TPR : Total Passed Reads , reads that PASS filtering 
        TDM : Total Duplicate Marked, number of reads that are duplicate marked
        VAF: Variant allele fraction
        TMMR: Number of multi-mapped reads
        MMF: multi-mapped read fraction of TPR
    '''
    
    if line[0] == "#":
        if re.match("#CHROM\t", line):
            # add header info line for the repeat annotation type
            # once you get to the #CHROM line, 
            #       write the new INFO annotations above the #CHROM header line 
            outfile.write("##INFO=<ID=VPR,Number=1,Type=Integer,Description=\"Variant Passed Reads, reads that PASS filtering that contain the variation\">\n")
            outfile.write("##INFO=<ID=TPR,Number=1,Type=Integer,Description=\"Total Passed Reads , reads that PASS filtering\">\n")
            outfile.write("##INFO=<ID=TDM,Number=1,Type=Integer,Description=\"Total Duplicate Marked, number of reads that are duplicate marked \">\n")
            outfile.write("##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele fraction (VPR / TPR) \">\n")
            outfile.write("##INFO=<ID=TMMR,Number=1,Type=Integer,Description=\"Total multi-mapped reads at site\">\n")
            outfile.write("##INFO=<ID=MMF,Number=1,Type=Float,Description=\"Multi-mapped read fraction (TMMR / TPR) \">\n")

        outfile.write(line)

def check_reverse_strand(sam_flag):
    #----------------------------------------------------
    # Check SAM flag
    #----------------------------------------------------
    # set the reverse indicator to 0
    revers_strand = 0 
    # revers_strand = 1 if (alignment & 16)
    # check the SAM flag, (1000 means this is a revers strand)
    binary_flag = bin(int(sam_flag))
    if len(binary_flag) >= 5:
        rev_test = binary_flag[-5]
        if int(rev_test) == 1:
            revers_strand = 1

    return revers_strand


def check_duplicate_marked(sam_flag):

    binary_flag = bin(int(sam_flag))
    
    # set the multi-mapping indicator to 0
    duplicateMarked = 0 
    # duplicateMarked = 1 if (alignment & 1024)
    # check the SAM flag, (10000000000 means this is multi-mapping)
    if len(binary_flag) >= 11:
        duplicate_test = binary_flag[-11] 
        if int(duplicate_test) == 1:
            duplicateMarked = 1

    return duplicateMarked
    


def evaluate_PASS_reads(i, bamFile):

    vals = i.split("\t")

    #sys.stderr.write("[{}]\n".format(vals[0] + "::" + vals[1]))

    #--------------
    # Constants 
    #--------------
    quality_offset = 33
    minimal_base_quality = 25
    minimum_mismatch = 1
    not_filtered, removed= 0, 0

    newmismatch = 0
    mismatchreadcount = 0
    newcov, newmismatch = 0, 0 
    basequalFail, readPosFail = 0, 0
    duplicateMarked = 0
    multimappedReadCount = 0
    output =""
    output_failed =""
    
    # process the input line 
    lstr_vcfline = i.split("\t")
    editnuc = lstr_vcfline[4]
    chrom, position = lstr_vcfline[0], lstr_vcfline[1]
    bamposition = chrom + ':' + position + '-' + position

    lstr_outvcfline = lstr_vcfline

    #------------------
    # Run Samtools view
    #------------------ 
    # Run Samtools view on the BAM file with the given location

    cmd = "samtools view {} {}".format(bamFile, bamposition)
    sam_output = subprocess.check_output(cmd, shell=True, encoding='utf8')
    
    # separate the Samtools view output by lines (rstrip to remove last \n)
    sam_output = sam_output.rstrip().split("\n")

    for line in sam_output:

        if not re.match("\w", line):
            continue # sometimes get a blank line for some reason

        bamfields = line.split("\t")
        
        if len(bamfields) < 11:
            print("-warning: line[{}] has insufficient fields... skipping".format(line), file=sys.stderr)
            continue
        
        # separate the output
        samflag, readstart, cigar, sequencebases, qualscores = bamfields[1], bamfields[3], bamfields[5], bamfields[9], bamfields[10]

        if check_duplicate_marked(samflag):
            duplicateMarked += 1
            continue # not evaluating duplicate-marked reads.
        
        # get the current position 
        currentpos, readpos = int(readstart), 1
        base_readpos = []
        
        # letters 
        cigarletters = re.findall(r'[MIDNSHP=X]',cigar)
        # numbers
        cigarnums = re.split('[MIDNSHP]', cigar)
        
        '''
        CIGAR letters:
            I = insertion into the reference
            S = Soft clipping, clipped sequence is SEQ
            D = Deletion from reference 
            N = Skipped region from reference 
            M = Alignment Match, can be sequence mismatch or match 
        '''
        
        for k in range(len(cigarletters)):
            position = int(position)
            if currentpos > position:
                break
            if cigarletters[k] == "I" or cigarletters[k] == "S":
                readpos = readpos + int(cigarnums[k])

            elif cigarletters[k] == "D" or cigarletters[k] == "N": 
                currentpos = currentpos + int(cigarnums[k])
            
            elif cigarletters[k] == "M":
                for j in range(int(cigarnums[k])):
                    if currentpos == position:
                        base_readpos = readpos
                    currentpos += 1
                    readpos += 1

        read_end_pos = readpos
        
        if base_readpos:
            #----------------------------------------------------
            # check is within 6 bases of the ends and quality score 
            #----------------------------------------------------
            # If the base position is within the 6 base pairs of either side of the sequence -> Pass
            if (int(base_readpos) > 6) and (base_readpos < read_end_pos - 5):
                # If quality score is greater than or equal to the cutoff  --> PASS

                if ord(str(qualscores)[base_readpos-1]) >= minimal_base_quality + quality_offset:
                    ## counting as a PASS read, contributes to total coverage (newcov)

                    newcov += 1

                    ## see if it's a multi-mapped read  NH:i:(x)  where (x) > 1
                    m = re.search("\tNH:i:(\d+)", line)
                    if m:
                        num_mappings = int(m.group(1))
                        if num_mappings > 1:
                            multimappedReadCount += 1
                    
                    # check if the read base is the variant
                    if (sequencebases[base_readpos-1] == editnuc):
                        newmismatch+=1
                else:
                    ## Not a PASS read
                    basequalFail=1
            else:
                readPosFail=1

    #-----------------------
    # output lines 
    #-----------------------
    # VPR : Variant Passed Reads, reads that PASS filtering that contain the variation 
    # TPR : Total Passed Reads , reads that PASS filtering 
    # TDM : Total Duplicate Marked, number of reads that are duplicate marked
    lstr_outvcfline[7] += ";VPR={}".format(newmismatch)
    lstr_outvcfline[7] += ";TPR={}".format(newcov)
    lstr_outvcfline[7] += ";TDM={}".format(duplicateMarked)
    lstr_outvcfline[7] += ";VAF={:0.3f}".format(newmismatch/newcov if newcov > 0 else 0)
    lstr_outvcfline[7] += ";TMMR={}".format(multimappedReadCount)
    lstr_outvcfline[7] += ";MMF={:0.3f}".format(multimappedReadCount/newcov if newcov > 0 else 0)
    
    # variant frequency if needed 
    # varfreq = (newmismatch/newcov)

    # newcov        : total number of PASS reads 
    # newmismatch   : number of PASS's that support the variant 
    # duplicateMarked : number of duplicate-marked reads 
    new_line = "\t".join(lstr_outvcfline)

    return(new_line)



# open gzipped or regular text vcf files
def open_file_for_reading(filename):
    if re.search("\.gz$", filename):
        return gzip.open(filename, 'rt') # t needed for python3 to look like regular text
    else:
        return open(filename, 'r') # regular text file
    



def createAnnotations(vcf, bamFile, cpu, output_vcf_filename):
    logging.info("Filtering out mismatches in first 6 bp of reads")

    ##############
    # Constants #
    ##############
    # counters 
    not_filtered, removed= 0, 0

    #-----------------------
    # Open and create files 
    #-----------------------
    # Open input and output Files 
    infile = open_file_for_reading(vcf)
    # create the output file 

    outfile = open(output_vcf_filename, "w")

    #-----------
    # VCF Header
    #-----------
    # counter to count number of variants 

    variant_count = 0

    # Process the head of the VCF first and print it to the new file
    for line in infile:
        line = str(line)
        #print(line)
        if line[0] == "#":
            processVCFHead(line, outfile)
        else:
            variant_count +=1
    infile.close()

    #-----------
    # VCF Variants 
    #-----------
    # now need to reopen the VCF file 
    infile = open_file_for_reading(vcf)
    
    # set up the parallelization 
    # pool = mp.Pool(mp.cpu_count())

    pool = mp.Pool(cpu)

    # Now process each variant in the VCF and print the variants that pass to the output VCF product 
    # results = [step3_processing(i, bamFile) for i in infile.readlines() if i[0][0] != "#"]

    vcf_lines = infile.readlines()

    results = [pool.apply(evaluate_PASS_reads, args=(vcfline, bamFile)) for vcfline in vcf_lines if vcfline[0] != "#"]
    
    # CHECK: check to ensure that the number of variants given in the input VCF equals the number of variants in the output VCF
    if variant_count != len(results):
        "The output VCF has a different number of variants than the input VCF"

    # resort records since they might now be out of order due to multithreading

    def chr_pos_retriever(result_line): # inner function for sorting vcf output by chr, pos
        vals = result_line.split("\t")
        chr_val = vals[0]
        pos_val = int(vals[1])

        if len(chr_val) < 5:
            chr_val = chr_val.replace("chr", "chr0") # ensure chr08 comes before chr12, etc.
        
        return(chr_val, pos_val)

    # sort it
    results = sorted(results, key=chr_pos_retriever)

    # output results
    for i in results:
        outfile.write(i)

    pool.close()
    
    infile.close()
    outfile.close()



def main():
    ###########################################
    # Gather the arguments passed to the SNPiR script in command line 
    ###########################################

    ## Input Arguments
    # Description 
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Add annotations to a given VCF file.\n")
    # mandatory arguments 
    parser.add_argument('--vcf', required = True,
        help="input vcf.")
    parser.add_argument('--bam', required = True,
            help="input BAM file.")
    parser.add_argument("--output_vcf", required=True, 
        help="output vcf file containing annotations")
    parser.add_argument("--threads", required=False, type = int, default = 4, 
        help="Number of threads to run on.")
    
    # Parse the given arguments 
    args = parser.parse_args()

    vcf = args.vcf
    bam = args.bam
    cpu = args.threads
    output_vcf_filename = args.output_vcf

    createAnnotations(vcf, bam, cpu, output_vcf_filename)

    sys.exit(0)

if __name__ == "__main__":

    main()
