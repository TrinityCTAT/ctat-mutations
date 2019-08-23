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
        TMM : Total Multi-Mapping, number of reads that are multi-mapped 
    '''
    if line[0] == "#":
        if re.match("#CHROM\t", line):
            # add header info line for the repeat annotation type
            # once you get to the #CHROM line, 
            #       write the new INFO annotations above the #CHROM header line 
            outfile.write("##INFO=<ID=VPR,Type=Integer,Description=\"Variant Passed Reads, reads that PASS filtering that contain the variation\">\n")
            outfile.write("##INFO=<ID=TPR,Type=Integer,Description=\"Total Passed Reads , reads that PASS filtering\">\n")
            outfile.write("##INFO=<ID=TMM,Type=Integer,Description=\"Total Multi-Mapping, number of reads that are multi-mapped \">\n")

        outfile.write(line)

def processSamFlag(sam_flag):
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
    
    # set the multi-mapping indicator to 0
    multi_mapping = 0 
    # multi_mapping = 1 if (alignment & 1024)
    # check the SAM flag, (10000000000 means this is multi-mapping)
    if len(binary_flag) >= 11:
        multimapping_test = binary_flag[-11] 
        if int(multimapping_test) == 1:
            multi_mapping = 1

    return multi_mapping
    

def step3_processing(i, bamFile):
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
    output =""
    output_failed =""
    
    # process the input line 
    line = i.split("\t")
    editnuc = line[4]
    chrom, position = line[0], line[1]
    bamposition = chrom + ':' + position + '-' + position

    #------------------
    # Run Samtools view
    #------------------ 
    # Run Samtools view on the BAM file with the given location

    cmd = "samtools view {} {}".format(bamFile, bamposition)
    sam_output = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8', shell=True).communicate()
    # separate the Samtools view output by lines (rstrip to remove last \n)
    sam_output = sam_output[0].rstrip().split("\n")
    
    multi_mapping = 0
    for j in sam_output:
        bamfields = j.split("\t")
        # separate the output
        alignment, readstart, cigar, sequencebases, qualscores = bamfields[1], bamfields[3], bamfields[5], bamfields[9], bamfields[10]

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
        
        if base_readpos:
            # Alignment is the SAM Flags
            multi_mapping += processSamFlag(alignment)

            
            #----------------------------------------------------
            # check is within 6 bases of the ends and quality score 
            #----------------------------------------------------
            # If the base position is within the 6 base pairs of either side of the sequence -> Pass
            if (int(base_readpos) > 6) and (base_readpos < readpos - 5):
                # If quality score is greater than or equal to the cutoff  --> PASS
                if ord(str(qualscores)[base_readpos-1]) >= minimal_base_quality + quality_offset:
                    newcov += 1
                    # check if the base is the edited nucleotide 
                    if (sequencebases[base_readpos-1] == editnuc):
                        newmismatch+=1
                else:
                    basequalFail=1
            else:
                readPosFail=1

    #-----------------------
    # output lines 
    #-----------------------
    # VPR : Variant Passed Reads, reads that PASS filtering that contain the variation 
    # TPR : Total Passed Reads , reads that PASS filtering 
    # TMM : Total Multi-Mapping, number of reads that are multi-mapped 
    line[7] += ";VPR={}".format(newmismatch)
    line[7] += ";TPR={}".format(newcov)
    line[7] += ";TMM={}".format(multi_mapping)

    # variant frequency if needed 
    # varfreq = (newmismatch/newcov)

    # newcov        : total number of PASS reads 
    # newmismatch   : number of PASS's that support the variant 
    # multi_mapping : number of multi-mapped reads 
    new_line = "\t".join(line) + "\n"
    return(new_line)



def createAnnotations(vcf, bamFile, cpu, output_path):
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
    infile = open(vcf, "r")
    # create the output file 
    outfile_path = "{}_annotated.txt".format(vcf)
    out_path = os.path.join(output_path, outfile_path)
    outfile = open(outfile_path, "w")


    #-----------
    # VCF Header
    #-----------
    # counter to count number of variants 
    variant_count = 0
    # Process the head of the VCF first and print it to the new file
    for i in infile.readlines():
        if i[0][0] == "#":
            processVCFHead(line = i , outfile = outfile)
        else:
            variant_count +=1
    infile.close()

    #-----------
    # VCF Variants 
    #-----------
    # now need to reopen the VCF file 
    infile = open(vcf, "r")
    
    # set up the parallelization 
    # pool = mp.Pool(mp.cpu_count())
    pool = mp.Pool(cpu)
    # Now process each variant in the VCF and print the variants that pass to the output VCF product 
    # results = [step3_processing(i, bamFile) for i in infile.readlines() if i[0][0] != "#"]
    results = [pool.apply(step3_processing, args=(i, bamFile)) for i in infile.readlines() if i[0][0] != "#"]

    # CHECK: check to ensure that the number of variants given in the input VCF equals the number of variants in the output VCF
    if variant_count != len(results):
        "The output VCF has a different number of variants than the input VCF"
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
    parser.add_argument("--output_dir", required=True, 
        help="output directory")
    parser.add_argument("--threads", required=False, type = int, default = 4, 
        help="Number of threads to run on.")
    
    # Parse the given arguments 
    args = parser.parse_args()

    vcf = args.vcf
    bam = args.bam
    cpu = args.threads
    output_path = args.output_dir

    createAnnotations(vcf, bam, cpu, output_path)

if __name__ == "__main__":

    main()
