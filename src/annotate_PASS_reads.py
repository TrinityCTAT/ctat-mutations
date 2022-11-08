#!/usr/bin/env python3

import os, sys, re
import logging
import argparse
import subprocess
import gzip
import multiprocessing
import time
import numpy as np

## Set up the logging  
logging.basicConfig(format='\n %(levelname)s : %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)


def processVCFHead(line, outfile):
    '''
    Add the annotation flag description to the head of the VCF
        VPR : Variant Passed Reads, reads that PASS filtering that contain the variation
        TPR : Total Passed Reads , reads that PASS filtering, position not in the terminal 6 bases of read
        TCR: Total Covered Reads, meet quality requirements, variant position anywhere in read.
        PctExtPos: Fraction of variant-supporting reads with variant positions in the first six bases of the reads
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
            outfile.write("##INFO=<ID=TCR,Number=1,Type=Integer,Description=\"TCR: Total Covered Reads, meet quality requirements, variant position anywhere in read.\">\n")
            outfile.write("##INFO=<ID=PctExtPos,Number=1,Type=Integer,Description=\"Fraction of variant-supporting reads with variant positions in the first six bases of the reads\">\n")
            outfile.write("##INFO=<ID=TDM,Number=1,Type=Integer,Description=\"Total Duplicate Marked, number of reads that are duplicate marked \">\n")
            outfile.write("##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele fraction (VPR / TPR) \">\n")
            outfile.write("##INFO=<ID=TMMR,Number=1,Type=Integer,Description=\"Total multi-mapped reads at site\">\n")
            outfile.write("##INFO=<ID=MMF,Number=1,Type=Float,Description=\"Multi-mapped read fraction (TMMR / TPR) \">\n")
            outfile.write("##INFO=<ID=VMMR,Number=1,Type=Float,Description=\"Total variant-supporting multi-mapped reads \">\n")
            outfile.write("##INFO=<ID=VMMF,Number=1,Type=Float,Description=\"Variant-supporting multi-mapped read fraction (VMMR / VPR) \">\n")

        outfile.write(line) #writes the CHROM line


def progress_bar(progress_percent):
    #-------------------------
    # Create the Progress bar
    #-------------------------
    # Create progress bar to monitor the progress of the multiprocessing
    ## Remove the line from before
    sys.stdout.write("\r")
    sys.stdout.flush()
    ## print the progress bar and percent
    if progress_percent == 100:
        sys.stdout.write("[{}{}]{}".format("*" * progress_percent, " "* (100-progress_percent), str(progress_percent)+"%\n"))
    else:
        sys.stdout.write("[{}{}]{}".format("*" * progress_percent, " "* (100-progress_percent), str(progress_percent)+"%"))





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


def evaluate_PASS_reads(vcf_lines, bamFile):

    try:
        return [worker_evaluate_PASS_reads(i.decode('ASCII'), bamFile) for i in vcf_lines]
    except Exception as e:
        traceback.print_exc()
        return("ERROR: " + str(e))



def worker_evaluate_PASS_reads(vcf_line, bamFile):

    vals = vcf_line.split("\t")

    #--------------
    # Constants
    #--------------
    quality_offset = 33
    minimal_base_quality = 25
    minimum_mismatch = 1


    total_covered_reads = 0  # meets quality criteria at pos.
    total_pass_reads = 0 # not in the terminal 6 bases of read
    total_pass_variant_reads = 0 #
    total_fail_variant_reads = 0 # in the terminal region.

    total_pass_multimapped_reads = 0
    total_pass_multimapped_variant_reads = 0


    total_duplicate_marked = 0

    # process the input line
    lstr_vcfline = vcf_line.split("\t")
    editnuc = lstr_vcfline[4]
    chrom, position = lstr_vcfline[0], lstr_vcfline[1]
    bamposition = chrom + ':' + position + '-' + position

    lstr_outvcfline = lstr_vcfline

    #------------------
    # Run Samtools view
    #------------------
    # Run Samtools view on the BAM file with the given location

    cmd = "samtools view {} {}".format(bamFile, bamposition)
    sam_output = subprocess.check_output(cmd, shell=True).decode()


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
            total_duplicate_marked += 1
            continue # not evaluating duplicate-marked reads.

        # get the current position
        currentpos, readpos = int(readstart), 1
        base_readpos = None

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

            if ord(str(qualscores)[base_readpos-1]) >= minimal_base_quality + quality_offset:
                ## counting as a PASS read, contributes to total coverage (newcov)

                total_covered_reads += 1

                # ----------------------------------------------------
                # check is within 6 bases of the ends and quality score
                # ----------------------------------------------------
                # If the base position is within the 6 base pairs of either side of the sequence -> Pass
                pass_central_read_pos =  (base_readpos > 6) and (base_readpos < read_end_pos - 5)
                # If quality score is greater than or equal to the cutoff  --> PASS

                if pass_central_read_pos:
                    # central pos, min qual base => PASS status
                    total_pass_reads += 1

                ## see if it's a multi-mapped read  NH:i:(x)  where (x) > 1
                is_multimapped_read = False

                m = re.search("\tNH:i:(\d+)", line)
                if m:
                    num_mappings = int(m.group(1))
                    if num_mappings > 1:
                        is_multimapped_read = True


                if pass_central_read_pos and is_multimapped_read:
                    total_pass_multimapped_reads += 1


                # check if the read base is the variant
                is_variant_containing_read = (sequencebases[base_readpos-1] == editnuc)
                if is_variant_containing_read:
                    if pass_central_read_pos:
                        total_pass_variant_reads += 1
                        if is_multimapped_read:
                            total_pass_multimapped_variant_reads += 1
                    else:
                        # in read terminus
                        total_fail_variant_reads += 1


    #-----------------------
    # output lines
    #-----------------------
    # VPR : Variant Passed Reads, reads that PASS filtering that contain the variation
    # TPR : Total Passed Reads , reads that PASS filtering
    # TDM : Total Duplicate Marked, number of reads that are duplicate marked


    lstr_outvcfline[7] += ";VPR={}".format(total_pass_variant_reads)
    lstr_outvcfline[7] += ";TPR={}".format(total_pass_reads)
    lstr_outvcfline[7] += ";TCR={}".format(total_covered_reads)
    lstr_outvcfline[7] += ";PctExtPos={:0.3f}".format( (total_fail_variant_reads / (total_fail_variant_reads + total_pass_variant_reads)) if (total_fail_variant_reads + total_pass_variant_reads) > 0 else 0)
    lstr_outvcfline[7] += ";TDM={}".format(total_duplicate_marked)
    lstr_outvcfline[7] += ";VAF={:0.3f}".format(total_pass_variant_reads/total_pass_reads if total_pass_reads > 0 else 0)
    lstr_outvcfline[7] += ";TMMR={}".format(total_pass_multimapped_reads)
    lstr_outvcfline[7] += ";VMMR={}".format(total_pass_multimapped_variant_reads)
    lstr_outvcfline[7] += ";MMF={:0.3f}".format(total_pass_multimapped_reads/total_pass_reads if total_pass_reads > 0 else 0)
    lstr_outvcfline[7] += ";VMMF={:0.3f}".format(total_pass_multimapped_variant_reads/total_pass_variant_reads  if total_pass_variant_reads > 0 else 0)

    # variant frequency if needed
    # varfreq = (newmismatch/newcov)

    # newcov        : total number of PASS reads
    # newmismatch   : number of PASS's that support the variant
    # duplicateMarked : number of duplicate-marked reads
    new_line = "\t".join(lstr_outvcfline)

    return(new_line)




class SplitVCF:
    '''
    Class to annotate a given vcf in a multiprocessing fashion 
    '''
    # initialize object
    def __init__(   self, 
                    input_vcf,
                    cpu,
                    bamFile,
                    chunks,
                    output_vcf
                    ): # arguments to class instantiation 
        
        self.input_vcf    = input_vcf
        self.cpu          = cpu
        self.bamFile      = bamFile
        self.chunks       = chunks
        self.output_vcf   = output_vcf

        # Print the inputs for reference 
        message_str = f"""
            {'Annotating VCF':^20s}\n\t{'VCF':11s} : {input_vcf}\n\t{'BAM':11s} : {bamFile}\n\t{'CPU count':11s} : {cpu}\n\t{'Chunking':11s} : {chunks}"""
        logger.info(message_str)

    def getIDs( self ):
        #####################
        # Get contig IDs
        #####################
        # Get the contig IDs that are present in the VCF file 

        message_str = f"\t\tGetting contig IDs."
        logger.info(message_str)

        # Get contig IDs from the VCF header 
        tmp = "#"
        vcf_head = []
        with gzip.open(self.input_vcf, 'rb') as vcf:
            while tmp == "#":
                line = next(vcf).decode('ASCII')
                tmp = line[0]
                vcf_head.append(line)

        r = re.compile("##contig=<ID=*")
        long_IDS = list(filter(r.match, vcf_head))

        IDs = [i.split(",")[0].split("##contig=<ID=")[1] for i in long_IDS]

        self.IDs  = IDs
        return self 

    def getStats( self ):
        test = 0
        with gzip.open(self.input_vcf, 'rb') as vcf:
            for line in vcf:
                test+=1
        self.stats = test
        return self 


    def AddAnnotaion( self ):
        # Apply annotaions to the VCF
        ## Split the VCF file by thier contig
        message_str = f"\tAdding Annotaions"
        logger.info(message_str)

        # Get the index values for the vcf rows that will be subsetted 
        ##  the VCf is subsetted based on the chunks wanted 
        idx_range = list(range(len(self.header), self.stats))
        idx_list = np.array_split(idx_range, self.chunks)
        
        results = []
        def logging_return(line):
            results.append(line)

        # Initiate the Pool
        pool = multiprocessing.Pool(self.cpu)
        
        
    
        # get the start time and print it for reference
        start_time = time.process_time()
        message_str = f"\t\tStart Time: {time.asctime( time.localtime(time.time()) )}"
        logger.info(message_str)


        # Loop over all the chunks, and send them to the CPUs
        for j,i in enumerate(idx_list):
            with gzip.open(self.input_vcf, 'r') as vcf:
                idx = list(i)
                # read in the specific lines 
                output_vcf_lines = vcf.readlines()[idx[0]:(idx[-1]+1)]
                pool.apply_async(evaluate_PASS_reads, args=(output_vcf_lines, self.bamFile), callback = logging_return)

        pool.close()

        #~~~~~~~~~~~~~~~~~~
        # Progress Bar
        #~~~~~~~~~~~~~~~~~~
        prev_results_len = 0;
        while len(results) < len(idx_list):
            
            curr_results_len = len(results)
            if curr_results_len > prev_results_len:
                ## check for errors
                for idx in range(prev_results_len, curr_results_len):
                    for str_line in results[idx]:
                        if str_line.startswith("ERROR: "):
                            # error detected
                            raise results[idx]
                prev_results_len = curr_results_len

            # get the total progress made as a percentage
            progress_percent = int(len(results)/len(idx_list) * 100)
            # Print the progress percentage to the terminal
            progress_bar(progress_percent)

            time.sleep(2)
        # make sure to finish the progress bar with 100%
        progress_bar(100)


        # Resutls are in a list of lists, so need to flatten (join list of lists )
        flaten_list = [item for sublist in results for item in sublist]
        self.results = flaten_list

        #~~~~~~~~~
        # CHECK: 
        #   check to ensure that the number of variants given in the input VCF equals the number of variants in the output VCF
        #~~~~~~~~~
        variant_count = list(range(len(self.header), self.stats))
        # len(variant_count)
        if len(variant_count) != len(self.results):
            message_str = f"The output VCF has a different number of variants than the input VCF \n\t Actual: {len(variant_count)} \n\t Given: {len(self.results)} : "
        logger.info(message_str)

        return self 


    def writeOutput( self ):


        # resort records since they might now be out of order due to multithreading
        def chr_pos_retriever(result_line): # inner function for sorting vcf output by chr, pos
            vals = result_line.split("\t")
            chr_val = vals[0]
            pos_val = int(vals[1])

            if len(chr_val) < 5:
                chr_val = chr_val.replace("chr", "chr0") # ensure chr08 comes before chr12, etc.

            return(chr_val, pos_val)
        # sort it
        results = sorted(self.results, key=chr_pos_retriever)


        # Write to output file 
        message_str = f"\tWriting to output vcf:  {output_vcf}"
        logger.info(message_str)
        outfile = open(self.output_vcf, "w")

        #-----------
        # VCF Header
        #-----------
        # counter to count number of variants
        variant_count = 0
        for header_line in self.header:
            processVCFHead(header_line, outfile)

        #~~~~~~~~~~~~
        # VCF Body 
        #~~~~~~~~~~~~
        for i in results:
            outfile.write(i)

        # close the output file 
        outfile.close()




    def getHeader( self ):
        #--------------------------------------------------------
        # Get the header for the vcf file and hold it as a string
        #-------------------------------------------------------- 
        vcf_head = []
        tmp = "#"
        with gzip.open(self.input_vcf, 'rb') as vcf:
            while tmp == "#":
                line = next(vcf).decode('ASCII')
                tmp = line[0]
                if tmp == "#": # had trouble with my while loop
                    vcf_head.append(line)
        self.header = vcf_head 
        return self 

    




def main():

    ####################
    # Parse the use supplied information 
    ####################
    # Set the variables to supply 
    parser = argparse.ArgumentParser(description="Rsplit VCF.", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--vcf", type=str, required=True, help="VCF of interest.")
    parser.add_argument("--output_vcf", type=str, required=False, help="output directory.", default = ".")
    parser.add_argument("--threads", type=int, required=False, help="Number of CPUs to use.", default = "8")
    parser.add_argument("--bam", type=str, required=False, help="input bam file.")
    parser.add_argument("--chunks", type=int, required=False, help="Number to divide the VCF into.", default = "1000")

    # Parse the variables given 
    args = parser.parse_args()
    VCF = args.vcf
    output_vcf = args.output_vcf
    cpu = args.threads
    bam = args.bam
    chunks = args.chunks

    # if output_path ==  ".":
    #     output_path = os.getcwd()

    message_str = "\n####################################################################################\n\tAnnotating\n####################################################################################"
    print(message_str)

    ##############################
    # Load Data
    ##############################
    # initiate the ViFi object 
    
    VCF = SplitVCF(input_vcf = VCF, cpu = cpu, bamFile = bam, chunks = chunks, output_vcf = output_vcf)
    VCF = VCF.getIDs()
    VCF = VCF.getStats()
    VCF = VCF.getHeader()


    VCF = VCF.AddAnnotaion()
    VCF.writeOutput()

    sys.exit(0)

if __name__ == "__main__":

    main()