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


def evaluate_PASS_reads(vcf_line, bamFile):

    try:
        return worker_evaluate_PASS_reads(vcf_line, bamFile)
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






def runAddAnnotator(input_vcf, idx, count, header, bamFile):
    '''
    Open the VCF file and extract the lines of interest
        write these lines to a new vcf (subsetted)
        run the expression analysis on this subsetted file 
        then remove the temp subsetted files 
    
    '''
    with gzip.open(input_vcf, 'r') as vcf:
        # read in the specific lines 
        output_vcf_lines = vcf.readlines()[idx[0]:(idx[-1]+1)]
        
        file_name = f"tmp{count}.vcf"
        with open(file_name, "w") as tmp:
            # Write the header to the file 
            for header_line in header:
                tmp.write(header_line)#.encode())
            # Write the output to the temp file 
            for line in output_vcf_lines:
                new_line = evaluate_PASS_reads(line.decode(),bamFile)
                tmp.write(new_line)

    


class SplitVCF:
    '''
    Class to annotate a given vcf in a multiprocessing fashion 
    '''
    # initialize object
    def __init__(   self, 
                    input_vcf,
                    cpu,
                    bamFile
                    ): # arguments to class instantiation 
        
        self.input_vcf    = input_vcf
        self.cpu          = cpu
        self.bamFile      = bamFile

        message_str = f"Annotating VCF:\n\tVCF : {input_vcf}\n\tCPU count : {cpu}\n\t BAM : {bamFile}"
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
        message_str = f"\tAdding Expression INFO"
        logger.info(message_str)

        idx_range = list(range(len(self.header), self.stats))
        idx_list = np.array_split(idx_range, self.cpu)

        # Initiate the Pool
        pool = multiprocessing.Pool(self.cpu)
        
        
    
        # get the start time
        start_time = time.process_time()
        message_str = f"\t\tStart Time: {time.asctime( time.localtime(time.time()) )}"
        logger.info(message_str)

        # get all the arguments into a list of lists to pass for multiprocessing 
        args = [[self.input_vcf, list(i), j, self.header, self.bamFile] for j,i in enumerate(idx_list) ]
        
        pool.starmap(runAddAnnotator, args)

        pool.close()

        # get the end time
        end_time = time.process_time()
        message_str = f"\t\tEnd Time: {time.asctime( time.localtime(time.time()) )}"
        logger.info(message_str)
        # get execution time
        time_total = end_time - start_time
        message_str = f"\t\tProcess time: {time_total}"
        logger.info(message_str)

        # save the arguments 
        file_list = []
        for count,i in enumerate(idx_list):
            file_list.append(f"tmp{count}.vcf")

        self.tmp_files = file_list

        return self 


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

    

    def mergeVCFs( self, file_list, output):
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        # Merge the VCFs created in expression annotaion step 
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        message_str = f"\t\tMerging VCFs."
        logger.info(message_str)

        vcfs = [f"I={i}" for i in file_list ]

        vcfs = " ".join(vcfs)

        cmd = f"""java -jar /usr/local/src/picard.jar MergeVcfs {vcfs} O={output}"""

        print(cmd.replace(" ", " \n"))

        subprocess.run(cmd, shell=True)


    def removeVCFs( self, file_list):
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        # Remove the unneeded VCFs created in expression annotaion step 
        #~~~~~~~~~~~~~~~~~~~~~~~~~
        message_str = f"\tRemoving unneeded VCFs."
        logger.info(message_str)

        # remove the intermediate files 
        vcfs_remove = " ".join(file_list)

        cmd = f"""rm {vcfs_remove}"""

        subprocess.run(cmd, shell=True)




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

    # Parse the variables given 
    args = parser.parse_args()
    VCF = args.vcf
    output_vcf = args.output_vcf
    cpu = args.threads
    bam = args.bam

    if output_path ==  ".":
        output_path = os.getcwd()

    message_str = "\n####################################################################################\n\tAnnotating Expression Information\n####################################################################################"
    print(message_str)

    ##############################
    # Load Data
    ##############################
    # initiate the ViFi object 
    
    VCF = SplitVCF(input_vcf = VCF, cpu = cpu, bamFile = bam)
    VCF = VCF.getIDs()
    VCF = VCF.getStats()
    VCF = VCF.getHeader()


    VCF.AddAnnotaion()
    VCF.mergeVCFs(file_list = VCF.tmp_files, output = output_vcf)
    # VCF.removeVCFs(file_list = VCF.tmp_files)

    sys.exit(0)

if __name__ == "__main__":

    main()