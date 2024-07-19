#!/usr/bin/env python3

import os, sys, re
import logging
import argparse
import subprocess
import gzip
import multiprocessing
import time
import numpy as np
import traceback

## Set up the logging

FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"

logging.basicConfig(level=logging.INFO, 
                    format=FORMAT,
                    datefmt='%H:%M:%S')


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


def evaluate_PASS_reads(vcf_lines, bamFile, sc_mode):

    try:
        #return [worker_evaluate_PASS_reads(i.decode('ASCII'), bamFile, sc_mode) for i in vcf_lines]
        return [worker_evaluate_PASS_reads(i.decode("utf-8"), bamFile, sc_mode) for i in vcf_lines]
    except Exception as e:
        traceback.print_exc()
        #return([["ERROR: " + str(e)]])
        raise e


def worker_evaluate_PASS_reads(vcf_line, bamFile, sc_mode):

    vals = vcf_line.split("\t")

    #--------------
    # Constants
    #--------------
    quality_offset = 33
    minimal_base_quality = 25

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

    ref_bases = lstr_vcfline[3]
    alt_bases = lstr_vcfline[4]

    if alt_bases == "*":
        # spanning deletion
        # the deletion would be tackled separately by a different vcf line with the specifics.
        # nothing to do here.
        return( [vcf_line, [], [] ] )

    variant_type = "M"
    if len(ref_bases) > len(alt_bases):
        variant_type = "D"
    elif len(ref_bases) < len(alt_bases):
        variant_type = "I"
    
    lstr_outvcfline = lstr_vcfline

    #------------------
    # Run Samtools view
    #------------------
    # Run Samtools view on the BAM file with the given location

    cmd = "samtools view {} {}".format(bamFile, bamposition)
    sam_output = subprocess.check_output(cmd, shell=True).decode()


    # separate the Samtools view output by lines (rstrip to remove last \n)
    sam_output = sam_output.rstrip().split("\n")

    reads_with_variant = []
    reads_without_variant = []
    
    for line in sam_output:
        if not re.match("\w", line):
            continue # sometimes get a blank line for some reason
        
        bamfields = line.split("\t")

        if len(bamfields) < 11:
            print("-warning: line[{}] has insufficient fields... skipping".format(line), file=sys.stderr)
            continue
        
        # separate the output
        readname, samflag, readstart, cigar, sequencebases, qualscores = bamfields[0], bamfields[1], bamfields[3], bamfields[5], bamfields[9], bamfields[10]

        readlen = len(sequencebases)

        if sequencebases == "*":
            continue

        if qualscores == "*":
            # single cell grouped deduped lack qual values - assume all hifi Q30, char '?'
            qualscores = "".join(["?" for x in range(len(sequencebases))])


        if check_duplicate_marked(samflag):
            total_duplicate_marked += 1
            if not sc_mode:
                continue # not evaluating duplicate-marked reads unless in single cell mode.

        logger.debug(f"\n\n//{vcf_line}")
        logger.debug("// readname: {}, readstart: {}, cigar: {}, sequencebases: {}".format(readname, readstart, cigar, sequencebases))
        
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
        
        position = int(position)
        num_cigarletters = len(cigarletters)
        num_cigarnums = len(cigarnums)

        adjacent_base_type = "M" #default
        
        for k in range(num_cigarletters):
            
            if currentpos > position:
                break

            if readpos > readlen:
                break

            try:

                cigarletter = cigarletters[k]
                cigarnumval = int(cigarnums[k])

                logger.debug("readpos: {}".format(readpos))
                
                if (cigarletter == "I") or (cigarletter == "S"):
                    readpos = readpos + cigarnumval
                    logger.debug("currentpos: {}, readpos: {}, cigar: {}, readbase: {}".format(currentpos, readpos, cigarletter, sequencebases[readpos-1])) 

                elif (cigarletter == "D") or (cigarletter == "N"):
                    currentpos = currentpos + cigarnumval
                    logger.debug("currentpos: {}, readpos: {}, cigar: {}, readbase: {}".format(currentpos, readpos, cigarletter, sequencebases[readpos-1])) 

                elif cigarletter == "M":
                    # iterate each base position of the match
                    num_aligned_bases = int(cigarnums[k])
                    segment_end_pos = readpos + num_aligned_bases -1
                    logger.debug("segment end pos: {}".format(segment_end_pos))
                    for j in range(num_aligned_bases):
                        logger.debug("currentpos: {}, readpos: {}, cigar: {}, readbase: {}".format(currentpos, readpos, cigarletter, sequencebases[readpos-1])) 
                        if currentpos == position:
                            base_readpos = readpos
                            logger.debug("** base_readpos: {} at variant site: {}".format(base_readpos, currentpos))
                            if base_readpos == segment_end_pos:
                                # see if its at an indel
                                logger.debug("k + 1 = {} and num cigar letters: {}".format(k+1, num_cigarletters))
                                if k + 1 < num_cigarletters:
                                    adjacent_base_type = cigarletters[k+1]
                                    logger.debug("adjacent_base_type: {}".format( adjacent_base_type))
                            #nuc = sequencebases[base_readpos-1]
                            #print("\t".join([readname, str(position), str(base_readpos), nuc]))

                        currentpos += 1
                        readpos += 1

                read_end_pos = readpos

                logger.debug("\t".join([ref_bases, alt_bases, adjacent_base_type, str(cigarletters)]))

            except IndexError:
                logger.error("Error processing vcf line: {}".format(vcf_line))
                raise IndexError
        

            
        if base_readpos:
            #print("Got read pos")
            base_qual = ord(str(qualscores)[base_readpos-1]) - quality_offset

            base_qual_ok = (base_qual >= minimal_base_quality)
            if base_qual_ok or sc_mode:
                ## counting as a PASS read, contributes to total coverage (newcov)

                if base_qual_ok:
                    total_covered_reads += 1

                # ----------------------------------------------------
                # check is within 6 bases of the ends and quality score
                # ----------------------------------------------------
                # If the base position is within the 6 base pairs of either side of the sequence -> Pass
                pass_central_read_pos =  (base_readpos > 6) and (base_readpos < read_end_pos - 5)
                # If quality score is greater than or equal to the cutoff  --> PASS

                if pass_central_read_pos and base_qual_ok:
                    # central pos, min qual base => PASS status
                    total_pass_reads += 1

                ## see if it's a multi-mapped read  NH:i:(x)  where (x) > 1
                is_multimapped_read = False

                m = re.search("\tNH:i:(\d+)", line)
                if m:
                    num_mappings = int(m.group(1))
                    if num_mappings > 1:
                        is_multimapped_read = True


                if pass_central_read_pos and is_multimapped_read and base_qual_ok:
                    total_pass_multimapped_reads += 1


                # check if the read base is the variant

                is_variant_containing_read = False
                is_refbases_containing_read = False
                
                if ( (re.match(alt_bases, sequencebases[(base_readpos-1):]) is not None)
                     and
                     variant_type == adjacent_base_type ):
                    is_variant_containing_read = True

                elif (re.match(ref_bases, sequencebases[(base_readpos-1):]) is not None):
                    is_refbases_containing_read = True
                
                if is_variant_containing_read:
                    reads_with_variant.append(readname) #capture variant-containing read for single cells
                    if base_qual_ok:
                        if pass_central_read_pos:
                            total_pass_variant_reads += 1
                            if is_multimapped_read:
                                total_pass_multimapped_variant_reads += 1
                        else:
                            # in read terminus
                            total_fail_variant_reads += 1

                elif is_refbases_containing_read: # not variant containing
                    reads_without_variant.append(readname) #capture non-variant read for single cells
                    

    #-----------------------
    # output lines
    #-----------------------
    # VPR : Variant Passed Reads, reads that PASS filtering that contain the variation
    # TPR : Total Passed Reads , reads that PASS filtering
    # TDM : Total Duplicate Marked, number of reads that are duplicate marked


    if not sc_mode:
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

    return( [new_line, reads_with_variant, reads_without_variant ] )




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
                    output_vcf,
                    sc_mode
                    ): # arguments to class instantiation 
        
        self.input_vcf    = input_vcf
        self.cpu          = cpu
        self.bamFile      = bamFile
        self.chunks       = chunks
        self.output_vcf   = output_vcf
        self.sc_mode = sc_mode

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


    def AddAnnotation( self ):
        # Apply annotaions to the VCF
        ## Split the VCF file by thier contig
        message_str = f"\tAdding Annotations"
        logger.info(message_str)

        # Get the index values for the vcf rows that will be subsetted 
        ##  the VCf is subsetted based on the chunks wanted 
        idx_range = list(range(len(self.header), self.stats))
        idx_list = np.array_split(idx_range, self.chunks)
        ## Remove empty arrays 
        idx_list = [i for i in idx_list if len(i) != 0]

        #print("idx_list: {}".format(idx_list))
        
        results = []
        def logging_return(result):
            results.append(result)



    
        # Initiate the Pool
        if self.cpu > 1:
            pool = multiprocessing.Pool(self.cpu)


        def error_handler(error):
            # note, without this, process will just hang forever
            logger.error("ERROR_HANDLER - CAUGHT: " + str(error))
            pool.terminate()
            pool.join()
            sys.exit(2)
            
                    
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

                if self.cpu > 1:
                    pool.apply_async(evaluate_PASS_reads, args=(output_vcf_lines, self.bamFile, self.sc_mode),
                                     callback = logging_return,
                                     error_callback = error_handler)
                else:
                    result = evaluate_PASS_reads(output_vcf_lines, self.bamFile, self.sc_mode)
                    logging_return(result)

        if self.cpu > 1:
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
                        chunk = results[idx]
                        for entry in chunk:
                            str_line = entry[0]
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


        # Results are in a list of lists, so need to flatten (join list of lists )
        ##
        ## looks like:
        ##   [
        ##     [ # chunk:
        ##        [line, var_reads_list, no_var_reads_list],
        ##        [line, var_reads_list, no_var_reads_list],
        ##        ...
        ##     ]
        ##  ]

        ##
        flattened_list = []

        for chunk in results:
            for line_n_reads_list in chunk:
                flattened_list.append(line_n_reads_list)

        # flaten_list = [item for sublist in results for item in sublist]
        
        self.results = flattened_list
        
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
            vals = result_line[0].split("\t")
            chr_val = vals[0]
            pos_val = int(vals[1])

            if len(chr_val) < 5:
                chr_val = chr_val.replace("chr", "chr0") # ensure chr08 comes before chr12, etc.

            return(chr_val, pos_val)

        # sort it
        results = sorted(self.results, key=chr_pos_retriever)


        # Write to output file 
        message_str = f"\tWriting to output vcf:  {self.output_vcf}"
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

        if self.sc_mode:
            var_reads_filename = self.output_vcf + ".sc_reads"
            sc_reads_ofh = open(var_reads_filename, "w")
            print("\t".join(["chr_pos_variant", "num_reads_with_variant", "reads_with_variant", "num_ref_matching_reads", "ref_matching_reads"]), file=sc_reads_ofh)

        
        for i in results:
            vcf_line = i[0]
            outfile.write(vcf_line)
            if self.sc_mode:
                reads_w_var_list = i[1]
                reads_wo_var_list = i[2]
                vcf_pts = vcf_line.split("\t")
                pos_token = ":".join([vcf_pts[0], vcf_pts[1], vcf_pts[3], vcf_pts[4] ])
                vcf_line = vcf_line.rstrip()
                print("\t".join([pos_token,
                                 str(len(reads_w_var_list)),
                                 ",".join(reads_w_var_list),
                                 str(len(reads_wo_var_list)),
                                 ",".join(reads_wo_var_list)]),
                      file=sc_reads_ofh)
                
        # close the output file 
        outfile.close()
        if self.sc_mode:
            sc_reads_ofh.close()




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
    parser.add_argument("--output_vcf", type=str, required=True, help="output vcf filename")
    parser.add_argument("--threads", type=int, required=False, help="Number of CPUs to use.", default = "8")
    parser.add_argument("--bam", type=str, required=True, help="input bam file.")
    parser.add_argument("--chunks", type=int, required=False, help="Number to divide the VCF into.", default = "1000")
    parser.add_argument("--sc_mode", action='store_true', default=False, help="single cell mode, so capture cell/variant info. Creates separate file --output_vcf + .sc_reads")
    parser.add_argument("--debug", "-d", action='store_true', default=False, help='debug mode, verbose')
    
    global DEBUG
    
    
    
    # Parse the variables given 
    args = parser.parse_args()
    VCF = args.vcf
    output_vcf = args.output_vcf
    cpu = args.threads
    bam = args.bam
    chunks = args.chunks
    sc_mode = args.sc_mode

    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    # if output_path ==  ".":
    #     output_path = os.getcwd()

    message_str = "\n####################################################################################\n\tAnnotating\n####################################################################################"
    print(message_str)

    ##############################
    # Load Data
    ##############################
    # initiate the ViFi object 
    
    VCF = SplitVCF(input_vcf = VCF, cpu = cpu, bamFile = bam, chunks = chunks, output_vcf = output_vcf, sc_mode=sc_mode)
    VCF = VCF.getIDs()
    VCF = VCF.getStats()
    VCF = VCF.getHeader()
    
    
    VCF = VCF.AddAnnotation()
    VCF.writeOutput()

    sys.exit(0)

if __name__ == "__main__":

    main()
