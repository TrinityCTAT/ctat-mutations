#!/usr/bin/env python3

import sys, os, re
import argparse
import pysam
import logging
import subprocess

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)



def main():

    parser = argparse.ArgumentParser(description="separate bam into strand-specific bam files", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--bam", type=str, required=True, help="input bam filename")

    parser.add_argument("--output_prefix", type=str, required=True, help="output prefix: files named ${output_prefix}.${strand}.bam")

    args = parser.parse_args()

    input_bam_filename = args.bam
    output_prefix = args.output_prefix
    
    bamfile_reader = pysam.AlignmentFile(input_bam_filename, "rb")

    top_strand_bam_filename = output_prefix + ".+.bam"
    bottom_strand_bam_filename = output_prefix + ".-.bam"
    
    top_strand_bamfile_writer = pysam.AlignmentFile(top_strand_bam_filename, "wb", template=bamfile_reader)
    bottom_strand_bamfile_writer = pysam.AlignmentFile(bottom_strand_bam_filename, "wb", template=bamfile_reader)

    output_bam_files = (top_strand_bam_filename, bottom_strand_bam_filename)

    num_records = 0
    num_forward = 0
    num_reverse = 0
    num_neither = 0
    for read in bamfile_reader:
        num_records += 1
        
        if read.is_forward:
            top_strand_bamfile_writer.write(read)
            num_forward += 1
            
        elif read.is_reverse:
            bottom_strand_bamfile_writer.write(read)
            num_reverse += 1
            
        else:
            num_neither += 1
    
    top_strand_bamfile_writer.close()
    bottom_strand_bamfile_writer.close()

    assert num_records > 0, "No records read from input bam file: {}".format(input_bam_filename)
    

    report = "\n".join(["Num input bam records: {}".format(num_records),
                        "Num top strand: {} = {:.1f}%".format(num_forward, num_forward/num_records * 100),
                        "Num bottom strand: {} = {:.1f}%".format(num_reverse, num_reverse/num_records * 100),
                        "Num neither strand and ignored: {} = {:.1f}%".format(num_neither, num_neither/num_records) ] )
    
    logger.info(report)

    for output_bam_file in output_bam_files:
        subprocess.check_call("samtools index {}".format(output_bam_file), shell=True)

    
    sys.exit(0)
     


if __name__=='__main__':
    main()
