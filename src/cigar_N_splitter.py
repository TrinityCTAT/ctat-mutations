#!/usr/bin/env python

import sys, os, re
import pysam


DEBUG = False

def main():
    
    usage = "\n\n\tusage: {} input.bam output.bam [DEBUG_flag]\n\n".format(sys.argv[0])

    if len(sys.argv) < 3:
        exit(usage)

    input_bam = sys.argv[1]
    output_bam = sys.argv[2]

    global DEBUG
    if len(sys.argv) > 3:
        DEBUG = True

    
    bam_reader = pysam.AlignmentFile(input_bam, "rb")
    bam_writer = pysam.AlignmentFile(output_bam, "wb", template=bam_reader)


    for read in bam_reader.fetch():
        if read.query_sequence is None:
            continue
        
        for split_read in split_cigar_N(read):
            bam_writer.write(split_read)


    bam_writer.close()
            
    sys.exit(0)




"""
M   BAM_CMATCH       0
I   BAM_CINS         1
D   BAM_CDEL         2
N   BAM_CREF_SKIP    3
S   BAM_CSOFT_CLIP   4
H   BAM_CHARD_CLIP   5
P   BAM_CPAD         6
=   BAM_CEQUAL       7
X   BAM_CDIFF        8
B   BAM_CBACK        9
"""
    

def split_cigar_N(read):


    if DEBUG:
        print("SPLITING_READ: {}".format(str(read)))
    
    cigar_tuples = read.cigartuples

    read_start = 0
    genome_start = read.reference_start -1
    read_pos = 0
    genome_pos = genome_start 
    cigar_tuples_include = list()

    read_seq = read.query_sequence
    try:
        read_quals = read.query_qualities
    except:
        read_quals = None
        
    #print("read_seq: {}".format(read_seq))
    
    split_reads = list()

    
    # initial soft or hard clipping:
    if cigar_tuples[0][0] in (4,5):
        if cigar_tuples[0][0] == 4:
            # soft-clipped
            read_start += cigar_tuples[0][1]
            read_pos += cigar_tuples[0][1]

        cigar_tuples.pop(0)

        
    for cigar_tuple in cigar_tuples:
        cigar_code, segment_length = cigar_tuple
        # M, =, X
        if cigar_code in (0, 7, 8):
            read_pos += segment_length 
            genome_pos += segment_length 
            cigar_tuples_include.append(cigar_tuple)

        # S, H at end
        elif cigar_code in (4,5):
            break

        # insertion
        elif cigar_code in (1,6):
            read_pos += segment_length
            cigar_tuples_include.append(cigar_tuple)

        # deletion
        elif cigar_code in (2, 9):
            genome_pos += segment_length
            cigar_tuples_include.append(cigar_tuple)

        # intron
        elif cigar_code == 3:
            # make a read alignment block
            split_read = make_split_read(read, read_seq, read_quals, read_start, read_pos, genome_start, genome_pos, cigar_tuples_include)
            split_reads.append(split_read)
            read_start = read_pos
            genome_start = genome_pos + segment_length
            genome_pos = genome_start
            cigar_tuples_include.clear()

    # get last segment
    if len(split_reads) > 0:
        split_read = make_split_read(read, read_seq, read_quals, read_start, read_pos, genome_start, genome_pos, cigar_tuples_include)
        split_reads.append(split_read)

    else:
        # no splitting, return original read
        split_reads = (read,)


    return split_reads

def make_split_read(read, read_seq, read_quals, read_start, read_pos, genome_start, genome_pos, cigar_tuples_include):

    read_start += 1
    genome_start += 1

    if DEBUG:
        print(f"Making split read: {read_start}-{read_pos} and {genome_start}-{genome_pos}")

    split_read_seq = read_seq[read_start-1:read_pos]

    if DEBUG:
        print("query_read_name: {}, split_read_seq: {}".format(read.query_name, split_read_seq))

    split_read_quals = None
    if read_quals is not None:
        split_read_quals = read_quals[read_start-1:read_pos] 

        if DEBUG:
            print("query_read_name: {}, query_read_quals: {}".format(read.query_name, split_read_quals))

        assert len(split_read_seq) == len(split_read_quals), "Error, split read seq: {} differs in length from read quals: {}".format(split_read_seq, split_read_quals)

    if DEBUG:
        print("query len: {}, cigars: {}".format(len(split_read_seq), cigar_tuples_include)) 
    

    
    a = pysam.AlignedSegment()
    a.query_name = read.query_name
    a.query_sequence = split_read_seq
    a.flag = read.flag
    a.reference_id = read.reference_id
    a.reference_start = genome_start
    a.mapping_quality = read.mapping_quality
    a.cigar = cigar_tuples_include #((0,10), (2,1), (0,25))
    a.next_reference_id = read.next_reference_id
    a.next_reference_start= read.next_reference_start
    a.template_length= read.template_length

    if split_read_quals is not None:
        a.query_qualities = split_read_quals   #pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
    #a.tags = (("NM", 1),
    #          ("RG", "L1"))
    
    return a




if __name__=='__main__':
    main()
