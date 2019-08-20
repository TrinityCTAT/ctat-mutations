#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "usage: $0 paired_reads_aligned.bam\n\n";

my $paired_bam_filename = $ARGV[0] or die $usage;


main: {

    ## flags we care about for identifying and cleansing pair info.
    my $IS_PAIRED_MASK = 0x0001;
    my $IS_PROPPER_PAIR_MASK = 0x0002;
    my $IS_MATE_UNMAPPED_MASK = 0x0008;
    my $STRAND_OF_MATE_MASK = 0x0020;
    my $IS_FIRST_IN_PAIR_MASK = 0x0040;
    my $IS_SECOND_IN_PAIR_MASK = 0x0080;
    
    
    my $cmd = "samtools view -h $paired_bam_filename | ";
    open(my $fh, $cmd) or die $!;
    while(my $line = <$fh>) {
        if ($line =~ /^\@/) {
            # header line
            print $line;
            next;
        }
        
        my @x = split(/\t/, $line);
        my $flag = $x[1];
        
        if (!  $flag & $IS_PAIRED_MASK ) {
            # single end.
            print $line;
            next;
            
        }

        my $read_name = $x[0];

        $flag &= ~$IS_PAIRED_MASK; # remove paired flag
        $flag &= ~$IS_PROPPER_PAIR_MASK; # remove this too.
        $flag &= ~$IS_MATE_UNMAPPED_MASK;
        $flag &= ~$STRAND_OF_MATE_MASK;
        
        # make single-end
        if ($flag & $IS_FIRST_IN_PAIR_MASK) {
            $read_name .= "/1";
            $flag &= ~$IS_FIRST_IN_PAIR_MASK; # remove bit setting.
        }
        elsif ($flag & $IS_SECOND_IN_PAIR_MASK) {
            $read_name .= "/2";
            $flag &= ~$IS_SECOND_IN_PAIR_MASK; # remove bit setting
        }
        else {
            die "Error, should be paired read, but the 1st or 2nd read pair bit flag is not set";
        }

        $x[0] = $read_name;
        $x[1] = $flag;
        
        $x[6] = "*";
        $x[7] = 0;
        $x[8] = 0;

        print join("\t", @x);
    }

    exit(0);
}


