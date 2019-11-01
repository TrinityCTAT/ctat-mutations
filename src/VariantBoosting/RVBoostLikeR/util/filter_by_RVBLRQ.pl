#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "\n\tusage: $0 rvblrq.vcf [minQscore=0.05]\n\n";

my $vcf_file = $ARGV[0] or die $usage;
my $min_qscore = $ARGV[1];

unless (defined $min_qscore) {
    $min_qscore = 0.05;
}


main: {
    
    open(my $fh, $vcf_file) or die "Error, cannot open file $vcf_file";
    
    while(<$fh>) {
        if (/^\#/) {
            print;
            next;
        }
        /RVBLRQ=([\d\.]+)/ or die "Error, no RVBLRQ found for $_";
        my $qscore = $1;
        if ($qscore >= $min_qscore) {
            print;
        }
    }

    exit(0);
}


        
        
    
