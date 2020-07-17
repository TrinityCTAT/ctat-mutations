#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(max);


my $ATT_LIST = "DJ,ReadPosRankSum,QD,FS,ED,VPR,VAF,VMMF,SPLICEADJ,RPT,Homopolymer,Entropy,RNAEDIT,RS,INDEL";

my $usage = "\n\n\tusage: $0 variants_annotated.vcf atts_list=$ATT_LIST\n\n";

my $vcf_file = $ARGV[0] or die $usage;
if ($ARGV[1]) {
    $ATT_LIST = $ARGV[1];
}

my @ATTRIBUTE_LIST = split(/,/, $ATT_LIST);

my %ATTRIBUTE_HASH = map { + $_ => 1 } @ATTRIBUTE_LIST;

main: {

    # print header for output
    print join("\t", "chrpos", @ATTRIBUTE_LIST) . "\n";

    my $fh;
    if ($vcf_file =~ /\.gz$/) {
        open($fh, "gunzip -c $vcf_file | ") or die "Error, cannot gunzip file $vcf_file";
    }
    else {
        open($fh, $vcf_file) or die "Error, cannot open file: $vcf_file";
    }
    
    while(<$fh>) {
        if (/^\#/) { next; }
        chomp;
        my @x = split(/\t/);
        
        my $info = $x[7];
        
        my @annotations = split(/;/, $info);
        my %annot_hash;
        foreach my $annot (@annotations) {
            my ($key, $val) = split(/=/, $annot, 2);
            $annot_hash{$key} = $val;
        }
        
        if ($ATTRIBUTE_HASH{INDEL}) {
            my $ref_base = $x[3];
            my $alt_bases = $x[4];
            $annot_hash{INDEL} = &get_indel_length($ref_base, $alt_bases);
        }
        
        my $chr_pos = join(":", $x[0], $x[1], $x[3], $x[4]);
        
        my @vals = ($chr_pos);
        foreach my $att (@ATTRIBUTE_LIST) {
            my $val = (exists $annot_hash{$att}) ? $annot_hash{$att} : "NA";
            push (@vals, $val);
        }
        print join("\t", @vals) . "\n";
    }
    close $fh;

    exit(0);
}


####
sub get_indel_length {
    my ($ref_base, $alt_bases) = @_;

    my @lengths;

    push (@lengths, length($ref_base));
    
    my @alt_bases_list = split(/,/, $alt_bases);
    foreach my $alt_base_entry (@alt_bases_list) {
        push (@lengths, length($alt_base_entry));
    }

    my $max_len = max(@lengths);

    return($max_len - 1);
}

