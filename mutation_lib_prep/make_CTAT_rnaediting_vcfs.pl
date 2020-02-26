#!/usr/bin/env perl

use strict;
use warnings;

main: {
    
    open(my $ofh_37, ">GRCh37.RNAediting.vcf") or die $!;
    print $ofh_37 
        "##fileformat=VCFv4.2\n"
        . "##reference=GRCh37\n"
        . "##source=Rediportal\n"
        . "##INFO=<ID=RNAEDIT,Type=String,Description=\"A known or predicted RNA-editing site\">\n"
        . "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    

    open (my $ofh_38, ">GRCh38.RNAediting.vcf") or die $!;
    print $ofh_38 
        "##fileformat=VCFv4.2\n"
        . "##reference=GRCh38\n"
        . "##source=Rediportal\n"
        . "##INFO=<ID=RNAEDIT,Type=String,Description=\"A known or predicted RNA-editing site\">\n"
        . "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    
    print STDERR "-working...";
    
    open(my $fh, "TABLE1_hg19_hg38_download.txt") or die $!;
    my $header = <$fh>;
    while(<$fh>) {
        chomp;
        my @x = split(/\t/);
        
        my $chr_hg38 = $x[0];
        my $pos_hg38 = $x[1];
        
        my $chr_hg19 = $x[2];
        my $pos_hg19 = $x[3];

        my $ref_allele = $x[4];
        my $alt_allele = $x[5];
        my $orient = $x[6];

        my $db = $x[7];
        
        my $dblisting = &parse_db_listing($db);

        print $ofh_37 join("\t", $chr_hg19, $pos_hg19, ".", $ref_allele, $alt_allele, ".", ".", "RNAEDIT=$dblisting\n");

        print $ofh_38 join("\t", $chr_hg38, $pos_hg38, ".", $ref_allele, $alt_allele, ".", ".", "RNAEDIT=$dblisting\n");
        
    }

    close $fh;

    close $ofh_37;
    close $ofh_38;

    print STDERR "\ndone, now indexing.";


    foreach my $vcf_file ("GRCh37.RNAediting.vcf", "GRCh38.RNAediting.vcf") {

        &process_cmd("set -eou pipefail; cat $vcf_file | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | \"sort -k1,1 -k2,2n\"}' > tmp.sorted.vcf");
        &process_cmd("mv tmp.sorted.vcf $vcf_file");
        
        &process_cmd("bgzip -f $vcf_file");
        &process_cmd("bcftools index $vcf_file.gz");
    }
        

    exit(0);
    
}

sub process_cmd {
    my ($cmd) = @_;
    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
}

####
sub parse_db_listing {
    my ($db) = @_;
    
    # DARNED, RADAR, and ATLAS
    # in rediportal paper: Currently, RNA editing events are annotated in three main databases: DARNED (http://darned.ucc.ie/) (10), RADAR (http://rnaedit.com/) (11) and REDIdb (http: //srv00.recas.ba.infn.it/py script/REDIdb/)1"
    # "In contrast, REDIportal includes more than 4.5 millions of A-to-I changes obtained merging RNA editing positions from our Inosinome ATLAS (6)"

    my %db_conversion = ( 'R' => 'RADAR',
                          'A' => 'ATLAS',
                          'D' => 'DARNED');

    my @dbs;
    for my $dbchar (split(/,/, $db)) {
        
        push(@dbs, $db_conversion{$dbchar});
    }

    my $dblisting = join(",", sort @dbs);

    return($dblisting);
    
}


        
