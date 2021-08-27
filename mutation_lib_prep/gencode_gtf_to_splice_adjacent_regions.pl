#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GTF_utils;

use Carp;

# gencode_gtf_to_splice_adjacent_regions.pl $CTAT_GENOME_LIB/ref_annot.gtf > ref_annot.splice_adj.bed

my $usage = "usage: $0 file.GTF > file.splice_adjacent.bed\n\n";

my $gtf_file = $ARGV[0] or die $usage;
my $exon_splice_dist = 10;

main: {
            

    my %exon_splice_regions;
    

    my $gene_obj_indexer = {};

    print STDERR "-parsing GTF file: $gtf_file\n";
    my $asmbl_id_to_gene_list_href = &GTF_utils::index_GTF_gene_objs_from_GTF($gtf_file, $gene_obj_indexer);

    foreach my $asmbl_id (sort keys %$asmbl_id_to_gene_list_href) {
        
        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};
        
        #print "ASMBL: $asmbl_id, gene_ids: @gene_ids\n";
        
        foreach my $gene_id (@gene_ids) {
            
            my $gene_obj_ref = $gene_obj_indexer->{$gene_id} or die "Error, no gene obj for $gene_id";
            
            
            foreach my $gene_obj ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {

                my @intron_coordinates = $gene_obj_ref->get_intron_coordinates();
                
                foreach my $coordset (@intron_coordinates) {
                    my ($intron_lend, $intron_rend) = sort {$a<=>$b} @$coordset;
                    
                    my $intron_span_left = join(":", $asmbl_id, "L", $intron_lend, $intron_lend + $exon_splice_dist -1);
                    
                    my $intron_span_right = join(":", $asmbl_id, "R", $intron_rend - $exon_splice_dist + 1, $intron_rend);
                    
                    $exon_splice_regions{$intron_span_left} = { name => $intron_span_left,
                                                                asmbl_id => $asmbl_id,
                                                                lend => $intron_lend,
                                                                rend => $intron_lend + $exon_splice_dist -1,
                    };
                    
                    $exon_splice_regions{$intron_span_right} = { name => $intron_span_right,
                                                                 asmbl_id => $asmbl_id,
                                                                 lend => $intron_rend - $exon_splice_dist + 1,
                                                                 rend => $intron_rend,
                    };
                }
            }
        }
    }
    
    my @entries = values %exon_splice_regions;
    
    @entries = sort {$a->{asmbl_id} cmp $b->{asmbl_id}
                     ||
                         $a->{lend} <=> $b->{lend} } @entries;


    foreach my $entry (@entries) {
        print join("\t", $entry->{asmbl_id}, $entry->{lend}, $entry->{rend}, $entry->{name}) . "\n";
    }

    
    exit(0);
}

