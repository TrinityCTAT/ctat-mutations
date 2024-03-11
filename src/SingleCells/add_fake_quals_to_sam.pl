#!/usr/bin/env perl

use strict;
use warnings;


my $seen_header = 0;

while(my $line = <>) {
    unless ($seen_header) {
        unless ($line =~ /^\@/) {
            die "Error, no header observed ";
        }
        $seen_header = 1;
    }
    if ($line =~ /^\@/) { 
        print $line;
        next;
    }

    my @pts = split(/\t/, $line);
    if ($pts[9] eq "*") {
        next;
    }

    if ($pts[10] eq "*") {
        # adding fake qual values at Q30
        $pts[10] = "?" x length($pts[9]);
    }

    print join("\t", @pts);

}


exit(0);


        
