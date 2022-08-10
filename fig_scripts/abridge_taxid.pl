#!/usr/bin/perl -w

### This program will take as input a large taxon id file
### and output a much smaller one with genera, families, and orders only

use strict;
use warnings;

### from this read name / taxon name file
my $taxa_file = "../fig_data/taxonomy.tab";

open INPUT, $taxa_file || die "Can't open taxa file\n";
open OUTPUT, ">abridged_taxonomy.tab" || die "Can't open taxa file\n";

while (<INPUT>) {
	chomp;
	if(my ($node, $level, $name) = /(\d+)\t\d+\t\d+\t(no\srank|phylum|class|order|family|genus)\t(.+)/) {
 		print OUTPUT "$node\t$level\t$name\n";
 	}
	#else { print "line not recognized\n$_\n"; }
}
close INPUT;
close OUTPUT;
exit;