#!/usr/bin/perl -w

### This program will take as input slightly modified MEGAN file of the leaves of a 
### certain taxon depth (e.g. "Class"), and create a table in which the first column
### is the name of the taxon (e.g. "Ascomycota"), and the second is the name of the read
### this can then be parsed easily in R

use strict;
use warnings;

my $infile= "../fig_data/MEGAN_phylum.txt";

open INPUT, $infile || die "Can't open taxa $infile\n";
my $outfile = $infile;
$outfile =~ s/MEGAN/MEGAN_parse/;

open OUTPUT, ">$outfile" || die "Can't open $outfile\n";
while(<INPUT>) {
	if(my ($taxon, $reads) = /(\w+)\t(channel.+)/) {
		my @reads = split('\t', $reads);			
		for(@reads) {
			print OUTPUT "$_\t$taxon\n";
		}
	}
	else { print "line not recognized\n$_\n"; }
}
close INPUT;
close OUTPUT;
exit;