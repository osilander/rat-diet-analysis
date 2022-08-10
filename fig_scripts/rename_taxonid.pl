#!/usr/bin/perl -w

### This program will take as input a list of read names
### with associated orders, families, etc., as well as a list as Krona-
### compatible taxa. Output is a new Krona-compatible MEGAN file.
### This can be turned into a Krona plot by typing:
### ktImportTaxonomy name_of_krona_megan_file.txt at the command
### prompt

use strict;
use warnings;

### from this read name / taxon name file
my $taxon_name_file = "../fig_data/taxa_level_changes_fp.txt";
my $outfile = $taxon_name_file;
$outfile =~ s/\.txt/_taxonid\.txt/;

my $taxa_file = "../fig_data/taxonomy.tab";

my %taxa_id;

open INPUT, $taxa_file || die "Can't open taxa file\n";
while (<INPUT>) {
	chomp;
	if(my ($node, $name) = /(\d+)\t\d+\t\d+\t.+\t(.+)/) {
		$taxa_id{$name} = $node;
	}
	else { print "line not recognized\n$_\n"; }
}
close INPUT;

open INPUT, $taxon_name_file || die "Can't open name file\n";
open OUTPUT, ">$outfile" || die "Can't open output file\n";
while (<INPUT>) {
	if(my ($read, $taxa) = /\d+\t(\S+)\t(\w+)/) {
		if(exists($taxa_id{$taxa})) { print OUTPUT "$read\t$taxa_id{$taxa}\n"; }
	}
	else { print "taxon not recognized\n$_\n"; }
}
close INPUT;
close OUTPUT;
exit;