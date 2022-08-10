#!/usr/bin/perl -w

### This program will take as input a list of read names
### with associated orders, families, etc., as well as a list
### of read names and orders, families, etc. from a MEGAN file.
### All read names that match the input read names will be
### reassigned to the 
### order, family, etc. from the input list, but as Krona-
### compatible taxa. Output is a new Krona-compatible MEGAN file.
### This can be turned into a Krona plot by typing:
### ktImportTaxonomy name_of_krona_megan_file.txt at the command
### prompt

use strict;
use warnings;

### from this krona file
my $krona_file = "../fig_data/FinalDataSetAll_krona.txt";

### to make new krona files that exclude this list of taxa
my $exclude_file = "../fig_data/taxa_to_remove.txt";

my $outfile = $krona_file;
$outfile =~ s/\.txt/_removed\.txt/;

my $taxa_file = "../fig_data/taxonomy.tab";
my %taxa_children;
my %exclude_names;
my %exclude_taxa_hash;
my @exclude_nums;
my @all_exclude_taxa;

open INPUT, $exclude_file || die "Can't open exclude file\n";
while (<INPUT>) {
	chomp;
	$exclude_names{$_}=1;
}
close INPUT;

print "\nExcluding taxa:\n";
open INPUT, $taxa_file || die "Can't open taxa file\n";
while (<INPUT>) {
	chomp;
	if(my ($node, $unk, $parent, $level, $name) = /(\d+)\t(\d+)\t(\d+)\t(.+)\t(.+)/) {
		push @{$taxa_children{$parent}}, $node;
		if(exists($exclude_names{$name})) {
			push @exclude_nums, $node;
			print "$name\t$node\n";
		}
	}
	else { print "line not recognized\n$_\n"; }
}
close INPUT;

push @all_exclude_taxa, @exclude_nums;

print "\n\n";

### This is a bit silly, but the easiest way to implement
### tree traversal: collect all the children of a node, then all the children of 
### those children, etc.
for my $parent (@exclude_nums) {
	if(exists($taxa_children{$parent})) {
		#print "parent1: $parent\n";
		my @children = @{$taxa_children{$parent}};
		#print "@children ", "\n\n";
		push @all_exclude_taxa, @children;
		for my $parent (@children) {
			if(exists($taxa_children{$parent})) {
				#print "parent2: $parent\n";
				my @children = @{$taxa_children{$parent}};
				#print "@children ", "\n\n";
				push @all_exclude_taxa, @children;
				for my $parent (@children) {
					if(exists($taxa_children{$parent})) {
						#print "parent3: $parent\n";
						my @children = @{$taxa_children{$parent}};
						#print "@children ", "\n\n";
						push @all_exclude_taxa, @children;
						for my $parent (@children) {
							if(exists($taxa_children{$parent})) {
								#print "parent4: $parent\n";
								my @children = @{$taxa_children{$parent}};
								#print "@children ", "\n\n";
								push @all_exclude_taxa, @children;
								for my $parent (@children) {
									if(exists($taxa_children{$parent})) {
										#print "parent5: $parent\n";
										my @children = @{$taxa_children{$parent}};
										#print "@children ", "\n\n";
										push @all_exclude_taxa, @children;
										for my $parent (@children) {
											if(exists($taxa_children{$parent})) {
												#print "parent6: $parent\n";
												my @children = @{$taxa_children{$parent}};
												#print "@children ", "\n\n";
												push @all_exclude_taxa, @children;
												for my $parent (@children) {
													if(exists($taxa_children{$parent})) {
														#print "parent7: $parent\n";
														my @children = @{$taxa_children{$parent}};
														#print "@children ", "\n\n";
														push @all_exclude_taxa, @children;
														for my $parent (@children) {
															if(exists($taxa_children{$parent})) {
																print "parent8: $parent\n";
																print "warning: still more children\n";
																my @children = @{$taxa_children{$parent}};
																print "@children ", "\n\n";
																push @all_exclude_taxa, @children;
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

%exclude_taxa_hash = map { $_ => 1 } @all_exclude_taxa;

open INPUT, $krona_file || die "Can't open krona file\n";
open OUTPUT, ">$outfile" || die "Can't open output file\n";

while (<INPUT>) {
	if(my ($read, $taxa) = /(\S+)\t(\-?\d+)/) {
		unless(exists($exclude_taxa_hash{$taxa})) { print OUTPUT "$read\t$taxa\n"; }
	}
	else { print "line not recognized\n$_\n"; }
}
close INPUT;
close OUTPUT;
exit;