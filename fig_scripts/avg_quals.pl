#!/usr/bin/perl -w

### This program will take as input a directory of fastq and outputs the 
### average quality of each sequence into a renamed file

use strict;
use warnings;
use POSIX;

### change record separator
#$/ = "^@";

### from this directory
my $fastq_dir = "../fig_data/2018_06_14_albacore/2017_01_24_albacore_fastq";

opendir(DIR, $fastq_dir) || die "Can't open directory\n";
my @fastq = grep(/\.fastq$/, readdir(DIR));

open NUM_OUTPUT, ">../fig_data/read_numbers.txt" || die "Can't open numbers file\n";
print NUM_OUTPUT "file\ttotal_reads\ttotal_bp\tmean_length\n";


foreach my $infile (glob("$fastq_dir/*.fastq")) {
	open INPUT, $infile || die "Can't open fastq $infile\n";
	my $outfile = $infile;
	$outfile =~ s/\.fastq/_qual\.txt/;
	open OUTPUT, ">$outfile" || die "Can't open fastq $outfile\n";
	print "$infile\t$outfile\n";
	my $i = 0;
	my $total_bp = 0;
	my @slurp = <INPUT>;
	#while(<INPUT>) {
	### this is stupid, but simple.
	while (my @fastq_seq = splice(@slurp, 0, 4)) {
		#if(my ($name, $seq, $qual) = /(channel.+fast5)\n([ACGT]+)\n\+\n(\S+)\n/) {
		$i++;
		chomp(my $fname = $fastq_seq[0]);
		chomp(my $fseq = $fastq_seq[1]);
		chomp(my $fqual = $fastq_seq[3]);
		#print $fname, "\n", $fastq_seq[1], "\n", $fastq_seq[2], "\n", $fqual, "\n\n";
		my @qual = split('', $fqual);			
		my $numq=0;
		my @sortq=();		
		for(@qual) {
			$numq+= (1 - 10**(-((ord($_) - 33)/10)));
			push @sortq, ord($_);
		}
		$numq/=length($fseq);
		$total_bp+=length($fseq);
		@sortq = sort {$a <=> $b} @sortq;
		my $medq = (1 - 10**(-(($sortq[floor(length($fseq)/2)] - 33)/10)));
		print OUTPUT "$fname\t", length($fseq), "\t$numq\t$medq\n";
		#exit;

		#}
		#else { print "line not recognized\n$_\n"; }
	}
	my $bcname = $infile;
	$bcname =~ s/\.fastq//;
	$bcname =~ s/.+barcode/NB/;
	print NUM_OUTPUT "$bcname\t$i\t", $total_bp/1000000,"\t",$total_bp/$i,"\n";
	close INPUT;
	close OUTPUT;
}
exit;