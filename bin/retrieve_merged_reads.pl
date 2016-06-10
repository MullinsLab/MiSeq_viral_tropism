#!/usr/bin/perl

##########################################################################################
# program: template_primer_fasta.pl
# From before and after trimmed files to locate trimmed sequences, output as fasta file
# usage: perl template_primer_fasta.pl before_trim_fastq after_trim_fastq
# Author: Wenjie Deng
# Date: 2016-05-04
##########################################################################################


use strict;

my $usage = "usage: perl retrieve_merged_reads.pl pear_merged_fastq min_merge_length_to_retrieve outFastq\n";
my $infile = shift or die $usage;
my $minlen = shift or die $usage;
my $outfile = shift or die $usage;
my $count = my $passcutoffcount = 0;
open IN, $infile or die "couldn't open $infile: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
while (my $name = <IN>) {
	chomp $name;
	my $seq = <IN>;
	chomp $seq;
	my $plus = <IN>;
	my $qual = <IN>;
	++$count;
	if (length $seq >= $minlen) {
		print OUT "$name\n$seq\n$plus$qual";
		++$passcutoffcount;
	}
}
close IN;
close OUT;

print "\n* retrieve_merged_reads.pl: total $count merged reads, $passcutoffcount merged reads >= $minlen bp\n\n";