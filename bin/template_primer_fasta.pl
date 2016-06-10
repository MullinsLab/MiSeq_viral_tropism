#!/usr/bin/perl

##########################################################################################
# program: template_primer_fasta.pl
# From before and after trimmed files to locate trimmed sequences, output as fasta file
# usage: perl template_primer_fasta.pl before_trim_fastq after_trim_fastq
# Author: Wenjie Deng
# Date: 2016-05-04
##########################################################################################


use strict;

my $usage = "usage: perl template_primer_fasta.pl before_trim_fastq after_trim_fastq\n";
my $beforefile = shift or die $usage;
my $afterfile = shift or die $usage;
my $outfile = $afterfile;
$outfile =~ s/\.fastq$/_tnp\.fasta/;
my $count = my $fastacount = 0;
my %afternameSeq = ();
open AFTER, $afterfile or die "couldn't open $afterfile: $!\n";
while (my $name = <AFTER>) {
	chomp $name;
	my $seq = <AFTER>;
	chomp $seq;
	$afternameSeq{$name} = $seq;
	my $plus = <AFTER>;
	my $qual = <AFTER>;
	++$count;
}
close AFTER;

open BEFORE, $beforefile or die "couldn't open $beforefile: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
while (my $name = <BEFORE>) {
	chomp $name;
	my $seq = <BEFORE>;
	chomp $seq;
	my $plus = <BEFORE>;
	my $qual = <BEFORE>;
	if ($afternameSeq{$name}) {
		my $len = length $afternameSeq{$name};
		my $tnpseq = substr($seq, $len);		
		print OUT ">$name\n$tnpseq\n";
		++$fastacount;
	}
}
close BEFORE;
close OUT;
print "\n* template_primer_fasta.pl: total $count reads after, write $fastacount to fasta\n\n";