#!/usr/bin/perl -w

#########################################################################################################
# Program: sampleSummary.pl
# Purpose: go thorough each sample directory, read sample's .log file, output a summary file
# for all samples
# Output: .csv file
# Author: Wenjie Deng
# Date: 2019-06-26
###########################################################################################################

use strict;

my $usage = "perl sampleSummary.pl outSampleSummaryfile directory_where_samples_results\n";
my $outfile = shift or die $usage;
my $indir = shift || '.'; # default current directory
my $samplecount = 0;

opendir my $dh, $indir or die "couldn't open $indir: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
print OUT "Sample,Raw reads,After QC,Merged,Merged (>=300bp),Processed,Max UMI family size,fs cutoff,UMIs,variants,V3 variants\n";
while (my $name = readdir $dh) {
	if (-d $name) {
		opendir my $sub, $name or die "couldn't open $name: $!\n";
		while (my $file = readdir $sub) {
			if ($file =~ /$name.log/) {
				++$samplecount;
				my @values = ();
				push @values, $name;
				$file = $name.'/'.$file;
				open IN, $file or die "couldn't open $file: $!\n";
				while (my $line = <IN>) {
					chomp $line;
					next if $line =~ /^\s*$/;
					if ($line =~ /Total input FastQ records: (\d+)/) {
						push @values, $1;
					}elsif ($line =~ /FastQ paired records kept: (\d+)/) {
						push @values, $1;
					}elsif ($line =~ /retrieve_merged_reads.pl: total (\d+) merged reads, (\d+) merged reads >= 300 bp/) {
						push @values, $1, $2;
					}elsif ($line =~ /template_primer_fasta.pl: total \d+ reads after, write (\d+) to fasta/) {
						push @values, $1;
					}elsif ($line =~ /\d+ templates' family size >= (\d+)/) {
						push @values, $1;
					}elsif ($line =~ /(\d+) qualified for supermajority cutoff of/) {
						push @values, $1;
					}elsif ($line =~ /unique_consensus.pl: implemented \d+ consensus sequences, total (\d+) variants/) {
						push @values, $1;
					}elsif ($line =~ /(\d+) v3 rariants in the range of 1.2/) {
						push @values, $1;
					}
				}
				close IN;
				print OUT join(',', @values), "\n";
			}
		}
		closedir $sub;
	}
}
closedir $dh;
close OUT;

print "Total $samplecount samples\n";

