#!/usr/bin/perl

use strict;

my $usage = "perl unique_consensus.pl inconsfile outvariantfile\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my %seqCount = ();
my $count = my $idx = 0;
open IN, $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	unless ($line =~ /^>/) {
		$line =~ s/\-//g;
		++$seqCount{$line};
		++$count;
	}
}
close IN;

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
foreach my $seq (sort{$seqCount{$b} <=> $seqCount{$a}} keys %seqCount) {
	++$idx;
	my $name = "variant_$idx"."_$seqCount{$seq}";
	print OUT ">$name\n$seq\n";
}
close OUT;

print "\n* unique_consensus.pl: implemented $count consensus sequences, total $idx variants\n";
