#!/usr/bin/perl -w

use strict;

my $usage = "perl fragment_size_distribuion.pl inRTprimerTrimmedFastq outFile\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $totalcount = my $max = 0;
my %sizeCount = ();
open IN, $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	my $seq = <IN>;
	chomp $seq;
	my $plus = <IN>;
	my $qual = <IN>;
	my $len = length $seq;
	++$sizeCount{$len};
	++$totalcount;
}
close IN;

foreach my $len (sort{$a <=> $b} keys %sizeCount) {
	if ($len > $max) {
		$max = $len;
	}
}

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
foreach (0 .. $max) {
	if (!$sizeCount{$_}) {
		$sizeCount{$_} = 0;
	}
	my $percentage = $sizeCount{$_} / $totalcount;
	print OUT "$_\t$sizeCount{$_}\t$percentage\n";
}
close OUT;

print "\n* fragment_size_distribution.pl: total $totalcount merged reads\n";
