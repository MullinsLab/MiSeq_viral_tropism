#!/usr/bin/perl -w

use strict;

my $usage = "perl template_size_distribuion.pl inR2TemplateFile outFile\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $totalcount = 0;
my %sizeCount = ();
open IN, $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if ($line =~ /^\s*$/ || $line =~ /^templateid/);
	my ($id, $count) = split /\t/, $line;
	my $idlength = length $id;
	++$sizeCount{$idlength};
	++$totalcount;
}
close IN;

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
foreach my $len (sort{$a <=> $b} keys %sizeCount) {
	my $percentage = $sizeCount{$len} / $totalcount;
	print OUT "$len\t$sizeCount{$len}\t$percentage\n";
}
close OUT;

print "\n* template_size_distribution.pl: total $totalcount template ids\n";
