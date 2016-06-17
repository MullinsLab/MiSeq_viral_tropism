#!/usr/bin/perl -w

use strict;

my $usage = "perl template_family_size_distribution.pl inR2TemplateFile outFile\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $min = shift;
my $max = shift;
my $totalcount = 0;
my %countCount = ();
open IN, $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if ($line =~ /^\s*$/ || $line =~ /^templateid/);
	my ($id, $count) = split /\t/, $line;
	if (defined $min and defined $max) {
		if (length $id >= $min and length $id <= $max) {
			++$countCount{$count};
			++$totalcount;
		}
	}else {		
		++$countCount{$count};
		++$totalcount;
	}	
}
close IN;

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
foreach my $count (sort{$a <=> $b} keys %countCount) {
	my $percentage = $countCount{$count} / $totalcount;
	print OUT "$count\t$countCount{$count}\t$percentage\n";
}
close OUT;

print "\n* template_family_size_distribution.pl: total $totalcount template ids with length between $min and $max bps\n";
