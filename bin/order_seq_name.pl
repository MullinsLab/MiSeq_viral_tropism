#!/usr/bin/perl -w

# from a sequence fasta file, ordered the sequences by sequence names

use strict;

my $usage = "perl order_seq_name.pl inFasta outFasta\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my $name = "";
my $count = 0;
my (@names, %nameSeq);
open IN, $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		push @names, $name;
		++$count;
	}else {
		$line =~ s/\-//g;
		$nameSeq{$name} .= $line;
	}
}
close IN;
open OUT, ">$outfile" or die "couldn't open $outfile: $!\n";
foreach my $name (sort{$a cmp $b} @names) {
	print OUT ">$name\n$nameSeq{$name}\n";
}
close OUT;

print "total $count sequences\n";

