#!/usr/bin/perl -w

use strict;
use File::Path;

my $usage = "perl extract_v3.pl invariantsfile inampreffile outv3file\n";
my $infile = shift or die $usage;
my $reffile = shift or die $usage;
my $outfile = shift or die $usage;
my $ref = "";
open REF, $reffile or die "couldn't open $reffile: $!\n";
while (my $line = <REF>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	unless ($line =~ /^>/) {
		$ref .= $line;
	}
}
close REF;

my $tmpfile = "tmp.fas";
my $alnvariantfile = $infile;
$alnvariantfile =~ s/\.fasta/\.aln/;
open TMP, ">", $tmpfile or die "couldn't open $tmpfile: $!\n";
print TMP ">reference\n$ref\n";
open IN, $infile or die "couldn't open $infile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	print TMP "$line\n";
}
close IN;
close TMP;

system("muscle -in $tmpfile -out $alnvariantfile -quiet");

unlink($tmpfile);
my $name = "";
my $varcount = my $v3varcount = 0;
my (%nameSeq, @names, %nameCount, %v3Count);
open ALN, $alnvariantfile or die "couldn't open $alnvariantfile: $!\n";
while (my $line = <ALN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		unless ($name eq "reference") {
			push @names, $name;
			if ($name =~ /_(\d+)$/) {
				$nameCount{$name} = $1;
			}else {
				die "no duplicate info\n";
			}
			++$varcount;
		}
	}else {
		$nameSeq{$name} .= $line;
	}
}
close ALN;
my $alnlen = my $idx = my $v3start = my $v3end = my $coverv3 = my $notcoverv3 = 0;
if ($nameSeq{"reference"}) {
	$alnlen = length $nameSeq{"reference"};
	my @refnts = split //, $nameSeq{"reference"};
	for (my $i = 0; $i < scalar @refnts; $i++) {
		if ($refnts[$i] =~ /[ACGT]/) {
			++$idx;
		}
		if ($idx == 132) {	# there are v3 index in reference sequence
			$v3start = $i;
		}elsif ($idx == 239) {
			$v3end = $i;
			last;
		}
	}
}

foreach my $name (@names) {
	my $start = my $end = 0;
	my @nts = split //, $nameSeq{$name};
	if (length $nameSeq{$name} != $alnlen) {
		die "sequences not aligned\n";
	}
	for (my $i = 0; $i < $alnlen; $i++) {
		if ($nts[$i] =~ /[A-Z]/) {
			$start = $i;
			last;
		}
	}
	for (my $i = $alnlen - 1; $i >= 0 ; $i--) {
		if ($nts[$i] =~ /[A-Z]/) {
			$end = $i;
			last;
		}
	}
	if ($start <= $v3start and $end >= $v3end) {
		my $v3 = substr($nameSeq{$name}, $v3start, $v3end - $v3start + 1);
		$v3 =~ s/\-//g;
		$v3Count{$v3} += $nameCount{$name};
		++$coverv3;
	}else {
		++$notcoverv3;
	}
}

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
$idx = 0;
foreach my $v3 (sort{$v3Count{$b} <=> $v3Count{$a}} keys %v3Count) {
	++$idx;
	my $name = "v3_$idx"."_$v3Count{$v3}";
	print OUT ">$name\n$v3\n";
}
close OUT;

print "\n* extract_v3.pl: implemented $varcount variants, $coverv3 cover v3 loop, $notcoverv3 not, extract $idx v3 variants\n";

