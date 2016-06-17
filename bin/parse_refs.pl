#!/usr/bin/perl

use strict;

my $fpfile = shift;
my $rpfile = shift;
my $rtfile = shift;
my $fp = my $rp = my $rt = "";
open FP, $fpfile;
while (my $line = <FP>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	unless ($line =~ /^>/) {
		$fp = $line;
		last;
	}
}
close FP;

open RP, $rpfile;
while (my $line = <RP>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	unless ($line =~ /^>/) {
		$rp = $line;
		last;
	}
}
close RP;

open RT, $rtfile;
while (my $line = <RT>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	unless ($line =~ /^>/) {
		$rt = $line;
		$rt =~ tr/ACGT/TGCA/;
		$rt = reverse $rt;
		last;
	}
}
close RT;

my $len = length $fp;
if ($len > length $rp) {
	$len = length $rp;
}

my $minoverlap = int($len * 0.9);
my $RTminoverlap = int(0.9 * length $rt);

print "$rt\n$fp\n$rp\n$minoverlap\n$RTminoverlap\n";
