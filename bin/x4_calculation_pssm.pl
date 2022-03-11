#!/usr/bin/perl

use strict;
use warnings;
use v5.10;


# parse webpssm txt file to calculate the percentage of X4 usage

my $usage = "x4_calculation_pssm.pl inFile outFile\n";
my $inFile = shift or die $usage;
my $outFile = shift or die $usage;
my $ipIdx = shift or die $usage;
my $file = shift || "";

my $total_downloads = my $uniq_users = 0;
my %monthIpStatus = my %ipStatus = ();
open IN, $inFile;
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my @fields = split /\s+/, $line;
	if ($file) {
		if ($line =~ /viroblast\-(.*?)\.tar.gz/) {
			++$total_downloads;
			my $ip = $fields[$ipIdx-1];
			if (!$ipStatus{$ip}) {
				$ipStatus{$ip} = 1;
				++$uniq_users;
			}
		}
	}else {
		++$total_downloads;
		my $ip = $fields[$ipIdx-1];
		if (!$ipStatus{$ip}) {
			$ipStatus{$ip} = 1;
			++$uniq_users;
		}
	}
}
close IN;

#my $meanUsg = int($total_downloads / 12);
#my $meanUniqUsg = int($uniq_users / 12);
open OUT, ">", $outFile or die "couldn't open $outFile: $!\n";
print OUT "total $uniq_users users, $total_downloads downloads\n";
print "total $uniq_users users, $total_downloads downloads\n";