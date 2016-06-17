#!/usr/bin/perl -w

use strict;
use File::Path;
use File::Basename;

my $usage = "perl retrieve_sequences_per_template_family_gt2.pl inTemplateFile inRTtrimmedFastq templateidlength\n";
my $intemplate = shift or die $usage;
my $inRTtrimmed = shift or die $usage;
my $tempidlen = shift or die $usage;
my $dir = dirname($intemplate);
my $totalcount = my $gt2count = my $passlencut = my $notpasslencut = 0;
my (%nameFullname, %nameSeq, %nameQual);
open FASTQ, $inRTtrimmed or die "couldn't open $inRTtrimmed: $!\n";
while (my $line = <FASTQ>) {
	chomp $line;
	$line =~ /^(\S+)/;
	my $name = $1;
	$nameFullname{$name} = $line;
	my $seq = <FASTQ>;
	chomp $seq;
	my $plus = <FASTQ>;
	my $qual = <FASTQ>;
	chomp $qual;
	$nameSeq{$name} = $seq;
	$nameQual{$name} = $qual;
}
close FASTQ;

my $outdir = "$dir/template_sequences";
rmtree($outdir) if (-e $outdir);
mkdir $outdir;

open IN, $intemplate or die "couldn't open $intemplate: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if ($line =~ /^\s*$/ || $line =~ /^templateid/);
	my ($id, $count, $name) = split /\t/, $line;
	my @names = ();	
	if (length $id == $tempidlen) {
		if ($count >= 3) {
			my $totallen = my $avglen = 0;
			@names = split /,/, $name;
			foreach my $n (@names) {
				$totallen += length $nameSeq{$n};
			}
			$avglen = $totallen / $count;
			if ($avglen > 250) {
				my $fastq = $outdir."/".$id."_$count.fastq";
				my $fasta = $outdir."/".$id."_$count.fasta";
				open FASTQ, ">", $fastq or die "couldn't open $fastq: $!\n";
				open FASTA, ">", $fasta or die "couldn't open $fasta: $!\n";
				foreach my $n (@names) {
					if ($nameSeq{$n}) {
						print FASTQ "$nameFullname{$n}\n$nameSeq{$n}\n+\n$nameQual{$n}\n";
						print FASTA ">$nameFullname{$n}\n$nameSeq{$n}\n";
					}else {
						die "something wrong\n";
					}
				}
				close FASTQ;
				my $alignfile = $outdir."/".$id."_$count.aln";
				system("bin/muscle -in $fasta -out $alignfile -quiet");
				++$passlencut;
			}else {
				++$notpasslencut;
			}			
			++$gt2count;
		}
		++$totalcount;	
	}
}
close IN;

my $logfile = "$dir/template_sequences.log";
open LOG, ">", $logfile or die "couldn't open $logfile: $!\n";
print LOG "retrieve_sequences_per_template_family_gt2.pl: Total $totalcount template ids, $gt2count templates' family size >= 3, $passlencut templates average fragment size > 250bp, $notpasslencut not\n";
close LOG;

print "\n* retrieve_sequences_per_template_family_gt2.pl: Total $totalcount template ids, $gt2count templates' family size >= 3, $passlencut templates average fragment size > 250bp, $notpasslencut not\n";
