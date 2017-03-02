#!/usr/bin/perl -w

use strict;
use File::Path;
use File::Copy;
use File::Basename;

my $usage = "perl retrieve_sequences_per_template_family_cutoff.pl inTemplateFile inRTtrimmedFastq templateidlength errorrate(default 0.005) templatefamilysizecutoff(default 0)\n";
my $intemplate = shift or die $usage;
my $inRTtrimmed = shift or die $usage;
my $tempidlen = shift or die $usage;
my $errorrate = shift || 0.005;	# overall error rate used to calculate template family size cutoff
my $familycut = shift || 0; # if 0, using cutoff model to determine template family size cutoff
my $dir = dirname($intemplate);
my $totalcount = my $gtcutoff = my $tempidlencount = 0;
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

my $flag = 0;
open IN, $intemplate or die "couldn't open $intemplate: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if ($line =~ /^\s*$/ || $line =~ /^templateid/);
	my ($id, $count, $name) = split /\t/, $line;
	my @names = ();	
	if (length $id == $tempidlen) {
		if (!$flag and $familycut == 0) {
			$flag = 1;
			if ($errorrate == 0.02) {
				if ($count <= 10) {
					$familycut = 2;
				}elsif ($count <= 8500) {
					$familycut = -1.24*10**-21*$count**6 + 3.53*10**-17*$count**5 - 3.90*10**-13*$count**4 + 2.12*10**-9*$count**3 - 6.06*10**-6*$count**2 + 1.80*10**-2*$count + 3.15;
				}else {
					$familycut = 0.0079 * $count + 9.4869;
				}				
			}elsif ($errorrate == 0.01) {
				if ($count <= 10) {
					$familycut = 2;
				}else {
					$familycut = 1.09*10**-26*$count**6 + 7.82*10**-22*$count**5 - 1.93*10**-16*$count**4 + 1.01*10**-11*$count**3 - 2.31*10**-7*$count**2 + 0.00645*$count + 2.872;
				}
			}else { # default for error rate of 0.005
				if ($count <= 10) {
					$familycut = 2;
				}else {
					$familycut = -9.59*10**-27*$count**6 + 3.27*10**-21*$count**5 - 3.05*10**-16*$count**4 + 1.2*10**-11*$count**3 - 2.19*10**-7*$count**2 + 0.004044*$count + 2.273;
				}
			}
			$familycut = int($familycut + 0.5);
			if ($familycut < 3) {
				$familycut = 3;
			}
		}
		if ($count >= $familycut) {
			@names = split /,/, $name;
			my %seqCount = my %seqNames = ();
			my $uniqfasta = $outdir."/".$id."_$count"."_uniq.fasta";
			my $namefile = $outdir."/".$id."_$count"."_uniq_name.txt";
			foreach my $n (@names) {
				my $seq = $nameSeq{$n};
				if ($seq) {		
					++$seqCount{$seq};
					push @{$seqNames{$seq}}, $n;
				}else {
					die "something wrong\n";
				}
			}
			my $idx = 0;
			my %uniqidxNames = ();
			open UFASTA, ">", $uniqfasta or die "couldn't open $uniqfasta: $!\n";
			open NFILE, ">", $namefile or die "couldn't open $namefile: $!\n";
			foreach my $seq (sort{$seqCount{$b} <=> $seqCount{$a}} keys %seqCount) {
				++$idx;
				my $name = "uniq_$idx"."_$seqCount{$seq}";
				print UFASTA ">$name\n$seq\n";
				print NFILE "$name\t", join(',', @{$seqNames{$seq}}), "\n";
				push @{$uniqidxNames{$name}}, @{$seqNames{$seq}};
			}			
			close UFASTA;
			close NFILE;
			# align unique sequences
			my $uniqalignfile = $uniqfasta;
			$uniqalignfile =~ s/\.fasta/\.aln/;
			if (keys %uniqidxNames == 1) {	# one unique sequence, don't need to align
				copy($uniqfasta, $uniqalignfile) or die "copy failed: $!\n";
			}else {
				system("bin/muscle -in $uniqfasta -out $uniqalignfile -quiet");
			}			
			++$gtcutoff;
		}
		++$tempidlencount;			
	}
	++$totalcount;
}
close IN;

my $logfile = "$dir/template_sequences.log";
open LOG, ">", $logfile or die "couldn't open $logfile: $!\n";
print LOG "retrieve_sequences_per_template_family_cutoff.pl: Total $totalcount template ids, $tempidlencount with template id length of $tempidlen, $gtcutoff templates' family size >= $familycut\n";
close LOG;

print "\n* retrieve_sequences_per_template_family_cutoff.pl: Total $totalcount template ids, $tempidlencount with template id length of $tempidlen, $gtcutoff templates' family size >= $familycut\n";
