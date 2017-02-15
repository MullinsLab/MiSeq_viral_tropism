#!/usr/bin/perl -w

# calucates the consensus from an aligned sequence fasta file, requires the frequency of 
# supermojority nucleotide (include gap) at each position greater than a cutoff user 
# privides. If there is one position where the supermojority less than cutoff, the
# consensus will not be calculated and discard. 

use strict;

my $usage = "perl cal_consensus_supermajority_from_dir.pl indir templateLogFile outfile frequency_cutoff(default 0)\n";

my $indir = shift or die $usage;
my $logfile = shift or die $usage;
my $outfile = shift or die $usage;
my $cutoff = shift || 0;	# default 0 to calculate consensus with simple majority
my $tfamilycut = my $filecount = my $qualified = my $unqualified = 0;
open LOG, $logfile or die "couldn't open $logfile: $!\n";
while (my $line = <LOG>) {
	chomp $line;
	if ($line =~ /family size >= (\d+)/) {
		$tfamilycut = $1;
	}
}
close LOG;
die "No template family size cutoff\n" if (!$tfamilycut);

open OUT, ">$outfile" or die "couldn't open $outfile: $!\n";
opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $name = readdir DIR) {
	if ($name =~ /(.*?)\.aln/) {
		my $inFASTA = $indir."/$name";
		my $consname = $1."_consensus";		
		my $count = my $alignlen = 0;
		my (%fastanameSeq, $seqName, @seqNames);
		open FASTA, $inFASTA or die "couldn't open $inFASTA: $!\n";
		while (my $line = <FASTA>) {
			chomp $line;
			next if ($line =~ /^\s*$/);
			if ($line =~ /^>(.*)$/) {
				if ($count == 1) {
					$alignlen = length $fastanameSeq{$seqName};
				}elsif ($count > 1) {
					if (length $fastanameSeq{$seqName} != $alignlen) {
						die "sequences not aligned\n";
					}
				}
				$seqName = $1;
				push @seqNames, $seqName;
				$count++;
			}else {
				if (!$fastanameSeq{$seqName}) {
					$fastanameSeq{$seqName} = '';
				}
				$fastanameSeq{$seqName} .= $line;		
			}
		}
		# last sequence
		if (length $fastanameSeq{$seqName} != $alignlen) {
			die "sequences not aligned\n";
		}
		close FASTA;
		my %posCoverage = my %posNtcoverage = ();
		foreach my $name (@seqNames) {
			my $alignseq = $fastanameSeq{$name};
			my $start = my $end = 0;
			my @alignnts = split //, $alignseq;
			for (my $i = 0; $i < $alignlen; $i++) {
				if ($alignnts[$i] =~ /[ACGT]/) {
					$start = $i;
					last;
				}
			}
			for (my $i = $alignlen-1; $i >= 0; $i--) {
				if ($alignnts[$i] =~ /[ACGT]/) {
					$end = $i;
					last;
				}
			}
			for (my $i = $start; $i <= $end; $i++) {
				++$posCoverage{$i};
				my $nt = $alignnts[$i];
				++$posNtcoverage{$i}{$nt};
			}
		}
		my %consNt = ();
		my $flag = 1;
		for (my $i = 0; $i < $alignlen; $i++) {
			if ($posCoverage{$i} >= $tfamilycut) {
				my $cons = my $firstnt = "";
				my $count = my $firstfreq = 0;
				foreach my $nt (sort{$posNtcoverage{$i}{$b} <=> $posNtcoverage{$i}{$a}} keys %{$posNtcoverage{$i}}) {
					++$count;
					my $freq = $posNtcoverage{$i}{$nt} / $posCoverage{$i};
					if ($cutoff == 0) {	# simple majority
						if ($count == 1) {
							if ($freq > 0.5) {
								$cons = $nt;
								last;
							}else {
								$firstfreq = $freq;
								$firstnt = $nt;
							}						
						}else {
							if ($firstfreq > $freq) {
								$cons = $firstnt;
							}
							last;
						}
					}else {	# super majority
						if ($cutoff == 0.5) {
							if ($freq > $cutoff) {
								$cons = $nt;
							}
						}elsif ($cutoff > 0.5) {
							if ($freq >= $cutoff) {
								$cons = $nt;
							}
						}else {
							die "cutoff must be >= 0.5\n";
						}
						last;
					}					
				}
				if ($cons) {
					$consNt{$i} = $cons;
				}else {
					%consNt = ();
					$flag = 0;
					++$unqualified;
					last;
				}
			}
		}
		if ($flag) {
			my $consseq = '';
			for (my $i = 0; $i < $alignlen; $i++) {
				if (!$consNt{$i}) {
					$consNt{$i} = '';
				}
				$consseq .= $consNt{$i};
			}
			print OUT ">$consname\n";
			print OUT $consseq,"\n";
			++$qualified;
		}		
		++$filecount;
	}
}
close OUT;

print "\n* cal_consensus_supermajority_from_dir.pl: implemented $filecount alignment files, $qualified qualified for supermajority cutoff of $cutoff, template family size cutoff of $tfamilycut, $unqualified not\n";

