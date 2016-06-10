#!/usr/bin/perl -w

# calucates the consensus from an aligned sequence fasta file and a fastq file
# if gap is majority (> 50%) at one position, the consensus will be gap
# else calculate sum of quality score for each nucleotide at a position, the nucleotide
# having max sum of quality is the consensus at the position

use strict;

my $usage = "perl cal_consensus_qual_from_dir.pl indir outfile\n";

my $indir = shift or die $usage;
my $outfile = shift or die $usage;
my $filecount = 0;
open OUT, ">$outfile" or die "couldn't open $outfile: $!\n";
opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $name = readdir DIR) {
	if ($name =~ /(.*?)\.aln/) {
		my $inFASTA = $indir."/$name";
		my $inFASTQ = $indir."/$1.fastq";
		my $consname = $1."_consensus";		
		my $count = my $alignlen = 0;
		my (@seqNameNseq, %fastqnameSeq, %fastqnameQual, %fastanameSeq, %fastanameQual, $seqName, @seqNames);

		open FASTQ, $inFASTQ or die "couldn't open $inFASTQ: $!\n";
		while (my $name = <FASTQ>) {
			chomp $name;
			my $seq = <FASTQ>;
			chomp $seq;
			my $plus = <FASTQ>;
			my $qual = <FASTQ>;
			chomp $qual;
			$fastqnameSeq{$name} = $seq;
			$fastqnameQual{$name} = $qual;
		}
		close FASTQ;

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
		my %posCoverage = my %posGapcoverage = my %posNtQual = ();
		foreach my $name (@seqNames) {
			my $alignseq = my $seq = $fastanameSeq{$name};
			$seq =~ s/\-//g;
			if (length $seq > length $fastqnameSeq{$name}) {
				die "impossible: $name\n";
			}
			if ($seq ne $fastqnameSeq{$name}) {
				my $idx = index($fastqnameSeq{$name}, $seq);
				if ($idx >= 0) { 
					$fastqnameQual{$name} = substr($fastqnameQual{$name}, $idx, length $seq);
				}else {
					die "something wrong\n";
				}
			}
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
			my @ntquals = split //, $fastqnameQual{$name};
			for (my $i = $start; $i <= $end; $i++) {
				++$posCoverage{$i};
				my $nt = $alignnts[$i];
				if ($nt eq '-') {
					++$posGapcoverage{$i};
				}else {
					my $qual = shift @ntquals;
					$qual = ord($qual) - 33;
					$posNtQual{$i}{$nt} += $qual;
				}
			}
			die "something wrong in @ntquals\n" if (@ntquals);
		}
		my %consNt = ();
		for (my $i = 0; $i < $alignlen; $i++) {
			if ($posCoverage{$i} >= 3) {	# requires consensus of at least 3 coverage;
				if ($posGapcoverage{$i} and $posGapcoverage{$i} / $posCoverage{$i} > 0.5) {
					$consNt{$i} = '-';
				}else {
					my $maxqual = 0;
					my $cons = '';
					foreach my $nt (sort{$a cmp $b} keys %{$posNtQual{$i}}) {
#						if ($i == 274) {
#							print "nt: $nt, posNtQual{$i}{$nt}: $posNtQual{$i}{$nt}\n";
#						}
				
						if ($posNtQual{$i}{$nt} > $maxqual) {
							$cons = $nt;
							$maxqual = $posNtQual{$i}{$nt};
						}elsif ($posNtQual{$i}{$nt} == $maxqual) {
							$cons .= $nt;
						}
					}
#					if ($i == 274) {
#						print "Cons: $cons\n";
#					}
					if ($cons) {
						if (length $cons == 1) {
							$consNt{$i} = $cons;
						}else {
							$consNt{$i} = Ambiguitiescons($cons);
							die "no ambiguities consensus for $cons\n" if (!$consNt{$i});
						}
					}else {
						die "cannot get cons\n";
					}			
				}
			}
		}

		my $consseq = '';
		for (my $i = 0; $i < $alignlen; $i++) {
			if (!$consNt{$i}) {
				$consNt{$i} = '';
			}
			$consseq .= $consNt{$i};
		}
		print OUT ">$consname\n";
		print OUT $consseq,"\n";
		++$filecount;
	}
}
close OUT;

print "\n* cal_consensus_qual_from_dir.pl: implemented $filecount alignment files\n";


sub Ambiguitiescons {
	my $nts = shift;
	my %nts2Ambig = (
		"AG"   => "R",
		"CT"   => "Y",
		"GT"   => "K",
		"AC"   => "M",
		"CG"   => "S",
		"AT"   => "W",
		"CGT"  => "B",
		"AGT"  => "D",
		"ACT"  => "H",
		"ACG"  => "V",
		"ACGT" => "N",		
	);
	my $ambig = $nts2Ambig{$nts};
	return $ambig;
}
