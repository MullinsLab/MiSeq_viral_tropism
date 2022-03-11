#!/usr/bin/perl -w

##########################################################################################
# Program: merge_files.pl
# Purpose: In a directory of files, merges files based on file name
# Author: Wenjie Deng
# Date: 2021-03-01
##########################################################################################

use strict;
use warnings;
use v5.10;
use Getopt::Long;

my %option = (
	'id' => '.',
);

my $usage = "\nusage: merge_files.pl [-option value]

options:  
-id		input directory with fasta files (default: . )

";

GetOptions (\%option, 'id=s');

my $indir = $option{'id'} or die $usage;
my $outdir = $indir."/outputs";
my @programs = ("g2p", "PSSM");
my (%idProgramStatus, @idprograms, %idprogramnamestatus, %idprogramnames, %idheader, %idprogramlines, %idprogramexclude, %idprogramtotaltemplates, %idprogramtemplates, %idprogramfprless2, %idprogramfprless575, %idprogramx4);
opendir DIR, $indir or die "couldn't open $indir: $!\n";
while (my $file = readdir DIR) {
	if ($file =~ /(.*)\.csv|txt$/) {		
		if ($file =~ /^(.*?)_(.*)$/) {
			my $id = $1;
			my $rest = $2;			
			if ($id =~ /(.*?)\-\d+/) {
				$id = $1;
			}
			my $program = "";	
			if ($rest =~ /g2p/) {
				$program = "g2p";
			}elsif ($rest =~ /PSSM/) {
				$program = "PSSM";
			}else {
				die "no program for $rest\n";
			}
			if (!$idProgramStatus{$id}{$program}) {
				$idProgramStatus{$id}{$program} = 1;
				my $idprogram = $id."_".$program;
				push @idprograms, $idprogram;
			}
			$file = $indir."/".$file;
			open IN, $file or die "couldn't open $file: $!\n";
			while (my $line = <IN>) {
				$line =~ s/\R$//;
				next if ($line =~ /^\s*$/ or $line =~ /^WebPSSM/ or $line =~ /^name/ or $line =~ /^header/);
				if ($line =~ /^v3_/) {
					my $name = "";
					if ($program eq "g2p") {
						my @fields = split /,/, $line;
						$name = $fields[0];
						my $comment = $fields[5];
						if ($comment =~ /stop codons/ or $comment =~ /Problems/ or $comment =~ /no C at/) {
							$idprogramexclude{$id}{$program}{$name} = 1;
						}
					}elsif ($program eq "PSSM") {
						my @fields = split /\t/, $line;
						$name = $fields[0];
						if ($fields[1] =~ /\*/ or $fields[3] eq "") {
							$idprogramexclude{$id}{$program}{$name} = 1;
						}
					}else {
						die "No program: $program\n";
					}
					if (!$idprogramnamestatus{$id}{$program}{$name}) {
						$idprogramnamestatus{$id}{$program}{$name} = 1;
						push @{$idprogramnames{$id}{$program}}, $name;
						push @{$idprogramlines{$id}{$program}}, $line;
					}else {
						print "duplicate name: $id, $program, $name\n";
					}
					
				}else {
					die "unformatted entry: $id, $program, $line\n";
				}
			}
			close IN;
		}
	}
}
closedir DIR;

foreach my $id (sort{$a cmp $b} keys %idprogramnames) {
	my @g2pnames = sort(@{$idprogramnames{$id}{'g2p'}});
	my @pssmnames = sort(@{$idprogramnames{$id}{'PSSM'}});
	my $g2pvariants = scalar @g2pnames;
	my $pssmvariants = scalar @pssmnames;
	if ($g2pvariants != $pssmvariants) {
		die "variant number is different between g2p and pssm: $id, $g2pvariants vs $pssmvariants\n";
	}
	while (@g2pnames) {
		my $g2pname = shift @g2pnames;
		my $pssmname = shift @pssmnames;
		if ($g2pname ne $pssmname) {
			die "name is different: id, $g2pname from g2p vs $pssmname from pssm\n";
		}
	}
}

unless (-e $outdir) {
	mkdir $outdir;
}

my $g2pexcludefile = $outdir."/g2p_exclude.csv";
my $pssmexcludefile = $outdir."/pssm_exclude.txt";
open GEXCL, ">", $g2pexcludefile or die "couldn't open $g2pexcludefile: $!\n";
open PEXCL, ">", $pssmexcludefile or die "couldn't open $pssmexcludefile: $!\n";
print GEXCL "sample,header,V3-LOOP,SUBTYPE,FPR,PERCENTILE,COMMENT\n";
print PEXCL "sample\tname\ttranslated amino acid sequence\tscore\tpred\tx4.pct\tr5.pct\tgeno\tpos.chg\tnet.chg\tpercentile\n";

foreach my $id (sort {$a cmp $b} keys %idprogramlines) {
	foreach my $program (keys $idprogramlines{$id}) {
		foreach my $line (@{$idprogramlines{$id}{$program}}) {
			my $name = "";
			if ($program eq "g2p") {
				($name) = split /,/, $line;
				if ($idprogramexclude{$id}{$program}{$name}) {
					print GEXCL "$id,$line\n";
				}				
			}elsif ($program eq "PSSM") {
				($name) = split /\t/, $line;
				if ($idprogramexclude{$id}{$program}{$name}) {
					print PEXCL "$id\t$line\n";
				}
			}
			my ($v3, $idx, $duplicates) = split /\_/, $name;
			$idprogramtotaltemplates{$id}{$program} += $duplicates;
			if ($idprogramexclude{$id}{"g2p"}{$name} or $idprogramexclude{$id}{"PSSM"}{$name}) {
#				print "$id, $program, $line\n";				
			}else {
				if ($program eq "g2p") {
					my @fields = split /,/, $line;
					my ($v3, $idx, $duplicates) = split /\_/, $fields[0];
					$idprogramtemplates{$id}{$program} += $duplicates;
					my $fpr = $fields[3];
					if ($fpr < 2) {
						$idprogramfprless2{$id}{$program} += $duplicates;
					}
					if ($fpr < 5.75) {
						$idprogramfprless575{$id}{$program} += $duplicates;
					}
				}elsif ($program eq "PSSM") {
					my @fields = split /\t/, $line;
					my ($v3, $idx, $duplicates) = split /\_/, $fields[0];
					$idprogramtemplates{$id}{$program} += $duplicates;
					if ($fields[3] == 1) {
						$idprogramx4{$id}{$program} += $duplicates;
					}
				}
			}
		}		
	}
}
close GEXCL;
close PEXCL;

my $count = 0;
my $outfile = $outdir."/x4_calculation.csv";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
print OUT "sample,total_templates,valid_templates,pssm_x4_count,pssm_x4%,pssm_tropism_call,g2p_FPR2_x4_count,g2p_FPR2_x4%,g2p_FPR2_tropism_call,g2p_FPR575_x4_count,g2p_FPR575_x4%,g2p_FPR575_tropism_call,g2p_tropism_call\n";
foreach my $id (sort {$a cmp $b} keys %idprogramlines) {
	++$count;
	if ($idprogramtotaltemplates{$id}{"g2p"} != $idprogramtotaltemplates{$id}{"PSSM"}) {
		die "total templates are different. $id, $idprogramtotaltemplates{$id}{'g2p'} vs $idprogramtotaltemplates{$id}{'PSSM'}\n";
	}elsif ($idprogramtemplates{$id}{"g2p"} != $idprogramtemplates{$id}{"PSSM"}) {
		die "valid templates are different. $id, $idprogramtemplates{$id}{'g2p'} vs $idprogramtemplates{$id}{'PSSM'}\n";
	}
	if (!$idprogramx4{$id}{'PSSM'}) {
		$idprogramx4{$id}{'PSSM'} = 0;
	}
	if (!$idprogramfprless2{$id}{'g2p'}) {
		$idprogramfprless2{$id}{'g2p'} = 0;
	}
	if (!$idprogramfprless575{$id}{'g2p'}) {
		$idprogramfprless575{$id}{'g2p'} = 0;
	}
	my $pssmx4pct = int ($idprogramx4{$id}{'PSSM'} / $idprogramtemplates{$id}{'PSSM'} * 100000 + 0.5) / 1000;
	my $fpr2x4pct = int ($idprogramfprless2{$id}{'g2p'} / $idprogramtemplates{$id}{'g2p'} * 100000 + 0.5) / 1000;
	my $fpr575x4pct = int ($idprogramfprless575{$id}{'g2p'} / $idprogramtemplates{$id}{'g2p'} * 100000 + 0.5) / 1000;
	my $pssmtropism = my $fpr2tropism = my $fpr575tropism = my $g2ptropism = "R5";
	if ($pssmx4pct > 2) {
		$pssmtropism = "X4";
	}
	if ($fpr2x4pct > 2) {
		$fpr2tropism = "X4";
	}
	if ($fpr575x4pct > 2) {
		$fpr575tropism = "X4";
	}
	if ($fpr2tropism eq "X4" and $fpr575tropism eq "X4") {
		$g2ptropism = "X4";
	}elsif ($fpr2tropism eq "R5" and $fpr575tropism eq "X4") {
		$g2ptropism = "INT";
	}
	print OUT "$id,$idprogramtotaltemplates{$id}{'PSSM'},$idprogramtemplates{$id}{'PSSM'},$idprogramx4{$id}{'PSSM'},$pssmx4pct,$pssmtropism,$idprogramfprless2{$id}{'g2p'},$fpr2x4pct,$fpr2tropism,$idprogramfprless575{$id}{'g2p'},$fpr575x4pct,$fpr575tropism,$g2ptropism\n";
}
close OUT;

print "total $count samples\n\n";

