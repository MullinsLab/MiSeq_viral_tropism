#!/usr/bin/make -f

SHELL := /bin/bash
BINPATH := ../../bin
AMPREF := ../env_amplicon_HXB2.fas
MLEN := 300
MIN := 12
MAX := 12

all : $(sample)_fragment_size.txt
$(sample)_R1_sickle.fastq $(sample)_R2_sickle.fastq : $(sample)_R1.fastq $(sample)_R2.fastq
	sickle pe -f $(sample)_R1.fastq -r $(sample)_R2.fastq -o $(sample)_R1_sickle.fastq -p $(sample)_R2_sickle.fastq -q 20 -w 10 -n -l 75 -t sanger -s $(sample)_sickle_single.fastq > $(sample).log
$(sample)_R1_sickle_cutadapt.fastq $(sample)_R2_sickle_cutadapt.fastq : $(sample)_R1_sickle.fastq $(sample)_R2_sickle.fastq
	cutadapt -g GCACAGTACAATGTACACATGGAATTA -G CCGCTCCGTCCGACGACTCACTATA --trimmed-only -o $(sample)_R1_sickle_cutadapt.fastq -p $(sample)_R2_sickle_cutadapt.fastq $^ >> $(sample).log
$(sample)_sickle_cutadapt_pear_all.assembled.fastq : $(sample)_R1_sickle_cutadapt.fastq $(sample)_R2_sickle_cutadapt.fastq
	pear -f $(sample)_R1_sickle_cutadapt.fastq -r $(sample)_R2_sickle_cutadapt.fastq -o $(sample)_sickle_cutadapt_pear_all -g 2 >> $(sample).log
	rm *.discarded.fastq *.unassembled.*
$(sample)_sickle_cutadapt_pear.fastq : $(sample)_sickle_cutadapt_pear_all.assembled.fastq
	$(BINPATH)/retrieve_merged_reads.pl $^ $(MLEN) $@ >> $(sample).log
$(sample)_sickle_cutadapt_pear_all_rt.fastq : $(sample)_sickle_cutadapt_pear.fastq
	cutadapt -a TTAATTGTGGAGGGGAATTTTTCTACTG --trimmed-only -o $@ $(sample)_sickle_cutadapt_pear_all.assembled.fastq >> $(sample).log
$(sample)_sickle_cutadapt_pear_rt.fastq : $(sample)_sickle_cutadapt_pear_all_rt.fastq
	cutadapt -a TTAATTGTGGAGGGGAATTTTTCTACTG --trimmed-only -o $@ $(sample)_sickle_cutadapt_pear.fastq >> $(sample).log
#$(sample)_sickle_cutadapt_pear.assembled.fastq : $(sample)_sickle_cutadapt_pear_all_rt.fastq
#	pear -f $(sample)_R1_sickle_cutadapt.fastq -r $(sample)_R2_sickle_cutadapt.fastq -o $(sample)_sickle_cutadapt_pear -g 2 -n 300 >> $(sample).log
$(sample)_sickle_cutadapt_pear_rt_tnp.fasta : $(sample)_sickle_cutadapt_pear.fastq $(sample)_sickle_cutadapt_pear_rt.fastq 
	$(BINPATH)/template_primer_fasta.pl $^ $@ >> $(sample).log
$(sample)_sickle_cutadapt_pear_rt_t.fasta : $(sample)_sickle_cutadapt_pear_rt_tnp.fasta
	cutadapt -g TTAATTGTGGAGGGGAATTTTTCTACTG --trimmed-only -o $@ $^ >> $(sample).log
$(sample)_templates.txt : $(sample)_sickle_cutadapt_pear_rt_t.fasta	
	$(BINPATH)/retrieve_templateid_from_tnpfile.pl $^ $@ >> $(sample).log
template_sequences : $(sample)_templates.txt $(sample)_sickle_cutadapt_pear_rt.fastq
	$(BINPATH)/retrieve_sequences_per_template_family_gt2.pl $^ >> $(sample).log
$(sample)_consensus.fasta : template_sequences
	$(BINPATH)/cal_consensus_qual_from_dir.pl $^ $@ >> $(sample).log
$(sample)_variants.fasta : $(sample)_consensus.fasta
	$(BINPATH)/unique_consensus.pl $^ $@ >> $(sample).log	
$(sample)_v3.fasta : $(sample)_variants.fasta $(AMPREF)
	$(BINPATH)/extract_v3.pl $^ $@ >> $(sample).log
#ifdef $(MIN)
#	ifdef $(MAX)
$(sample)_templates_family_size.txt : $(sample)_v3.fasta
	$(BINPATH)/template_family_size_distribution.pl $(sample)_templates.txt $@ $(MIN) $(MAX) >> $(sample).log
#else
#$(sample)_templates_family_size.txt : $(sample)_variants.fasta
#	$(BINPATH)/template_family_size_distribution.pl $(sample)_templates.txt $@
#endif
$(sample)_templates_size.txt : $(sample)_templates_family_size.txt
	$(BINPATH)/template_size_distribution.pl $(sample)_templates.txt $@ >> $(sample).log
$(sample)_fragment_size.txt : $(sample)_templates_size.txt
	$(BINPATH)/fragment_size_distribution.pl $(sample)_sickle_cutadapt_pear_all_rt.fastq $@ >> $(sample).log

cleanstat :
	rm *_size.txt *_variants.fasta *_v3.fasta
clean : 
	rm *_sickle*.* *.log *.fasta *.txt *.aln
	rm -r template_sequences
