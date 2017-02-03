SHELL := /bin/bash
export SHELLOPTS:=errexit:pipefail
.DELETE_ON_ERROR:

BIN := bin
DATA := data
SDATA := $(DATA)/$(sample)
AMPREF := refs/amplicon.fas
FPFILE := refs/fprimer.fas
RPFILE := refs/rprimer.fas
RTFILE := refs/RTprimer.fas
MINLEN := 300
TEMPID := 12
MIN := 12
MAX := 12

TEXT := $(shell bin/parse_refs.pl $(FPFILE) $(RPFILE) $(RTFILE))
RT := $(word 1, $(TEXT))
FP := $(word 2, $(TEXT))
RP := $(word 3, $(TEXT))
MINOVLP := $(word 4, $(TEXT))
RTOVLP := $(word 5, $(TEXT))

all : $(SDATA)/$(sample)_fragment_size.txt
# sickle quality trimming
$(SDATA)/$(sample)_R1_sickle.fastq $(SDATA)/$(sample)_R2_sickle.fastq : $(SDATA)/$(sample)_R1.fastq.gz $(SDATA)/$(sample)_R2.fastq.gz
	$(BIN)/sickle pe -q 20 -w 10 -l 75 -n -t sanger \
	-f <(gunzip -c $(SDATA)/$(sample)_R1.fastq.gz) \
	-r <(gunzip -c $(SDATA)/$(sample)_R2.fastq.gz) \
	-o $(SDATA)/$(sample)_R1_sickle.fastq \
	-p $(SDATA)/$(sample)_R2_sickle.fastq \
	-s $(SDATA)/$(sample)_sickle_single.fastq \
	> $(SDATA)/$(sample).log
# cutadapt trimming 2nd round primers 
$(SDATA)/$(sample)_R1_sickle_cutadapt.fastq $(SDATA)/$(sample)_R2_sickle_cutadapt.fastq : $(SDATA)/$(sample)_R1_sickle.fastq $(SDATA)/$(sample)_R2_sickle.fastq
	$(BIN)/cutadapt --trimmed-only -O $(MINOVLP) \
	-g $(FP) \
	-G $(RP) \
	-o $(SDATA)/$(sample)_R1_sickle_cutadapt.fastq \
	-p $(SDATA)/$(sample)_R2_sickle_cutadapt.fastq \
	$^ >> $(SDATA)/$(sample).log
# pear merging trimmed paired reads
$(SDATA)/$(sample)_sickle_cutadapt_pear_all.assembled.fastq : $(SDATA)/$(sample)_R1_sickle_cutadapt.fastq $(SDATA)/$(sample)_R2_sickle_cutadapt.fastq
	$(BIN)/pear -g 2 \
	-f $(SDATA)/$(sample)_R1_sickle_cutadapt.fastq \
	-r $(SDATA)/$(sample)_R2_sickle_cutadapt.fastq \
	-o $(SDATA)/$(sample)_sickle_cutadapt_pear_all \
	>> $(SDATA)/$(sample).log	
# retrieve merged reads with length >= MINLEN
$(SDATA)/$(sample)_sickle_cutadapt_pear.fastq : $(SDATA)/$(sample)_sickle_cutadapt_pear_all.assembled.fastq
	$(BIN)/retrieve_merged_reads.pl $^ $(MINLEN) $@ >> $(SDATA)/$(sample).log
# trim RT+templatedid in all merged reads
$(SDATA)/$(sample)_sickle_cutadapt_pear_all_rt.fastq : $(SDATA)/$(sample)_sickle_cutadapt_pear.fastq
	cutadapt --trimmed-only -O $(RTOVLP) \
	-a $(RT) \
	-o $@ \
	$(SDATA)/$(sample)_sickle_cutadapt_pear_all.assembled.fastq >> $(SDATA)/$(sample).log
# trim RT+templatedid in merged reads with length >= MINLEN
$(SDATA)/$(sample)_sickle_cutadapt_pear_rt.fastq : $(SDATA)/$(sample)_sickle_cutadapt_pear_all_rt.fastq
	cutadapt --trimmed-only -O $(RTOVLP) \
	-a $(RT) \
	-o $@ \
	$(SDATA)/$(sample)_sickle_cutadapt_pear.fastq >> $(SDATA)/$(sample).log
# retrieve only RT+templatedid sequences
$(SDATA)/$(sample)_sickle_cutadapt_pear_rt_tnp.fasta : $(SDATA)/$(sample)_sickle_cutadapt_pear.fastq $(SDATA)/$(sample)_sickle_cutadapt_pear_rt.fastq 
	$(BIN)/template_primer_fasta.pl $^ $@ >> $(SDATA)/$(sample).log
# trim RT from RT+templateid sequences
$(SDATA)/$(sample)_sickle_cutadapt_pear_rt_t.fasta : $(SDATA)/$(sample)_sickle_cutadapt_pear_rt_tnp.fasta
	cutadapt --trimmed-only -O $(RTOVLP) \
	-g $(RT) \
	-o $@ \
	$^ >> $(SDATA)/$(sample).log
# collaps each template id and link merged reads to templated id
$(SDATA)/$(sample)_templates.txt : $(SDATA)/$(sample)_sickle_cutadapt_pear_rt_t.fasta	
	$(BIN)/retrieve_templateid_from_tnpfile.pl $^ $@ >> $(SDATA)/$(sample).log
# retrieve merged reads for each template id with the length of template id = TEMPID (here 12) and family size of merged reads >= 3
$(SDATA)/template_sequences.log : $(SDATA)/$(sample)_templates.txt $(SDATA)/$(sample)_sickle_cutadapt_pear_rt.fastq
	$(BIN)/retrieve_sequences_per_template_family_gt2.pl $^ $(TEMPID) >> $(SDATA)/$(sample).log
# calculate consensus sequence for each template id
$(SDATA)/$(sample)_consensus.fasta : $(SDATA)/template_sequences.log
	$(BIN)/cal_consensus_qual_from_dir.pl $(SDATA)/template_sequences $@ >> $(SDATA)/$(sample).log
# collaps template consensus to variants
$(SDATA)/$(sample)_variants.fasta : $(SDATA)/$(sample)_consensus.fasta
	$(BIN)/unique_consensus.pl $^ $@ >> $(SDATA)/$(sample).log	
# extract V3 loop from variants and collaps to V3 variants	
$(SDATA)/$(sample)_v3.fasta : $(SDATA)/$(sample)_variants.fasta $(AMPREF)
	$(BIN)/extract_v3.pl $^ $@ >> $(SDATA)/$(sample).log
# calculate template family size distributions for template id length between MIN and MAX (here 12)
$(SDATA)/$(sample)_templates_family_size.txt : $(SDATA)/$(sample)_v3.fasta
	$(BIN)/template_family_size_distribution.pl $(SDATA)/$(sample)_templates.txt $@ $(MIN) $(MAX) >> $(SDATA)/$(sample).log
# calculate template id length distributions
$(SDATA)/$(sample)_templates_size.txt : $(SDATA)/$(sample)_templates_family_size.txt
	$(BIN)/template_size_distribution.pl $(SDATA)/$(sample)_templates.txt $@ >> $(SDATA)/$(sample).log
# calculate distribution of all merged reads
$(SDATA)/$(sample)_fragment_size.txt : $(SDATA)/$(sample)_templates_size.txt
	$(BIN)/fragment_size_distribution.pl $(SDATA)/$(sample)_sickle_cutadapt_pear_all_rt.fastq $@ >> $(SDATA)/$(sample).log
	rm -rfv $(SDATA)/template_sequences/
	rm $(SDATA)/*_sickle*.*
	
clean : 
	rm $(SDATA)/*.log $(SDATA)/*.fasta $(SDATA)/*.txt $(SDATA)/*.aln 
	

