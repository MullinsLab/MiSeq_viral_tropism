SHELL := /bin/bash
export SHELLOPTS:=errexit:pipefail
.DELETE_ON_ERROR:

BIN := bin
DATA := data
SDATA := $(DATA)/$(sample)
OUTPUT := Outputs
AMPREF := refs/amplicon.fas
FPFILE := refs/fprimer.fas
RPFILE := refs/rprimer.fas
RTFILE := refs/RTprimer.fas
MINLEN := 300
TEMPID := 12
MIN := 12
MAX := 12
TFS ?= 0	# template family size cutoff - default to use cutoff model to determine template family size cutoff to calculate consensus
SM ?= 0.5		# supermajority cutoff - default to calculate consensus using super majority > 0.5
ER ?= 0.01	# overall error rate cutoff (default 0.01) - used to calculate template family size cutoff (the other allowed values are 0.005 and 0.02)
V3 ?= 1.2		# V3 sequences within range of V3 loop size to keep. 0 means no range, including all (value between 1 and 1.5)

TEXT := $(shell bin/parse_refs.pl $(FPFILE) $(RPFILE) $(RTFILE))
RT := $(word 1, $(TEXT))
FP := $(word 2, $(TEXT))
RP := $(word 3, $(TEXT))
MINOVLP := $(word 4, $(TEXT))
RTOVLP := $(word 5, $(TEXT))

all : $(OUTPUT)/$(sample)/$(sample)_fragment_size.txt
# sickle quality trimming
$(OUTPUT)/$(sample)/$(sample)_R1_sickle.fastq $(OUTPUT)/$(sample)/$(sample)_R2_sickle.fastq : $(DATA)/$(sample)_R1.fastq.gz $(DATA)/$(sample)_R2.fastq.gz | mkdir-output-$(sample)
	$(BIN)/sickle pe -q 20 -w 10 -l 75 -n -t sanger \
	-f <(gunzip -c $(DATA)/$(sample)_R1.fastq.gz) \
	-r <(gunzip -c $(DATA)/$(sample)_R2.fastq.gz) \
	-o $(OUTPUT)/$(sample)/$(sample)_R1_sickle.fastq \
	-p $(OUTPUT)/$(sample)/$(sample)_R2_sickle.fastq \
	-s $(OUTPUT)/$(sample)/$(sample)_sickle_single.fastq \
	> $(OUTPUT)/$(sample)/$(sample).log
# cutadapt trimming 2nd round primers 
$(OUTPUT)/$(sample)/$(sample)_R1_sickle_cutadapt.fastq $(OUTPUT)/$(sample)/$(sample)_R2_sickle_cutadapt.fastq : $(OUTPUT)/$(sample)/$(sample)_R1_sickle.fastq $(OUTPUT)/$(sample)/$(sample)_R2_sickle.fastq
	cutadapt --trimmed-only -O $(MINOVLP) \
	-g $(FP) \
	-G $(RP) \
	-o $(OUTPUT)/$(sample)/$(sample)_R1_sickle_cutadapt.fastq \
	-p $(OUTPUT)/$(sample)/$(sample)_R2_sickle_cutadapt.fastq \
	$^ >> $(OUTPUT)/$(sample)/$(sample).log
# pear merging trimmed paired reads
$(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_all.assembled.fastq : $(OUTPUT)/$(sample)/$(sample)_R1_sickle_cutadapt.fastq $(OUTPUT)/$(sample)/$(sample)_R2_sickle_cutadapt.fastq
	$(BIN)/pear -g 2 \
	-f $(OUTPUT)/$(sample)/$(sample)_R1_sickle_cutadapt.fastq \
	-r $(OUTPUT)/$(sample)/$(sample)_R2_sickle_cutadapt.fastq \
	-o $(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_all \
	>> $(OUTPUT)/$(sample)/$(sample).log	
# retrieve merged reads with length >= MINLEN
$(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear.fastq : $(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_all.assembled.fastq
	$(BIN)/retrieve_merged_reads.pl $^ $(MINLEN) $@ >> $(OUTPUT)/$(sample)/$(sample).log
# trim RT+templatedid in all merged reads
$(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_all_rt.fastq : $(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear.fastq
	cutadapt --trimmed-only -O $(RTOVLP) \
	-a $(RT) \
	-o $@ \
	$(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_all.assembled.fastq >> $(OUTPUT)/$(sample)/$(sample).log
# trim RT+templatedid in merged reads with length >= MINLEN
$(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_rt.fastq : $(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_all_rt.fastq
	cutadapt --trimmed-only -O $(RTOVLP) \
	-a $(RT) \
	-o $@ \
	$(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear.fastq >> $(OUTPUT)/$(sample)/$(sample).log
# retrieve only RT+templatedid sequences
$(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_rt_tnp.fasta : $(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear.fastq $(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_rt.fastq 
	$(BIN)/template_primer_fasta.pl $^ $@ >> $(OUTPUT)/$(sample)/$(sample).log
# trim RT from RT+templateid sequences
$(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_rt_t.fasta : $(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_rt_tnp.fasta
	cutadapt --trimmed-only -O $(RTOVLP) \
	-g $(RT) \
	-o $@ \
	$^ >> $(OUTPUT)/$(sample)/$(sample).log
# collaps each template id and link merged reads to templated id
$(OUTPUT)/$(sample)/$(sample)_templates.txt : $(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_rt_t.fasta	
	$(BIN)/retrieve_templateid_from_tnpfile.pl $^ $@ >> $(OUTPUT)/$(sample)/$(sample).log
# retrieve merged reads for each template id with the length of template id = TEMPID (here 12) and family size of merged reads >= $(TFS)
$(OUTPUT)/$(sample)/template_sequences.log : $(OUTPUT)/$(sample)/$(sample)_templates.txt $(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_rt.fastq
	$(BIN)/retrieve_sequences_per_template_family_cutoff.pl $^ $(TEMPID) $(ER) $(TFS) >> $(OUTPUT)/$(sample)/$(sample).log
# calculate consensus sequence for each template id
$(OUTPUT)/$(sample)/$(sample)_consensus.fasta : $(OUTPUT)/$(sample)/template_sequences.log
	$(BIN)/cal_consensus_supermajority_from_dir.pl $(OUTPUT)/$(sample)/template_sequences $@ $(SM) >> $(OUTPUT)/$(sample)/$(sample).log
# collaps template consensus to variants
$(OUTPUT)/$(sample)/$(sample)_variants.fasta : $(OUTPUT)/$(sample)/$(sample)_consensus.fasta
	$(BIN)/unique_consensus.pl $^ $@ >> $(OUTPUT)/$(sample)/$(sample).log	
# extract V3 loop from variants and collaps to V3 variants	
$(OUTPUT)/$(sample)/$(sample)_v3.fasta : $(OUTPUT)/$(sample)/$(sample)_variants.fasta $(AMPREF)
	$(BIN)/extract_v3.pl $^ $@ $(V3) >> $(OUTPUT)/$(sample)/$(sample).log
# calculate template family size distributions for template id length between MIN and MAX (here 12)
$(OUTPUT)/$(sample)/$(sample)_templates_family_size.txt : $(OUTPUT)/$(sample)/$(sample)_v3.fasta
	$(BIN)/template_family_size_distribution.pl $(OUTPUT)/$(sample)/$(sample)_templates.txt $@ $(MIN) $(MAX) >> $(OUTPUT)/$(sample)/$(sample).log
# calculate template id length distributions
$(OUTPUT)/$(sample)/$(sample)_templates_size.txt : $(OUTPUT)/$(sample)/$(sample)_templates_family_size.txt
	$(BIN)/template_size_distribution.pl $(OUTPUT)/$(sample)/$(sample)_templates.txt $@ >> $(OUTPUT)/$(sample)/$(sample).log
# calculate distribution of all merged reads
$(OUTPUT)/$(sample)/$(sample)_fragment_size.txt : $(OUTPUT)/$(sample)/$(sample)_templates_size.txt
	$(BIN)/fragment_size_distribution.pl $(OUTPUT)/$(sample)/$(sample)_sickle_cutadapt_pear_all_rt.fastq $@ >> $(OUTPUT)/$(sample)/$(sample).log

	rm $(OUTPUT)/$(sample)/*_sickle*.*

mkdir-output-%:
	mkdir -p $(OUTPUT)/$*
	
clean : 
	rm -rfv $(OUTPUT)/
	

