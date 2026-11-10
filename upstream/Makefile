# Focused entry points. Run after `bash run.sh prepare-upstream` from the staged workspace.
SHELL := /bin/bash
.SHELLFLAGS := -euo pipefail -c
ROOT := $(abspath .)
PYTHON ?= python3
RSCRIPT ?= Rscript
CORES ?= 16
MEM_MB ?= 64000
SMK ?= snakemake
.PHONY: dryrun primary references de qc signature programme bootstrap-check bootstrap-task bootstrap-summary programme-figures
# `envs/rnaseq` must exist at the staged root: Snakefile calls tools from that prefix.
dryrun:
	$(SMK) -s workflow/Snakefile -n --cores $(CORES) primary
primary:
	$(SMK) -s workflow/Snakefile --cores $(CORES) --resources slots=$(CORES) mem_mb=$(MEM_MB) primary
references:
	bash scripts/download_references.sh "$(ROOT)"
	bash scripts/reference_checks.sh "$(ROOT)"
de:
	$(RSCRIPT) scripts/de_analysis.R config/config.yaml
	$(RSCRIPT) scripts/de_plots.R config/config.yaml
qc:
	$(RSCRIPT) scripts/qc_plots.R config/config.yaml
	$(PYTHON) scripts/focal_locus_plot.py --root "$(ROOT)"
signature:
	VVECPG1_UPSTREAM_ROOT="$(ROOT)" $(RSCRIPT) scripts/embryogenic_signature.R
programme:
	VVECPG1_UPSTREAM_ROOT="$(ROOT)" $(RSCRIPT) scripts/focal_24h_programme.R
bootstrap-check:
	VVECPG1_UPSTREAM_ROOT="$(ROOT)" $(RSCRIPT) scripts/focal_24h_programme_boot.R 0 1
# 50 tasks x 4 refits = original 200 replicates; task 0 is the unresampled check.
bootstrap-task:
	@test -n "$(TASK)" || { echo 'Supply TASK=1 through 50'; exit 2; }
	VVECPG1_UPSTREAM_ROOT="$(ROOT)" $(RSCRIPT) scripts/focal_24h_programme_boot.R $(TASK) 4
bootstrap-summary:
	VVECPG1_UPSTREAM_ROOT="$(ROOT)" $(RSCRIPT) scripts/focal_24h_programme_boot_summary.R
programme-figures:
	VVECPG1_UPSTREAM_ROOT="$(ROOT)" $(PYTHON) scripts/focal_24h_programme_figures.py

# Requires restored databases/genomes identified in source_addenda/E.
.PHONY: homology-blastp homology-domains homology-supporting homology-summary
homology-blastp:
	VVECPG1_UPSTREAM_ROOT="$(ROOT)" bash scripts/vvecpg1_homology_blastp.sh
homology-domains:
	VVECPG1_UPSTREAM_ROOT="$(ROOT)" bash scripts/vvecpg1_homology_step3_domains.sh
homology-supporting:
	VVECPG1_UPSTREAM_ROOT="$(ROOT)" $(PYTHON) scripts/vvecpg1_homology_supporting.py
homology-summary:
	VVECPG1_UPSTREAM_ROOT="$(ROOT)" $(PYTHON) scripts/vvecpg1_tblastn_summarise.py

.PHONY: technical-adjustment
technical-adjustment:
	VVECPG1_UPSTREAM_ROOT="$(ROOT)" $(RSCRIPT) scripts/rebuild_technical_adjustment.R
