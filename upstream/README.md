# Upstream execution

From the repository root, restore the full data ZIP, then:

```sh
bash run.sh prepare-upstream
cd results/upstream_workspace
```

This stages a writable original-layout project with counts, metadata, QC, signature inputs, all-gene DE tables, fitted models, frozen homology outputs and a link to the Salmon quantifications. Existing staged results are preserved. The annotation-only lookup recovered from a saved DE table supports labels; it does not reconstruct the entire historical crosswalk project.

## Refit or reconstruct processed-input stages

```sh
make qc                    # All-library normalization/blind VST and QC figures
make technical-adjustment  # Original isolated all-126-library context-C fit
make de                    # Within-batch DESeq2 ~group, apeglm and diagnostics
make signature             # Signature discovery from atlas NumReads
make programme             # Programme/contrast/random-set calculations
make programme-figures PYTHON=python3
```

Fresh fitting requires the original compatible R/Bioconductor environment: DESeq2, apeglm, edgeR, limma and the other imports recorded in source/environment files. `technical-adjustment` writes `context_C_rebuilt.rds` and compares it with the frozen checkpoint. From the repository root, use the new fit with:

```sh
VVECPG1_ADJUSTED_RDS=results/upstream_workspace/results/pn40024_v5.1/focal_network/context_C_rebuilt.rds bash run.sh reproduce
```

The isolated transformation copies the original context-C calculations. Primary discovery uses the current 90-sample Pearson workflow, not the superseded bicor candidate calls. The historical full network script remains in addendum A of the data archive for provenance.

The optional original programme confidence-interval calculation requires expensive DESeq2/apeglm refits:

```sh
make bootstrap-check
# Run TASK=1 through TASK=50 once each: four refits per task = 200 replicates.
make bootstrap-task TASK=1
# After all tasks finish:
make bootstrap-summary
make programme-figures PYTHON=python3
```

Ten 24-h focal tests are time-matched. All 40 saved contrasts are retained for provenance/descriptive programme calculations. The final hormone plotting command at the repository root uses the ten-test family and within-batch replicate counts; the superseded plot with MLE labels is excluded from the working code.

## Raw-read processing

Restore raw FASTQs, vendor calibration BAMs, exact reference FASTAs/GFFs and compiled indexes using the manifests in `source_addenda/A_upstream_workflow/checksums/`. Original executed commands, software versions and reference source records are retained in A. Large raw data/index binaries are external resources.

Use the supplied Linux environment specs/locks in `upstream/envs/`; restore the RNA-seq tool environment at staged `envs/rnaseq` because the Snakefile invokes tools there. Run within an appropriate allocation:

```sh
make references             # Download/prepare references; verify against recorded hashes
make dryrun                 # Snakemake graph; requires restored raw/reference files
make primary CORES=16 MEM_MB=64000
```

The primary target runs fastp → STAR → UMIcollapse → featureCounts/QC and includes historical vendor-BAM calibration recounts, so those BAMs are needed for the exact target. The provided Makefile runs locally within its allocation; original SLURM profiles are retained. Only the documented primary targets are supported; unused historical optional genome branches in the original config are not manuscript analyses.

Atlas classification/selection and per-run Salmon scripts are retained. Use the final metadata under `inputs/atlas`, rather than re-querying ENA or silently replacing frozen curation. B's `salmon_index/BUILD_COMMANDS.sh` and manifests specify the exact decoy-aware reference/index. Final aggregation/plotting is done by `scripts/atlas/rebuild.py` and `scripts/figures/atlas.R`; the older duplicate aggregators are excluded. The added `ranged_fetch.sh` is an existing reference-download helper that was missing from the earlier package.

## Homology

From the repository root, `homology-features` reconstructs features from frozen HMM/SignalP/DeepTMHMM outputs. `homology-similarity` additionally repeats the global BLOSUM62 alignments using Biopython. `python3 scripts/homology/classify.py --rerun-hmm` repeats local HMM alignment with HMMER 3.4 and the archived PF28138 model, while preserving frozen predictor calls.

For full searches, use staged targets `homology-blastp`, `homology-domains`, `homology-summary`, `homology-supporting` in the step order documented in `source_addenda/E_vvecpg1_homology/ANALYSIS_NOTES.md`. Restore the recorded database/genome versions. Set `BLAST_BIN`, `HMMER_BIN` and `DATASETS_BIN` if needed; set `BLASTDB` to the recorded nr/refseq directory and `VVECPG1_LEGACY_FASTA` to the 12X.v2 FASTA. The domain script requires the recorded Pfam 39.0/CDD snapshots, not current downloads. Fresh SignalP/DeepTMHMM runs require their model software/access; historical launch scripts and frozen predictions remain in D.
