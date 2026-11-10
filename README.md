# VvECPG1 manuscript analysis code

Code for the grapevine *VvECPG1* RNA-seq, tissue atlas, hormone-response, co-expression, functional-context and protein-homology results. Focal locus: **Vitvi05_01chr01g05750**, PN40024 T2T v5.1. Original code is licensed under [MIT](LICENSE).

This repository contains the manuscript analysis code and its dependencies. Superseded exploratory networks, duplicate plotting implementations, internal report generators and packaging history are excluded. The complete reproducibility ZIP retains the original inputs, checkpoints, database/version records and historical provenance separately.

## Install and check

Python **3.11+** is required for integrity/restoration utilities. Main R analysis dependencies are data.table, Matrix, ggplot2, patchwork and ggrepel. Optional Python figure/atlas/alignment dependencies are listed in [requirements-optional.txt](requirements-optional.txt). Original Linux upstream environment records are in `upstream/envs/`.

```sh
# No data or R installation needed for repository integrity checking.
python3 scripts/check_repository.py

# Also check the main R environment.
bash run.sh check
```

## Obtain the analysis data

Download the full ZIP supplied with the release, then restore it from the repository root:

```sh
python3 scripts/restore_data.py "/path/to/VvECPG1_submission_code_complete_2026-11-10.zip"
bash run.sh check-data
```

The expected archive size and SHA-256 are pinned in [data_archive.json](data_archive.json). Restoration checks the ZIP and all 2,301 archived input/source files, extracts only `inputs/`, `reference_results/` and `source_addenda/`, and preserves repository code and existing different files. Those data directories and generated `results/` are ignored by Git. The full archive contains supporting provenance beyond the minimum inputs for a single analysis; these files do not need to be committed to the repository. Allow approximately 1 GB for the data plus working space.

No public repository/release URL has been assigned. **For publication:** upload the contents of this small folder as the repository, and attach the existing 840-MB full ZIP and its `.sha256` sidecar to a GitHub release. Add the release link here and in `data_archive.json`, and add the manuscript citation/DOI when assigned. Git or GitHub Desktop preserves the supplied executable permissions. No upload has been performed by this packaging task.

## Reproduce the manuscript analyses

After data restoration:

```sh
bash run.sh verify-reference     # Numerical co-expression/functional/hormone checks
bash run.sh figures-reference   # Final atlas and combined network/function panels
bash run.sh hormone-figure      # Ten 24-h focal tests and normalized-count plot
bash run.sh atlas-rebuild       # Atlas quantities from individual Salmon archives
bash run.sh programme-rebuild   # Programme/random-set statistics from saved DE tables
bash run.sh homology-features   # EA-box features from archived HMM/predictor outputs
bash run.sh reproduce           # Full processed-input discovery/functional permutations
```

Optional commands are `paralog-reference` and `homology-similarity` (Biopython). Set `PYTHON` or `RSCRIPT` to a compatible executable. Fresh raw-read processing, DESeq2/signature fitting, technical adjustment, programme bootstrap and homology searches are described in [upstream/README.md](upstream/README.md). [MANUSCRIPT_CODE_MAP.md](MANUSCRIPT_CODE_MAP.md) maps results to entry points.

`metadata/` provides small, browsable final sample/atlas metadata copies. Execution uses the checksum-pinned copies restored under `inputs/`. `provenance/input_manifest.tsv` is required by the restoration/checking utilities. [validation/README.md](validation/README.md) distinguishes completed checks from unexecuted fresh analyses.

## Scientific choices

- **Atlas:** 624 runs / 591 BioSamples / 48 studies in the primary view; 630 / 597 / 49 with six inferred-callus sensitivity samples. SRR39573501–SRR39573512 are all **berry / fruit**. BioSample medians precede study medians and medians of study medians; cross-study differences are descriptive.
- **Controls/hormones:** three biological untreated controls collected at 24 h were prepared in both batches. Only the ten time-matched 24-h focal Wald tests form the primary focal BH family. Other times are descriptive. Genome-wide BH, legacy 40-test focal BH, MLE fold changes and apeglm estimates remain distinct. The plot uses within-batch normalized counts.
- **Discovery:** technical adjustment fitted across 126 libraries; primary within-group Pearson discovery uses 90 treated 6/12/24-h samples. Selection retains the original 1,000 permutations, 200 group bootstraps, 80% stability and signed leave-one-hormone-out threshold. Independent evaluation uses 27 samples/nine conditions, 10,000 within-condition directional permutations and BH across fixed measurable candidates.
- **Interpretation:** 106 partners replicate, all positive; broad independent positive correlations limit specificity. The positive embryogenic signature is central; pair-specific conditional scores exclude focal/tested genes. No GO BP or KEGG set passes corrected cross-dataset conjunction testing. DNA-maintenance enrichment depends on background, and transferred annotations do not establish function. Partner-to-partner edges are descriptive.
- **Protein:** EA-box features and negative localization predictions do not establish peptide function or cellular location. The feature-class heuristic agrees with 11/16 published representatives, with disagreements retained. Focal and paralog correlation measures must be distinguished across datasets and methods.

This repository covers these computational analyses. Transformation, pigment scoring, editing and qRT-PCR source data/statistics, and the separate seven-cultivar RNA-seq project's original pipeline, remain separate submission components. Third-party data/software retain their existing terms; the MIT code license does not assign a license to research datasets.
