# Manuscript analysis-to-code map

| Manuscript content | Reproduction entry point | Data / provenance |
| --- | --- | --- |
| Hormone-read QC, alignment, UMI deduplication and integer counting | Staged `make primary` | A driver, sample/reference manifests and executed commands; external raw reads/reference binaries |
| Within-batch differential expression | Staged `make de` | C complete 40-table checkpoint, model objects, size factors; primary counts and metadata |
| Ten 24-h focal hormone responses / Supplementary hormone panel | `bash run.sh hormone-figure` | C labelled within-batch focal counts and Wald tests; ten-test BH computed separately from genome-wide and legacy 40-test families |
| Curated tissue atlas / focal_atlas panels | `bash run.sh atlas-rebuild`, then existing `figures-reference` for manuscript panel layout | B individual quantifications plus final `inputs/atlas` metadata; `figures-reference` intentionally reads frozen final summaries; reconstructed summaries are independently checked against them |
| Embryogenic signature construction | Staged `make signature` | B transcript NumReads/tx2gene, PRJNA552475 and PRJNA980588, original DESeq2 designs |
| Positive-signature programme context and sensitivities | `bash run.sh programme-rebuild`; staged programme/bootstrap targets | C all-gene shrunken estimates, fixed positive signature, adjusted matrix/blind VST; optional bootstrap refits require DESeq2/apeglm |
| Technical adjustment of all 126 libraries | Staged `make technical-adjustment` | Original isolated context-C fit, primary counts, QC and metadata; comparison with archived adjusted matrix |
| Primary discovery (90 samples), fixed independent validation (27 samples), conditional sensitivities | `bash run.sh reproduce` | `scripts/coexpression/run.R`; frozen processed inputs, original seeds/permutation/bootstrap counts |
| Replicated functional interpretation / combined S4 B–C | `scripts/functional/run.R`, `scripts/figures/combined_S4.R` via runner | Frozen GO/KEGG/homolog annotations; original backgrounds and corrected conjunction analysis |
| Replicated focal neighborhood / combined S4 A | `scripts/figures/network.R` via `figures-reference` or `reproduce` | Replicated candidates ranked by smaller positive correlation; post-selection descriptive partner edges |
| Paralog expression/correlation | `bash run.sh paralog-reference` | Separately labelled legacy bicor and primary discovery/independent Pearson results |
| EA-box protein features/classes | `bash run.sh homology-features` / `homology-similarity` | D frozen HMM alignments and SignalP/DeepTMHMM outputs; Biopython for global alignment |
| Full-length homology, legacy locus, tandem paralog, Vitaceae comparisons | Staged homology targets; E `ANALYSIS_NOTES.md` | E original BLAST/Pfam/CDD outputs and corrected supporting tables; fresh searches require recorded external databases/genomes |

This package covers the RNA-seq, atlas, co-expression, functional-context and homology analyses. It does not reconstruct transformation, pigment scoring, editing or qRT-PCR experiments from figures, nor the separate seven-cultivar RNA-seq project's original pipeline. Those authors' source data/workflows remain separate submission components.
