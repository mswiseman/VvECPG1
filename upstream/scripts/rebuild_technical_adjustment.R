#!/usr/bin/env Rscript
# Isolated context-C transformation from focal_coexpression.R, unchanged math.
# No WGCNA modules, annotations, or legacy candidate calling are required.
# Requires DESeq2 and limma; run inside the staged upstream workspace.
suppressPackageStartupMessages({library(data.table);library(DESeq2);library(limma)})
ROOT <- Sys.getenv("VVECPG1_UPSTREAM_ROOT", normalizePath("."));setwd(ROOT)
set.seed(20261016)
cm <- fread("results/pn_v51/counts/pn_v51.primary.counts.tsv")
mat <- as.matrix(cm[, grep("^PPPJ8P_", names(cm)), with = FALSE]); rownames(mat) <- cm$gene_id
mat <- mat[!grepl("^MT:", rownames(mat)), ]; storage.mode(mat) <- "integer"
s <- fread("metadata/samples.tsv"); setkey(s, sample_uid)
s[, th := suppressWarnings(as.numeric(time_hours))]
s[, grp := ifelse(is_control, "control", condition)]
s[, grp_c := ifelse(is_control, control_label, condition)]          # bridging pairs share a level
s[, hormone := ifelse(is_control, "control", hormone)]
q <- fread("results/pn_v51/qc/pn_v51.per_sample_qc.tsv")
TECH <- c("log10_dedup_alignments", "dedup_retained_frac", "star_unique_frac", "frac_MT", "frac_chloroplast",
          "fc_sense_assigned_frac", "strand_sense_frac", "frac_within_500bp_upstream_of_end",
          "frac_downstream_of_annotated_end")
q[, log10_dedup_alignments := log10(dedup_output_alignments)]

make_ctx <- function(ids, combined = FALSE) {
  ss <- s[ids][order(sample_uid)]
  m <- mat[, ss$sample_uid]; m <- m[rowSums(m >= 10) >= 3, ]
  dgrp <- if (combined) ss$grp_c else ss$grp
  dds <- DESeqDataSetFromMatrix(m, data.frame(row.names = ss$sample_uid, g = factor(dgrp)), ~ g)
  dds <- estimateSizeFactors(dds)
  vs <- assay(vst(dds, blind = FALSE))
  design <- model.matrix(~ factor(dgrp))
  if (combined) vs <- removeBatchEffect(vs, batch = factor(ss$batch), design = design)
  tq <- merge(ss[, .(sample_uid)], q, by = "sample_uid", sort = FALSE)[match(ss$sample_uid, sample_uid)]
  tm <- as.matrix(tq[, ..TECH]); rownames(tm) <- ss$sample_uid
  tm <- cbind(tm, log_size_factor = log(sizeFactors(dds)), genes_detected = colSums(m > 0))
  if (combined) tm <- removeBatchEffect(t(tm), batch = factor(ss$batch), design = design) |> t()
  tres <- apply(tm, 2, function(v) resid(lm(v ~ factor(dgrp))))          # within-group technical variation only
  tres <- tres[, apply(tres, 2, sd) > 0, drop = FALSE]
  pc <- prcomp(tres, scale. = TRUE)
  ve <- cumsum(pc$sdev^2) / sum(pc$sdev^2); npc <- min(5, which(ve >= 0.8)[1])
  pcs <- pc$x[, seq_len(npc), drop = FALSE]
  adj <- removeBatchEffect(vs, covariates = pcs, design = design)
  # permutation blocks: replicate groups (bridging controls form one block in C)
  list(s = ss, grp = ss$grp, vs = vs, adj = adj, tech = tm, tres = tres, pcs = pcs, npc = npc, pc_var = ve[npc])
}

x <- make_ctx(s$sample_uid, combined = TRUE)
ref <- readRDS("results/pn40024_v5.1/focal_network/context_C_adjusted.rds")
stopifnot(identical(dimnames(x$adj),dimnames(ref$adj)))
delta <- max(abs(x$adj-ref$adj))
out <- "results/pn40024_v5.1/focal_network/context_C_rebuilt.rds"
saveRDS(list(adj=x$adj,s=x$s,technical_PCs=x$pcs),out)
writeLines(c(paste("max absolute matrix difference from frozen checkpoint:",delta),
             "The frozen checkpoint remains unchanged.",capture.output(sessionInfo())),
           "results/pn40024_v5.1/focal_network/technical_adjustment_rebuild.txt")
cat("Rebuilt context C; max absolute checkpoint discrepancy:",delta,"\n")
if(delta>1e-6)warning("Difference exceeds 1e-6; inspect package versions before substituting this matrix.")
