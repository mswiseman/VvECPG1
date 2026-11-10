#!/usr/bin/env Rscript
# Differential expression for PPPJ8P (PN40024 T2T v5.1 integer UMI-deduplicated molecule counts).
#
# Design (user, 2026-10-14, superseding the earlier PI note): C1, C2, C3 are the SAME three untreated
# callus samples in both batches. Plasmidsaurus re-prepared and re-sequenced them in batch 2 (OSIDs
# 124-126) as bridging samples to measure batch effects. So there are 3 biological controls, each with
# 2 libraries (one per batch). Within a batch the 3 controls are distinct biological samples.
# Controls were collected at 24 h. Only the ten 24-h contrasts are time-matched;
# comparisons at other times are descriptive, direction treated / control.
#
#  PRIMARY  : within-batch DESeq2 ~ group (batch 1: 31 groups vs controls 94-96; batch 2: 9 groups
#             vs controls 124-126). Wald test, BH across genes, apeglm-shrunken LFC. Each control
#             sample is used once per test, so the bridging design does not affect validity.
#  EDGER    : within-batch edgeR QL (glmQLFit / glmQLFTest) cross-check.
#  BRIDGE   : batch effect from the 3 bridging pairs: DESeq2 ~ control_label + batch on the 6
#             control libraries (paired; batch 2 / batch 1).
#  JOINT    : all 126 samples, limma-voom ~ batch + group with duplicateCorrelation, blocking on the
#             biological control sample (C1/C2/C3; treated libraries are their own blocks), so the 6
#             control libraries count as 3 correlated pairs, never as 6 independent replicates.
#             (An unblocked joint DESeq2 fit is invalid under this design and is not run.)
#  LEGACY FOCAL: all 40 contrasts retained for provenance; primary inference uses the
#             ten 24-h tests (scripts/hormone/summarize.py in the package). Historical family
#             from the per-contrast across-gene FDR); two-sided Wald, 95% Wald CI.
# Usage: Rscript scripts/de_analysis.R config/config.yaml
suppressPackageStartupMessages({
  library(yaml); library(data.table); library(DESeq2); library(edgeR); library(limma)
  library(BiocParallel); library(apeglm)
})
args <- commandArgs(trailingOnly = TRUE)
cfg <- yaml::read_yaml(ifelse(length(args), args[1], "config/config.yaml"))
stopifnot(isTRUE(cfg$design$controls_confirmed), cfg$design$control_model == "bridging_pairs")
setwd(cfg$project_root)
set.seed(20261014)
FOCAL <- "Vitvi05_01chr01g05750"
FDR <- cfg$design$fdr; LFC_T <- cfg$design$lfc_threshold
NCORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
register(MulticoreParam(NCORES))
out <- "results/pn40024_v5.1/de"; dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
for (d in c("joint_paired", "bridge_batch_effect")) dir.create(file.path(out, d), showWarnings = FALSE)
log_msg <- function(...) cat(format(Sys.time(), "%F %T"), ..., "\n")

# ---------------- inputs ----------------
cm <- fread("results/pn_v51/counts/pn_v51.primary.counts.tsv")
ann <- cm[, .(gene_id, Chr, Start, End, Strand, Length)]
mat <- as.matrix(cm[, grep("^PPPJ8P_", names(cm)), with = FALSE]); rownames(mat) <- cm$gene_id
stopifnot(all(mat == round(mat)), ncol(mat) == 126)
mode(mat) <- "integer"
mat <- mat[!grepl("^MT:", rownames(mat)), ]                       # organellar genes excluded from DE
s <- fread("metadata/samples.tsv"); setkey(s, sample_uid); s <- s[colnames(mat)]
stopifnot(identical(s$sample_uid, colnames(mat)))
s[, batch := factor(batch)]
s[, group := ifelse(is_control, "control", condition)]
wide <- fread("results/crosswalk/gene_id_lookup_wide.tsv")
ann <- merge(ann, wide[, .(gene_id = pn_v5_1_gene_id, symbol, descriptor, vendor_ensembl_r62_ids, vendor_mapping_status,
                            vitvi_v4_ids, vitvi_v3_ids, VIT_v1_ids, max_utr3_bp)], by = "gene_id", all.x = TRUE)
setkey(ann, gene_id)

label_of <- function(cond) {   # physical units
  x <- unique(s[condition == cond, .(hormone, dose_value, dose_unit, time_hours)])
  t <- ifelse(x$time_hours == 168, "7 d", paste0(x$time_hours, " h"))
  sprintf("%s %g %s, %s", x$hormone, x$dose_value, sub("uM", "µM", x$dose_unit), t)
}

prefilter <- function(m, min_count = 10, min_samples = 3) rowSums(m >= min_count) >= min_samples

summ <- list(); focal <- list(); allres <- list()

# ---------------- PRIMARY: within-batch DESeq2 ----------------
for (b in c("1", "2")) {
  sb <- s[batch == b]
  m <- mat[, sb$sample_uid]
  keep <- prefilter(m)
  log_msg("batch", b, ":", ncol(m), "samples;", sum(keep), "genes pass prefilter (>=10 counts in >=3 samples)")
  cd <- data.frame(row.names = sb$sample_uid, group = relevel(factor(sb$group), ref = "control"))
  dds <- DESeqDataSetFromMatrix(m[keep, ], cd, design = ~ group)
  dds <- DESeq(dds, parallel = TRUE)
  saveRDS(dds, file.path(out, sprintf("dds_within_batch%s.rds", b)))
  write.table(data.frame(sample_uid = colnames(dds), size_factor = sizeFactors(dds)),
              file.path(out, sprintf("size_factors_batch%s.tsv", b)), sep = "\t", quote = FALSE, row.names = FALSE)
  pdf(file.path(out, sprintf("dispersion_batch%s.pdf", b))); plotDispEsts(dds); dev.off()
  # edgeR QL cross-check on the same genes
  y <- DGEList(m[keep, ], group = cd$group); y <- calcNormFactors(y)
  des <- model.matrix(~ group, cd); y <- estimateDisp(y, des); fit <- glmQLFit(y, des, robust = TRUE)
  for (cond in setdiff(levels(cd$group), "control")) {
    coef <- paste0("group_", make.names(cond), "_vs_control")
    if (!coef %in% resultsNames(dds)) coef <- resultsNames(dds)[resultsNames(dds) == paste0("group_", cond, "_vs_control")]
    r <- results(dds, name = coef, alpha = FDR, parallel = TRUE)
    sh <- lfcShrink(dds, coef = coef, type = "apeglm", res = r, quiet = TRUE)
    nt <- sum(cd$group == cond); nc <- sum(cd$group == "control")
    tab <- data.table(gene_id = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
                      log2FC_shrunk_apeglm = sh$log2FoldChange, stat_wald = r$stat, pvalue = r$pvalue, padj_BH = r$padj,
                      ci95_low = r$log2FoldChange - 1.96 * r$lfcSE, ci95_high = r$log2FoldChange + 1.96 * r$lfcSE)
    qt <- glmQLFTest(fit, coef = paste0("group", make.names(cond)))
    et <- topTags(qt, n = Inf, sort.by = "none")$table
    tab[, `:=`(edgeR_logFC = et[gene_id, "logFC"], edgeR_F = et[gene_id, "F"], edgeR_PValue = et[gene_id, "PValue"],
               edgeR_FDR = et[gene_id, "FDR"])]
    tab <- merge(ann, tab, by = "gene_id", all.y = TRUE)
    tab[, `:=`(contrast = paste0(cond, "_vs_control_b", b), condition = cond, label = label_of(cond),
               n_treated = nt, n_control = nc, batch = b, model = "DESeq2 ~ group (within batch)",
               direction = "log2(treated/control)",
               control_samples = paste(sb[is_control == TRUE]$sample_uid, collapse = ","))]
    setorder(tab, pvalue, na.last = TRUE)
    fwrite(tab, file.path(out, "tables", paste0(cond, "_vs_control_b", b, ".all_genes.tsv")), sep = "\t")
    sig <- tab[!is.na(padj_BH) & padj_BH < FDR & abs(log2FC) >= LFC_T]
    fwrite(sig, file.path(out, "tables", paste0(cond, "_vs_control_b", b, ".significant.tsv")), sep = "\t")
    summ[[length(summ) + 1]] <- data.table(contrast = tab$contrast[1], label = tab$label[1], condition = cond, batch = b,
      hormone = s[condition == cond]$hormone[1], dose_value = s[condition == cond]$dose_value[1],
      dose_unit = s[condition == cond]$dose_unit[1], time_hours = s[condition == cond]$time_hours[1],
      n_treated = nt, n_control = nc, genes_in_model = nrow(tab), genes_tested = sum(!is.na(tab$padj_BH)),
      n_sig = nrow(sig), n_up = sum(sig$log2FC > 0), n_down = sum(sig$log2FC < 0),
      n_sig_edgeR = sum(tab$edgeR_FDR < FDR, na.rm = TRUE), fdr = FDR, lfc_threshold = LFC_T,
      model = "DESeq2 Wald ~ group (within batch); BH across genes; DESeq2 independent filtering")
    f <- tab[gene_id == FOCAL]
    focal[[length(focal) + 1]] <- if (nrow(f)) f else data.table(gene_id = FOCAL, contrast = tab$contrast[1],
      label = tab$label[1], condition = cond, batch = b, n_treated = nt, n_control = nc, note = "not in model (prefilter)")
    log_msg("  ", tab$contrast[1], ": sig", nrow(sig))
  }
}
summary_dt <- rbindlist(summ, fill = TRUE)
fwrite(summary_dt, file.path(out, "contrast_summary.tsv"), sep = "\t")
invisible(NULL)  # results/de -> pn40024_v5.1/de (symlink)

# ---------------- BRIDGE: batch effect from the 3 control pairs ----------------
sc <- s[is_control == TRUE]
stopifnot(nrow(sc) == 6, all(table(sc$control_label, sc$batch) == 1))
mc <- mat[, sc$sample_uid]; keepc <- prefilter(mc)
cdc <- data.frame(row.names = sc$sample_uid, control_label = factor(sc$control_label), batch = sc$batch)
ddc <- DESeq(DESeqDataSetFromMatrix(mc[keepc, ], cdc, design = ~ control_label + batch), parallel = TRUE)
rc <- results(ddc, name = "batch_2_vs_1", alpha = FDR, parallel = TRUE)
shc <- lfcShrink(ddc, coef = "batch_2_vs_1", type = "apeglm", res = rc, quiet = TRUE)
bt <- data.table(gene_id = rownames(rc), baseMean = rc$baseMean, log2FC_batch2_vs_batch1 = rc$log2FoldChange,
                 lfcSE = rc$lfcSE, log2FC_shrunk_apeglm = shc$log2FoldChange, pvalue = rc$pvalue, padj_BH = rc$padj,
                 model = "DESeq2 ~ control_label + batch (6 control libraries = 3 bridging pairs)")
bt <- merge(ann[, .(gene_id, symbol, descriptor)], bt, by = "gene_id", all.y = TRUE); setorder(bt, pvalue, na.last = TRUE)
fwrite(bt, file.path(out, "bridge_batch_effect", "batch2_vs_batch1_controls.all_genes.tsv"), sep = "\t")
fwrite(data.table(genes_tested = sum(!is.na(bt$padj_BH)), n_sig = sum(bt$padj_BH < FDR, na.rm = TRUE),
                  n_up_in_batch2 = sum(bt$padj_BH < FDR & bt$log2FC_batch2_vs_batch1 > 0, na.rm = TRUE),
                  n_down_in_batch2 = sum(bt$padj_BH < FDR & bt$log2FC_batch2_vs_batch1 < 0, na.rm = TRUE),
                  n_sig_abs_lfc_ge1 = sum(bt$padj_BH < FDR & abs(bt$log2FC_batch2_vs_batch1) >= 1, na.rm = TRUE),
                  median_abs_shrunk_lfc = median(abs(bt$log2FC_shrunk_apeglm), na.rm = TRUE),
                  focal_log2FC = bt[gene_id == FOCAL]$log2FC_batch2_vs_batch1, focal_lfcSE = bt[gene_id == FOCAL]$lfcSE,
                  focal_pvalue = bt[gene_id == FOCAL]$pvalue, focal_padj_BH = bt[gene_id == FOCAL]$padj_BH),
       file.path(out, "bridge_batch_effect", "summary.tsv"), sep = "\t")
log_msg("bridge batch effect: ", sum(bt$padj_BH < FDR, na.rm = TRUE), "genes at BH <", FDR)

# ---------------- JOINT (paired controls): limma-voom + duplicateCorrelation on biological control ----------------
cdj <- data.frame(row.names = s$sample_uid, batch = s$batch, group = relevel(factor(s$group), ref = "control"))
keep <- prefilter(mat)
X <- model.matrix(~ batch + group, cdj); stopifnot(qr(X)$rank == ncol(X))
log_msg("joint model full rank:", ncol(X), "columns")
blk <- ifelse(s$is_control, paste0("ctrl_", s$control_label), s$sample_uid)
y <- calcNormFactors(DGEList(mat[keep, ]))
v <- voom(y, X); dc <- duplicateCorrelation(v, X, block = blk)
v <- voom(y, X, block = blk, correlation = dc$consensus); dc2 <- duplicateCorrelation(v, X, block = blk)
fitb <- eBayes(lmFit(v, X, block = blk, correlation = dc2$consensus))
writeLines(sprintf("consensus within-pair correlation (same control sample, batch 1 vs batch 2 library): %.4f", dc2$consensus),
           file.path(out, "joint_paired", "consensus_correlation.txt"))
jsum <- list(); bf <- list()
for (cond in setdiff(levels(cdj$group), "control")) {
  tt <- topTable(fitb, coef = paste0("group", cond), number = Inf, sort.by = "none")
  tt <- data.table(gene_id = rownames(tt), tt, contrast = paste0(cond, "_vs_control_joint_paired"),
                   model = "limma-voom ~ batch + group; duplicateCorrelation block = biological control sample")
  fwrite(tt, file.path(out, "joint_paired", paste0(cond, "_vs_control_joint_paired.all_genes.tsv")), sep = "\t")
  jsum[[cond]] <- data.table(condition = cond, n_sig_joint_paired = sum(tt$adj.P.Val < FDR))
  bf[[cond]] <- tt[gene_id == FOCAL]
}
js <- rbindlist(jsum)
summary_dt <- merge(summary_dt, js, by = "condition", all.x = TRUE)
fwrite(summary_dt, file.path(out, "contrast_summary.tsv"), sep = "\t")

# ---------------- FOCAL gene family ----------------
ff <- rbindlist(focal, fill = TRUE)
ff[, padj_BH_across_40_focal_tests := p.adjust(pvalue, method = "BH")]
ff[, induction_evidence := fifelse(is.na(pvalue), "not_testable",
        fifelse(log2FC > 0 & padj_BH_across_40_focal_tests < FDR, "increased (FDR_focal<0.05)",
        fifelse(log2FC < 0 & padj_BH_across_40_focal_tests < FDR, "decreased (FDR_focal<0.05)", "no_significant_change")))]
bfd <- rbindlist(lapply(names(bf), function(k) bf[[k]][, .(condition = k, joint_paired_logFC = logFC, joint_paired_P = P.Value)]), fill = TRUE)
ff <- merge(ff, bfd, by = "condition", all.x = TRUE)
ff[, `:=`(joint_paired_padj_BH_across_40 = p.adjust(joint_paired_P, "BH"), edgeR_padj_BH_across_40 = p.adjust(edgeR_PValue, "BH"))]
dir.create("results/pn40024_v5.1/target_gene", showWarnings = FALSE, recursive = TRUE)
fwrite(ff, "results/pn40024_v5.1/target_gene/focal_gene_de_40_contrasts.tsv", sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
log_msg("DE complete")
