#!/usr/bin/env Rscript
# Embryogenic signature from the SRA atlas, scored in the PPPJ8P hormone experiment.
#
# Question: in our callus, does the focal gene (Vitvi05_01chr01g05750) move together with a
# general embryogenic state, or on its own?
#
# 1. Signature (atlas Salmon counts, first 10 M reads, PN40024 v5.1; gene = sum of transcripts):
#    S1 PRJNA552475 embryogenic vs non-embryogenic callus, ~ cultivar + state (6 vs 6)
#    S2 PRJNA980588 early somatic embryo (early induction, PG, globular 1-2) vs late
#       (cotyledonary 1-2, seedling, epicotyl, juvenile plant); heart/torpedo left out
#    signature_up   = BH < 0.01 and log2FC > 1 in BOTH; signature_down = same, < -1 in BOTH.
#    The focal gene is excluded from the signature (it is the thing being explained).
# 2. Contrast level (40 DESeq2 contrasts vs control, de/tables): score change =
#    mean shrunk log2FC(up genes) - mean shrunk log2FC(down genes); correlated with the focal
#    gene's shrunk log2FC across contrasts.
# 3. Sample level (technically adjusted, batch-removed matrix from focal_coexpression.R; unadjusted vst for comparison):
#    score = mean z(up) - mean z(down) (also up-only / down-only); correlation with focal
#    vst between group means and within groups (residuals after removing group means).
# Null for 2-3: 1000 random gene sets of the same sizes, matched to the signature genes'
# expression (baseMean deciles in our data), focal excluded.
# Outputs: results/sra_atlas/embryogenic_signature/
suppressPackageStartupMessages({ library(DESeq2); library(data.table) })
set.seed(20261107)
ROOT  <- Sys.getenv("VVECPG1_UPSTREAM_ROOT", normalizePath("."))
OUT   <- file.path(ROOT, "results/sra_atlas/embryogenic_signature"); dir.create(OUT, showWarnings = FALSE)
FOCAL <- "Vitvi05_01chr01g05750"
NPERM <- 1000

# ---------- 1. atlas DESeq2 ----------
runs <- fread(file.path(ROOT, "results/sra_atlas/focal_atlas_per_run.tsv"))[qc_ok == TRUE]
tx2g <- fread(file.path(ROOT, "references/salmon_v51/tx2gene.tsv"), header = FALSE, col.names = c("tx", "gene"))
read_counts <- function(ids) {
  m <- sapply(ids, function(r) {
    q <- fread(cmd = paste("gzip -dc", shQuote(file.path(ROOT, "results/sra_atlas/quant", r, "quant.sf.gz"))), select = c("Name", "NumReads"))
    q[, gene := tx2g$gene[match(Name, tx2g$tx)]]
    g <- q[, .(n = sum(NumReads)), by = gene]
    setNames(g$n, g$gene)
  })
  round(m)
}
de <- function(cnt, cd, design, coef) {
  keep <- rowSums(cnt >= 10) >= 3
  dds <- DESeqDataSetFromMatrix(cnt[keep, ], cd, design)
  dds <- DESeq(dds, quiet = TRUE)
  r <- results(dds, name = coef, alpha = 0.01)
  data.table(gene_id = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange, padj = r$padj)
}

a <- runs[study_accession == "PRJNA552475" & tissue_class %in% c("embryogenic callus", "non-embryogenic callus")]
a[, cultivar2 := substr(sample_title, 1, 2)]
a[, state := factor(ifelse(tissue_class == "embryogenic callus", "EC", "NEC"), levels = c("NEC", "EC"))]
s1 <- de(read_counts(a$run_accession), data.frame(cv = factor(a$cultivar2), state = a$state),
         ~ cv + state, "state_EC_vs_NEC")

b <- runs[study_accession == "PRJNA980588"]
b[, stage := sub("_rep\\d+$", "", sample_title)]
b <- b[stage %in% c("EI", "PG", "G1", "G2", "C1", "C2", "S", "EP", "JP")]
b[, phase := factor(ifelse(stage %in% c("EI", "PG", "G1", "G2"), "early", "late"), levels = c("late", "early"))]
s2 <- de(read_counts(b$run_accession), data.frame(phase = b$phase), ~ phase, "phase_early_vs_late")

sig <- merge(s1, s2, by = "gene_id", suffixes = c("_EC_vs_NEC", "_SEearly_vs_late"))
sig[, class := fifelse(padj_EC_vs_NEC < 0.01 & padj_SEearly_vs_late < 0.01 &
                         log2FC_EC_vs_NEC > 1 & log2FC_SEearly_vs_late > 1, "up",
               fifelse(padj_EC_vs_NEC < 0.01 & padj_SEearly_vs_late < 0.01 &
                         log2FC_EC_vs_NEC < -1 & log2FC_SEearly_vs_late < -1, "down", "none"))]
focal_atlas <- sig[gene_id == FOCAL]
sig[gene_id == FOCAL, class := "focal_excluded"]
ann <- fread(cmd = paste("gzip -dc", shQuote(file.path(ROOT, "results/pn40024_v5.1/focal_network/focal_coexpression_all_genes.tsv.gz"))),
             select = c("gene_id", "bicor_C", "robust", "direction", "v51_symbol", "v51_descriptor", "at_symbol", "at_description"))
sig <- merge(sig, ann, by = "gene_id", all.x = TRUE)
fwrite(sig[order(class, -log2FC_EC_vs_NEC)], file.path(OUT, "atlas_de_EC_vs_NEC_and_SE_early_vs_late.tsv.gz"), sep = "\t")
cat("atlas genes tested in both:", nrow(sig), "\n"); print(table(sig$class))
cat("focal gene in atlas DE:\n"); print(focal_atlas)
cat("agreement of the two atlas contrasts (genes sig in both): r =",
    round(sig[padj_EC_vs_NEC < 0.01 & padj_SEearly_vs_late < 0.01, cor(log2FC_EC_vs_NEC, log2FC_SEearly_vs_late)], 2), "\n")

# ---------- our data ----------
cs   <- fread(file.path(ROOT, "results/pn40024_v5.1/de/contrast_summary.tsv"))
tabs <- lapply(cs$contrast, function(cn) fread(file.path(ROOT, "results/pn40024_v5.1/de/tables", paste0(cn, ".all_genes.tsv")),
                                                select = c("gene_id", "baseMean", "log2FC_shrunk_apeglm")))
lfc  <- Reduce(function(x, y) merge(x, y, by = "gene_id"),
               Map(function(t, cn) setnames(t[, .(gene_id, log2FC_shrunk_apeglm)], 2, cn), tabs, cs$contrast))
lfcm <- as.matrix(lfc[, -1]); rownames(lfcm) <- lfc$gene_id; lfcm[is.na(lfcm)] <- 0
bm   <- tabs[[1]][match(rownames(lfcm), gene_id), baseMean]
names(bm) <- rownames(lfcm)

up   <- intersect(sig[class == "up", gene_id], rownames(lfcm))
down <- intersect(sig[class == "down", gene_id], rownames(lfcm))
cat("signature genes present in our data: up", length(up), "down", length(down), "\n")
sig[, in_our_data := gene_id %in% rownames(lfcm)]
fwrite(sig[class %in% c("up", "down")][order(class, -abs(log2FC_EC_vs_NEC))], file.path(OUT, "signature_genes.tsv"), sep = "\t")

# matched random sets: same size, same baseMean-decile composition
dec <- cut(rank(log10(bm + 1)), 10, labels = FALSE); names(dec) <- names(bm)
pool <- split(setdiff(names(bm), FOCAL), dec[setdiff(names(bm), FOCAL)])
draw <- function(g) unlist(lapply(split(g, dec[g]), function(x) sample(pool[[dec[x[1]]]], length(x))))
score_lfc <- function(u, d) colMeans(lfcm[u, , drop = FALSE]) - colMeans(lfcm[d, , drop = FALSE])


# Three versions of the score throughout: up - down ("both"), up genes only, down genes only (sign
# flipped so higher = more embryogenic). The two atlas contrasts agree on direction for the up set far
# more cleanly than for the down set (late SE samples include green plantlets), so "up" is the primary score.
versions <- list(both = c(TRUE, TRUE), up = c(TRUE, FALSE), down = c(FALSE, TRUE))
mkscore <- function(M, u, d, v) {
  s <- 0
  if (v[1]) s <- s + colMeans(M[intersect(u, rownames(M)), , drop = FALSE])
  if (v[2]) s <- s - colMeans(M[intersect(d, rownames(M)), , drop = FALSE])
  s
}
nullsum <- function(obs, null) c(obs = obs, null_median = median(null), null_q025 = unname(quantile(null, 0.025)),
                                 null_q975 = unname(quantile(null, 0.975)), p_upper = (1 + sum(null >= obs)) / (length(null) + 1))
res <- list()
rand_sets <- replicate(NPERM, list(u = draw(up), d = draw(down)), simplify = FALSE)

# ---------- 2. contrast level ----------
cs[, focal_lfc := lfcm[FOCAL, contrast]]
for (vn in names(versions)) {
  v <- versions[[vn]]
  sc <- mkscore(lfcm, up, down, v)[cs$contrast]
  cs[, (paste0("score_", vn)) := sc]
  nsc <- sapply(rand_sets, function(rs) mkscore(lfcm, rs$u, rs$d, v)[cs$contrast])
  # per-contrast: is the signature shift larger than random sets? two-sided, BH across the 40 contrasts
  p <- sapply(seq_along(sc), function(i) (1 + sum(abs(nsc[i, ] - median(nsc[i, ])) >= abs(sc[i] - median(nsc[i, ])))) / (NPERM + 1))
  cs[, (paste0("score_", vn, "_null_median")) := apply(nsc, 1, median)]
  cs[, (paste0("score_", vn, "_padj")) := p.adjust(p, "BH")]
  res[[paste0("contrast_r_", vn)]] <- nullsum(cor(cs$focal_lfc, sc), apply(nsc, 2, function(x) cor(cs$focal_lfc, x)))
}
fwrite(cs[, c("contrast", "label", "batch", "hormone", "dose_value", "dose_unit", "time_hours", "focal_lfc",
              grep("^score_", names(cs), value = TRUE)), with = FALSE], file.path(OUT, "contrast_scores.tsv"), sep = "\t")

# ---------- 3. sample level ----------
# Primary: the batch-removed, technically adjusted matrix from the focal co-expression analysis
# (focal_coexpression.R context C: removeBatchEffect via bridging design + 3 within-group technical PCs,
# treatment groups protected). Secondary (for comparison): blind vst, z-scored within batch, unadjusted.
cc  <- readRDS(file.path(ROOT, "results/pn40024_v5.1/focal_network/context_C_adjusted.rds"))
smp <- cc$s
mats <- list(adjusted = t(scale(t(cc$adj))))
vst <- fread(file.path(ROOT, "results/pn40024_v5.1/counts/pn_v51.vst_blind_visualization_only.tsv"))
vm  <- as.matrix(vst[, -1]); rownames(vm) <- vst$gene_id; vm <- vm[, smp$sample_uid]
for (bt in unique(smp$batch)) { j <- smp$batch == bt; vm[, j] <- t(scale(t(vm[, j]))) }
vm[!is.finite(vm)] <- 0
mats$unadjusted <- vm
grp <- paste(smp$batch, smp$group)
resid <- function(x) x - ave(x, grp)
for (mn in names(mats)) {
  M <- mats[[mn]]; M[!is.finite(M)] <- 0
  fz <- M[FOCAL, ]
  within_r  <- function(s) cor(resid(fz), resid(s))
  between_r <- function(s) { m <- tapply(s, grp, mean); f <- tapply(fz, grp, mean); cor(f, m[names(f)]) }
  smp[, (paste0("focal_z_", mn)) := fz]
  for (vn in names(versions)) {
    v <- versions[[vn]]
    s <- mkscore(M, up, down, v)
    smp[, (paste0("score_", vn, "_", mn)) := s]
    ns <- sapply(rand_sets, function(rs) { x <- mkscore(M, rs$u, rs$d, v); c(within_r(x), between_r(x)) })
    res[[paste0("sample_within_r_", vn, "_", mn)]]  <- nullsum(within_r(s), ns[1, ])
    res[[paste0("sample_between_r_", vn, "_", mn)]] <- nullsum(between_r(s), ns[2, ])
  }
}
fwrite(smp[, c("sample_uid", "batch", "group", "hormone", "dose_value", "dose_unit", "time_hours", "is_control",
               grep("^(focal_z|score)_", names(smp), value = TRUE)), with = FALSE], file.path(OUT, "sample_scores.tsv"), sep = "\t")

# ---------- marker genes ----------
markers <- c(AGL15 = "Vitvi05_01chr13g24750", SERK1 = "Vitvi05_01chr18g21040", FUS3 = "Vitvi05_01chr14g22320",
             WOX2 = "Vitvi05_01chr06g06770", BBM = "Vitvi05_01chr04g17020", PLT2 = "Vitvi05_01chr06g02840",
             ABI3 = "Vitvi05_01chr07g12060")
mk <- data.table(marker = names(markers), gene_id = markers)
mk <- merge(mk, sig[, .(gene_id, log2FC_EC_vs_NEC, padj_EC_vs_NEC, log2FC_SEearly_vs_late, padj_SEearly_vs_late, class, bicor_C)],
            by = "gene_id", all.x = TRUE)
A <- mats$adjusted
mk[, r_with_focal_contrasts := sapply(gene_id, function(g) if (g %in% rownames(lfcm)) cor(lfcm[g, cs$contrast], cs$focal_lfc) else NA)]
mk[, r_with_focal_within_group_adjusted := sapply(gene_id, function(g) if (g %in% rownames(A)) cor(resid(A[g, ]), resid(A[FOCAL, ])) else NA)]
fwrite(mk, file.path(OUT, "marker_genes.tsv"), sep = "\t"); print(mk[, -1])

summ <- rbindlist(lapply(names(res), function(k) as.data.table(as.list(c(test = k, round(res[[k]], 4))))))
summ <- rbind(data.table(test = c("atlas_up", "atlas_down", "up_in_ours", "down_in_ours"),
                         obs = c(sig[class == "up", .N], sig[class == "down", .N], length(up), length(down))),
              summ, fill = TRUE)
fwrite(summ, file.path(OUT, "summary.tsv"), sep = "\t"); print(summ)
print(cs[time_hours == 24, .(label, focal_lfc = round(focal_lfc, 2), up = round(score_up, 2), up_null = round(score_up_null_median, 2),
                             up_padj = signif(score_up_padj, 2), both = round(score_both, 2), both_padj = signif(score_both_padj, 2))])
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
