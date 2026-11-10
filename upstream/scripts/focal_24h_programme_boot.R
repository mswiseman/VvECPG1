#!/usr/bin/env Rscript
# Biological-sample bootstrap for the 24 h focal-vs-programme association (scripts/focal_24h_programme.R).
# Each replicate refits the ORIGINAL batch-1 DESeq2 model (~ group, all 96 batch-1 libraries, same prefiltered genes
# as results/pn40024_v5.1/de/dds_within_batch1.rds) after resampling:
#   - the 3 control libraries ONCE (with replacement), reused by all ten 24 h contrasts (they share the controls);
#   - the 3 treated libraries of each 24 h group, with replacement, within group;
#   - the other 21 batch-1 groups unchanged (they only inform the shared dispersion trend, as in the original).
# Then apeglm shrinkage per 24 h coefficient (as de_analysis.R), programme score/shift with the SAME 1000 random sets
# (signature_and_random_sets.rds), and correlations across the ten contrasts.
# Usage: Rscript focal_24h_programme_boot.R <task_id> <n_boot>   (task 0 = no resampling, pipeline check)
suppressPackageStartupMessages({ library(DESeq2); library(apeglm); library(data.table) })
args <- commandArgs(trailingOnly = TRUE); task <- as.integer(args[1]); nb <- as.integer(args[2])
ROOT <- Sys.getenv("VVECPG1_UPSTREAM_ROOT", normalizePath("."))
OUT  <- file.path(ROOT, "results/pn40024_v5.1/focal_24h_programme/bootstrap"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FOCAL <- "Vitvi05_01chr01g05750"
rs <- readRDS(file.path(ROOT, "results/pn40024_v5.1/focal_24h_programme/signature_and_random_sets.rds"))
dds0 <- readRDS(file.path(ROOT, "results/pn40024_v5.1/de/dds_within_batch1.rds"))
cnt0 <- counts(dds0); grp0 <- as.character(colData(dds0)$group)
G24 <- sort(unique(grp0[grepl("_24h$", grp0)])); stopifnot(length(G24) == 10)
set.seed(20261108 + task)

one <- function(resample) {
  cols <- split(colnames(cnt0), grp0)
  if (resample) for (g in c("control", G24)) cols[[g]] <- sample(cols[[g]], length(cols[[g]]), replace = TRUE)
  idx <- unlist(cols, use.names = FALSE); grp <- rep(names(cols), lengths(cols))
  m <- cnt0[, idx]; colnames(m) <- paste0(idx, "_b", seq_along(idx))
  cd <- data.frame(row.names = colnames(m), group = relevel(factor(grp), ref = "control"))
  dds <- DESeq(DESeqDataSetFromMatrix(m, cd, ~ group), quiet = TRUE)
  L <- sapply(G24, function(g) {
    coef <- paste0("group_", make.names(g), "_vs_control")
    lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)$log2FoldChange
  })
  rownames(L) <- rownames(dds); L[is.na(L)] <- 0
  sc   <- colMeans(L[rs$up, ])
  nmed <- apply(sapply(rs$rand_up, function(u) colMeans(L[u, , drop = FALSE])), 1, median)
  f <- L[FOCAL, ]; sh <- sc - nmed
  list(draw = vapply(cols[c("control", G24)], paste, "", collapse = ","),
       focal = f, score = sc, null_med = nmed, shift = sh,
       r_shift = cor(f, sh), rho_shift = cor(f, sh, method = "spearman"), r_raw = cor(f, sc))
}
t0 <- Sys.time()
res <- if (task == 0) list(one(FALSE)) else replicate(nb, one(TRUE), simplify = FALSE)
saveRDS(res, file.path(OUT, sprintf("boot_task%03d.rds", task)))
cat("task", task, "done:", length(res), "replicates in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
print(sapply(res, function(x) c(r_shift = x$r_shift, rho = x$rho_shift, r_raw = x$r_raw)))
