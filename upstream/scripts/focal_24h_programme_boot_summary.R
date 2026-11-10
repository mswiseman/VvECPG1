#!/usr/bin/env Rscript
# Collect bootstrap refits (scripts/focal_24h_programme_boot.R) -> percentile summaries.
# Outputs: focal_24h_bootstrap_summary.tsv, focal_24h_bootstrap_contrasts.tsv, focal_24h_bootstrap_replicates.tsv
suppressPackageStartupMessages(library(data.table))
D <- file.path(Sys.getenv("VVECPG1_UPSTREAM_ROOT", normalizePath(".")), "results/pn40024_v5.1/focal_24h_programme")
ct <- fread(file.path(D, "focal_24h_contrasts.tsv"))
key <- sub("_vs_control_b1$", "", ct$contrast)

# pipeline check: task 0 = refit without resampling must reproduce the saved tables
chk <- readRDS(file.path(D, "bootstrap/boot_task000.rds"))[[1]]
d_f <- max(abs(chk$focal[key] - ct$focal_log2FC_shrunk_apeglm)); d_s <- max(abs(chk$shift[key] - ct$programme_shift))
cat(sprintf("refit check: max |focal LFC diff| %.2e, max |shift diff| %.2e\n", d_f, d_s))

B <- unlist(lapply(setdiff(list.files(file.path(D, "bootstrap"), "^boot_task\\d+\\.rds$", full.names = TRUE),
                           file.path(D, "bootstrap/boot_task000.rds")), readRDS), recursive = FALSE)
nb <- length(B)
rep_t <- rbindlist(lapply(seq_len(nb), function(i) data.table(boot = i, pearson_r_shift = B[[i]]$r_shift,
                                                              spearman_rho_shift = B[[i]]$rho_shift, pearson_r_raw_score = B[[i]]$r_raw,
                                                              control_draw = B[[i]]$draw[["control"]],
                                                              n_distinct_controls = length(unique(strsplit(B[[i]]$draw[["control"]], ",")[[1]])))))
fwrite(rep_t, file.path(D, "focal_24h_bootstrap_replicates.tsv"), sep = "\t")
obs <- c(pearson_r_shift = cor(ct$focal_log2FC_shrunk_apeglm, ct$programme_shift),
         spearman_rho_shift = cor(ct$focal_log2FC_shrunk_apeglm, ct$programme_shift, method = "spearman"),
         pearson_r_raw_score = cor(ct$focal_log2FC_shrunk_apeglm, ct$programme_score))
sm <- rbindlist(lapply(names(obs), function(k) {
  x <- rep_t[[k]]
  data.table(metric = k, observed = obs[[k]], n_boot = nb, boot_median = median(x), q025 = quantile(x, .025), q975 = quantile(x, .975),
             q10 = quantile(x, .10), q90 = quantile(x, .90), frac_gt0 = mean(x > 0))
}))
sm[, refit_check_max_abs_diff_focal := d_f][, refit_check_max_abs_diff_shift := d_s]
fwrite(sm, file.path(D, "focal_24h_bootstrap_summary.tsv"), sep = "\t")
F  <- sapply(B, function(b) b$focal[key]); SH <- sapply(B, function(b) b$shift[key])
bc <- data.table(contrast = ct$contrast, label = ct$label,
                 focal_q025 = apply(F, 1, quantile, .025), focal_q975 = apply(F, 1, quantile, .975),
                 shift_q025 = apply(SH, 1, quantile, .025), shift_q975 = apply(SH, 1, quantile, .975))
fwrite(bc, file.path(D, "focal_24h_bootstrap_contrasts.tsv"), sep = "\t")
print(sm); print(bc)
cat("replicates with only 1 distinct control:", sum(rep_t$n_distinct_controls == 1), "of", nb, "\n")
