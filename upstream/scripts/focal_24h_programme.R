#!/usr/bin/env Rscript
# Follow-up (requested 2026-11-08): does the focal gene (VvECPG1 = Vitvi05_01chr01g05750 = Vitvi01g01965 =
# VIT_01s0011g04570) track the atlas embryogenic/proliferative signature SPECIFICALLY at 24 h, the only time with
# time-matched untreated controls?
#
# Reuses, unchanged, the inputs and calculations of scripts/embryogenic_signature.R (results there are not touched):
#   - signature: results/sra_atlas/embryogenic_signature/signature_genes.tsv (class "up", 964 atlas genes; focal excluded)
#   - per-contrast apeglm-shrunken log2FC: results/pn40024_v5.1/de/tables/<contrast>.all_genes.tsv (within-batch DESeq2
#     ~ group; NA -> 0, as in the original); gene universe = genes in all 40 tables
#   - random sets: 1000 sets matched to the signature on baseMean decile, focal excluded, drawn with the original seed
#     and draw order (verified below against the saved per-contrast null medians)
#   - programme score = mean shrunken log2FC of the up genes; programme shift = score - median score of the 1000 sets
#   - sample level: context-C technically adjusted matrix (focal_coexpression.R) z-scored per gene, and the blind VST
#     z-scored within batch (unadjusted); within-group association = Pearson r of residuals from (batch x group) means
# Outputs: results/pn40024_v5.1/focal_24h_programme/
suppressPackageStartupMessages({ library(data.table) })
set.seed(20261107)                       # same seed as embryogenic_signature.R
ROOT  <- Sys.getenv("VVECPG1_UPSTREAM_ROOT", normalizePath("."))
SIGD  <- file.path(ROOT, "results/sra_atlas/embryogenic_signature")
OUT   <- file.path(ROOT, "results/pn40024_v5.1/focal_24h_programme"); dir.create(OUT, showWarnings = FALSE)
FOCAL <- "Vitvi05_01chr01g05750"
NPERM <- 1000

# ---------- signature ----------
sig <- fread(file.path(SIGD, "signature_genes.tsv"))
stopifnot(!FOCAL %in% sig$gene_id)
sig_all <- fread(cmd = paste("gzip -dc", shQuote(file.path(SIGD, "atlas_de_EC_vs_NEC_and_SE_early_vs_late.tsv.gz"))), select = c("gene_id", "class"))
cat("focal gene class in atlas DE table:", sig_all[gene_id == FOCAL, class], "\n")

# ---------- shrunken log2FC, as in embryogenic_signature.R ----------
cs   <- fread(file.path(ROOT, "results/pn40024_v5.1/de/contrast_summary.tsv"))
tabs <- lapply(cs$contrast, function(cn) fread(file.path(ROOT, "results/pn40024_v5.1/de/tables", paste0(cn, ".all_genes.tsv")),
                                                select = c("gene_id", "baseMean", "log2FC_shrunk_apeglm")))
lfc  <- Reduce(function(x, y) merge(x, y, by = "gene_id"),
               Map(function(t, cn) setnames(t[, .(gene_id, log2FC_shrunk_apeglm)], 2, cn), tabs, cs$contrast))
lfcm <- as.matrix(lfc[, -1]); rownames(lfcm) <- lfc$gene_id
n_na_sig <- sum(is.na(lfcm[intersect(sig[class == "up", gene_id], rownames(lfcm)), ]))
lfcm[is.na(lfcm)] <- 0
bm   <- tabs[[1]][match(rownames(lfcm), gene_id), baseMean]; names(bm) <- rownames(lfcm)
# original order of the signature vectors (merge in the original sorted by gene_id)
up   <- intersect(sort(sig[class == "up", gene_id]), rownames(lfcm))
down <- intersect(sort(sig[class == "down", gene_id]), rownames(lfcm))
cat("signature up genes:", sig[class == "up", .N], "atlas;", length(up), "in our 40-contrast universe;",
    n_na_sig, "NA shrunken LFC values among them set to 0 (as original)\n")

# ---------- random sets: identical construction and draw order ----------
dec  <- cut(rank(log10(bm + 1)), 10, labels = FALSE); names(dec) <- names(bm)
pool <- split(setdiff(names(bm), FOCAL), dec[setdiff(names(bm), FOCAL)])
draw <- function(g) unlist(lapply(split(g, dec[g]), function(x) sample(pool[[dec[x[1]]]], length(x))))
rand_sets <- replicate(NPERM, list(u = draw(up), d = draw(down)), simplify = FALSE)
saveRDS(list(up = up, rand_up = lapply(rand_sets, `[[`, "u")), file.path(OUT, "signature_and_random_sets.rds"))

score <- function(M, g) colMeans(M[g, , drop = FALSE])
S_obs  <- score(lfcm, up)
S_null <- sapply(rand_sets, function(rs) score(lfcm, rs$u))          # contrasts x 1000
null_med <- apply(S_null, 1, median)
shift <- S_obs - null_med
focal <- lfcm[FOCAL, ]

# check against the saved original
orig <- fread(file.path(SIGD, "contrast_scores.tsv"))
chk <- data.table(contrast = names(S_obs), score = S_obs, null_med = null_med)[orig, on = "contrast"]
cat("max |score - original|:", signif(max(abs(chk$score - chk$score_up)), 3),
    "; max |null median - original|:", signif(max(abs(chk$null_med - chk$score_up_null_median)), 3), "\n")
stopifnot(max(abs(chk$score - chk$score_up)) < 1e-9)
null_reproduced <- max(abs(chk$null_med - chk$score_up_null_median)) < 1e-9

# ---------- helpers ----------
nullstats <- function(obs, null) list(obs = obs, null_median = median(null), null_q025 = unname(quantile(null, .025)),
                                      null_q975 = unname(quantile(null, .975)), p_upper = (1 + sum(null >= obs)) / (length(null) + 1),
                                      null_percentile = mean(null < obs))
summ <- list()
add <- function(analysis, metric, n_contrasts = NA, n_libraries = NA, n_biological = NA, n_groups = NA, resid_df = NA,
                obs, null = NULL, note = "") {
  r <- data.table(analysis = analysis, metric = metric, n_contrasts = n_contrasts, n_libraries = n_libraries,
                  n_biological_samples = n_biological, n_groups = n_groups, residual_df = resid_df, value = obs, note = note)
  if (!is.null(null)) r <- cbind(r, as.data.table(nullstats(obs, null)[-1]))
  summ[[length(summ) + 1]] <<- r
}

# ---------- 1a. all 40 contrasts: reproduce ----------
cor_null <- function(x, idx, M) apply(M[idx, , drop = FALSE], 2, function(s) cor(x[idx], s))
i40 <- seq_along(focal)
add("all 40 contrasts (reproduction)", "Pearson r, focal shrunken LFC vs raw programme score (mean shrunken LFC of up genes)", 40,
    obs = cor(focal, S_obs), null = cor_null(focal, i40, S_null),
    note = "the value reported as r = 0.48 in executive summary 8.12 / summary.tsv contrast_r_up")
add("all 40 contrasts (reproduction)", "Pearson r, focal shrunken LFC vs programme shift (score - random-set median)", 40,
    obs = cor(focal, shift), null = cor_null(focal, i40, S_null - null_med),
    note = "the quantity plotted on the y axis of embryogenic_signature figure panel B")
add("all 40 contrasts (reproduction)", "Spearman rho, focal vs programme shift", 40, obs = cor(focal, shift, method = "spearman"))

# ---------- 1b. the ten 24 h contrasts ----------
i24 <- which(cs$time_hours == 24)
stopifnot(length(i24) == 10, all(cs$batch[i24] == 1), all(grepl("_vs_control_b1$", cs$contrast[i24])))
f24 <- focal[i24]; sh24 <- shift[i24]; sc24 <- S_obs[i24]
add("24 h contrasts", "Pearson r, focal shrunken LFC vs programme shift", 10, obs = cor(f24, sh24),
    null = cor_null(focal, i24, S_null - null_med),
    note = "primary. Null: same r with each random set's shift (its score minus the per-contrast random-set median), same 10 contrasts")
add("24 h contrasts", "Spearman rho, focal shrunken LFC vs programme shift", 10, obs = cor(f24, sh24, method = "spearman"),
    null = apply(S_null[i24, ] - null_med[i24], 2, function(s) cor(f24, s, method = "spearman")))
add("24 h contrasts", "Pearson r, focal shrunken LFC vs raw programme score", 10, obs = cor(f24, sc24),
    null = cor_null(focal, i24, S_null),
    note = "same definition as the original all-40 r; raw score includes the generic shift shared with random gene sets")
add("24 h contrasts", "Pearson r, focal shrunken LFC vs random-set median score", 10, obs = cor(f24, null_med[i24]),
    note = "how much the focal LFC tracks the generic (non-signature) component of the response")
add("24 h contrasts", "range of programme shift (log2)", 10, obs = diff(range(sh24)), note = sprintf("%.3f to %.3f", min(sh24), max(sh24)))
add("24 h contrasts", "range of focal shrunken LFC (log2)", 10, obs = diff(range(f24)), note = sprintf("%.3f to %.3f", min(f24), max(f24)))
# 6 h and 12 h for context (not time-matched)
for (tt in c(6, 12, 168)) {
  ii <- which(cs$time_hours == tt)
  add(sprintf("%s contrasts (context; not time-matched)", ifelse(tt == 168, "7 d", paste(tt, "h"))), "Pearson r, focal shrunken LFC vs programme shift",
      length(ii), obs = cor(focal[ii], shift[ii]), null = cor_null(focal, ii, S_null - null_med))
}

# leave one contrast out
loo <- rbindlist(lapply(seq_along(i24), function(k) {
  j <- i24[-k]
  data.table(left_out = cs$label[i24[k]], contrast = cs$contrast[i24[k]],
             pearson_r_shift = cor(focal[j], shift[j]), spearman_rho_shift = cor(focal[j], shift[j], method = "spearman"),
             pearson_r_raw_score = cor(focal[j], S_obs[j]),
             null_p_upper_shift = { nn <- cor_null(focal, j, S_null - null_med); (1 + sum(nn >= cor(focal[j], shift[j]))) / (NPERM + 1) })
}))
loo[, change_in_r := pearson_r_shift - cor(f24, sh24)]
fwrite(loo, file.path(OUT, "focal_24h_leave_one_out.tsv"), sep = "\t")
add("24 h contrasts, leave one out", "minimum Pearson r (shift)", 9, obs = min(loo$pearson_r_shift), note = loo[which.min(pearson_r_shift), left_out])
add("24 h contrasts, leave one out", "maximum Pearson r (shift)", 9, obs = max(loo$pearson_r_shift), note = loo[which.max(pearson_r_shift), left_out])

# contrast table
ft <- fread(file.path(ROOT, "results/pn40024_v5.1/target_gene/focal_gene_de_40_contrasts.tsv"))
ct <- data.table(contrast = cs$contrast[i24], label = cs$label[i24], hormone = cs$hormone[i24], concentration = cs$dose_value[i24],
                 unit = sub("uM", "µM", cs$dose_unit[i24]), time_hours = 24L, batch = cs$batch[i24],
                 n_treated = cs$n_treated[i24], n_control = cs$n_control[i24],
                 focal_log2FC_shrunk_apeglm = f24, programme_score = sc24, random_set_median_score = null_med[i24],
                 random_set_score_q025 = apply(S_null[i24, ], 1, quantile, .025), random_set_score_q975 = apply(S_null[i24, ], 1, quantile, .975),
                 programme_shift = sh24)
ct <- merge(ct, ft[, .(contrast, focal_baseMean = baseMean, focal_log2FC_MLE = log2FC, focal_lfcSE = lfcSE, focal_ci95_low = ci95_low,
                       focal_ci95_high = ci95_high, focal_wald_stat = stat_wald, focal_pvalue = pvalue,
                       focal_padj_BH_genomewide = padj_BH, focal_padj_BH_across_40_focal_tests = padj_BH_across_40_focal_tests,
                       focal_edgeR_logFC = edgeR_logFC, focal_edgeR_PValue = edgeR_PValue)], by = "contrast", sort = FALSE)
ct <- merge(ct, orig[, .(contrast, programme_shift_padj_BH_across_40_original = score_up_padj)], by = "contrast", sort = FALSE)
p10 <- sapply(i24, function(i) (1 + sum(abs(S_null[i, ] - null_med[i]) >= abs(shift[i]))) / (NPERM + 1))
ct[, programme_shift_p_vs_random_sets := p10][, programme_shift_padj_BH_across_10 := p.adjust(p10, "BH")]
fwrite(ct, file.path(OUT, "focal_24h_contrasts.tsv"), sep = "\t")

# ---------- 3. within-group association (sample level) ----------
cc  <- readRDS(file.path(ROOT, "results/pn40024_v5.1/focal_network/context_C_adjusted.rds"))
smp <- cc$s
mats <- list(adjusted = t(scale(t(cc$adj))))
vst <- fread(file.path(ROOT, "results/pn40024_v5.1/counts/pn_v51.vst_blind_visualization_only.tsv"))
vm  <- as.matrix(vst[, -1]); rownames(vm) <- vst$gene_id; vm <- vm[, smp$sample_uid]
for (bt in unique(smp$batch)) { j <- smp$batch == bt; vm[, j] <- t(scale(t(vm[, j]))) }
vm[!is.finite(vm)] <- 0
mats$unadjusted <- vm
grp_all <- paste(smp$batch, smp$group)
subsets <- list(
  "all libraries (original)" = rep(TRUE, nrow(smp)),
  "treated only (controls excluded)" = !smp$is_control,
  "24 h treated groups only" = !smp$is_control & smp$time_hours == 24)
resid_rows <- list()
for (mn in names(mats)) {
  M <- mats[[mn]]; M[!is.finite(M)] <- 0
  gup <- intersect(up, rownames(M))
  s_obs <- colMeans(M[gup, , drop = FALSE])
  s_nul <- sapply(rand_sets, function(rs) colMeans(M[intersect(rs$u, rownames(M)), , drop = FALSE]))
  fz <- M[FOCAL, ]
  for (sn in names(subsets)) {
    k <- subsets[[sn]]; g <- grp_all[k]
    rs <- function(x) x - ave(x, g)
    fr <- rs(fz[k]); sr <- rs(s_obs[k])
    nul <- apply(s_nul[k, , drop = FALSE], 2, function(x) cor(fr, rs(x)))
    nlib <- sum(k); nbio <- nlib - sum(duplicated(smp$control_label[k][smp$is_control[k]]))
    add(paste0("within-group, ", mn), paste0("Pearson r of residuals from group means: ", sn), NA, nlib, nbio,
        length(unique(g)), nlib - length(unique(g)), obs = cor(fr, sr), null = nul,
        note = sprintf("%d signature genes in matrix; pooled residual correlation (not an average of per-group r)", length(gup)))
    if (mn == "adjusted" || sn != "all libraries (original)")
      resid_rows[[length(resid_rows) + 1]] <- data.table(matrix = mn, subset = sn, sample_uid = smp$sample_uid[k], group = g,
                                                         hormone = smp$hormone[k], time_hours = smp$time_hours[k],
                                                         focal_resid = fr, programme_resid = sr)
  }
  # sensitivity (treated-only and 24 h subsets): Spearman of residuals, MeJA groups removed, the three low-connectivity
  # MeJA libraries (PPPJ8P_19, _52, _81; WGCNA, kept in all primary analyses) removed, leave one group out
  for (sn in names(subsets)[2:3]) {
    base <- subsets[[sn]]
    rr <- function(k, method = "pearson") { g <- grp_all[k]; rs <- function(x) x - ave(x, g); cor(rs(fz[k]), rs(s_obs[k]), method = method) }
    an <- paste0("within-group sensitivity, ", mn)
    k <- base; add(an, paste0("Spearman rho of residuals: ", sn), NA, sum(k), sum(k), length(unique(grp_all[k])), sum(k) - length(unique(grp_all[k])), obs = rr(k, "spearman"))
    k <- base & smp$hormone != "MeJA"; add(an, paste0("Pearson r, MeJA groups excluded: ", sn), NA, sum(k), sum(k), length(unique(grp_all[k])), sum(k) - length(unique(grp_all[k])), obs = rr(k))
    k <- base & !smp$sample_uid %in% c("PPPJ8P_19", "PPPJ8P_52", "PPPJ8P_81")
    add(an, paste0("Pearson r, PPPJ8P_19/_52/_81 excluded: ", sn), NA, sum(k), sum(k), length(unique(grp_all[k])), sum(k) - length(unique(grp_all[k])), obs = rr(k),
        note = paste("excluded here:", paste(intersect(c("PPPJ8P_19", "PPPJ8P_52", "PPPJ8P_81"), smp$sample_uid[base]), collapse = ",")))
    lg <- sapply(unique(grp_all[base]), function(gg) rr(base & grp_all != gg))
    add(an, paste0("Pearson r, leave one group out (minimum): ", sn), NA, obs = min(lg), note = sub("^\\d+ ", "", names(which.min(lg))))
    add(an, paste0("Pearson r, leave one group out (maximum): ", sn), NA, obs = max(lg), note = sub("^\\d+ ", "", names(which.max(lg))))
  }
  cat(mn, ": signature genes in matrix", length(gup), "\n")
}
fwrite(rbindlist(resid_rows), file.path(OUT, "within_group_residuals.tsv"), sep = "\t")

S <- rbindlist(summ, fill = TRUE)
S[, random_sets_reproduced_exactly := null_reproduced]
fwrite(S, file.path(OUT, "focal_24h_summary.tsv"), sep = "\t")
print(S[, .(analysis, metric = substr(metric, 1, 70), n_contrasts, n_libraries, n_biological_samples, n_groups, residual_df,
            value = round(value, 3), null_median = round(null_median, 3), null_q975 = round(null_q975, 3), p_upper = round(p_upper, 4))])
print(ct[, .(label, focal = round(focal_log2FC_shrunk_apeglm, 3), mle = round(focal_log2FC_MLE, 3), se = round(focal_lfcSE, 3),
             shift = round(programme_shift, 3), padj10 = signif(programme_shift_padj_BH_across_10, 2))])
print(loo)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
