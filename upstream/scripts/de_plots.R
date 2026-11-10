#!/usr/bin/env Rscript
# Figures for DE results (run after scripts/de_analysis.R).
#  - DEG summary (up/down per contrast, faceted by batch)
#  - volcano + MA for selected contrasts (largest n_sig per hormone)
#  - focal gene: effect-size heatmap (hormone x dose rows, time columns; diverging blue-gray-red,
#    non-estimable cells marked 'NA'; significance shown ONLY by an explicit symbol, never by color),
#    forest plot with 95% Wald CIs, replicate-level dot/box plot with batch-matched controls
# Palette: validated reference instance (dataviz references/palette.md).
suppressPackageStartupMessages({ library(yaml); library(data.table); library(ggplot2); library(ggrepel); library(patchwork); library(DESeq2) })
args <- commandArgs(trailingOnly = TRUE)
cfg <- yaml::read_yaml(ifelse(length(args), args[1], "config/config.yaml")); setwd(cfg$project_root)
FOCAL <- "Vitvi05_01chr01g05750"; FDR <- cfg$design$fdr
de <- "results/pn40024_v5.1/de"; fig <- file.path(de, "figures"); dir.create(fig, recursive = TRUE, showWarnings = FALSE)
tg <- "results/pn40024_v5.1/target_gene"
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e0"; SURF <- "#fcfcfb"; NEU <- "#f0efec"
CAT <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
BLUE <- "#104281"; RED <- "#b3261e"
theme_set(theme_minimal(base_size = 10) + theme(text = element_text(colour = INK), axis.text = element_text(colour = INK2),
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
  plot.background = element_rect(fill = SURF, colour = NA), legend.position = "bottom"))
save3 <- function(p, dir, name, w, h) {
  for (ext in c("pdf", "svg")) ggsave(file.path(dir, paste0(name, ".", ext)), p, width = w, height = h, bg = SURF)
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = SURF)
}
tl <- function(h) { h <- suppressWarnings(as.numeric(h)); ifelse(h == 168, "7 d", paste0(h, " h")) }

cs <- fread(file.path(de, "contrast_summary.tsv"))
cs[, time_lab := factor(tl(time_hours), levels = c("6 h", "12 h", "24 h", "7 d"))]
cs[, dose_lab := paste(hormone, dose_value, sub("uM", "µM", dose_unit))]
long <- rbind(cs[, .(dose_lab, hormone, time_lab, batch, n = n_up, dir = "up (treated > control)")],
              cs[, .(dose_lab, hormone, time_lab, batch, n = -n_down, dir = "down")])
p <- ggplot(long, aes(time_lab, n, fill = dir)) + geom_col(width = 0.7) + geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
  facet_wrap(~ dose_lab, ncol = 5) + scale_fill_manual(values = c("up (treated > control)" = RED, down = BLUE), name = NULL) +
  labs(x = "Exposure time (control set = same batch)", y = sprintf("DE genes (BH FDR < %.2f; down plotted negative)", FDR),
       title = "Differentially expressed genes per hormone-dose-time group vs untreated callus",
       subtitle = "DESeq2 within-batch Wald tests. 7 d groups except ABA 5 mg/L are in batch 2 (different controls) - do not read as a time course")
save3(p, fig, "deg_summary", 13, 7)

# volcano/MA: the contrast with most DE genes per hormone
sel <- cs[order(-n_sig), .SD[1], by = hormone]
vp <- list(); mp <- list()
for (i in seq_len(nrow(sel))) {
  t <- fread(file.path(de, "tables", paste0(sel$contrast[i], ".all_genes.tsv")))
  t <- t[!is.na(padj_BH)]
  t[, cls := fifelse(padj_BH < FDR & log2FC > 0, "up", fifelse(padj_BH < FDR & log2FC < 0, "down", "ns"))]
  lab <- t[gene_id == FOCAL]
  vp[[i]] <- ggplot(t, aes(log2FC_shrunk_apeglm, -log10(pvalue), colour = cls)) + geom_point(size = 0.6, stroke = 0) +
    scale_colour_manual(values = c(up = RED, down = BLUE, ns = "#c9c8c2"), guide = "none") +
    geom_point(data = lab, colour = INK, shape = 21, size = 2.6, stroke = 0.8) +
    geom_text_repel(data = lab, aes(label = "focal gene"), colour = INK, size = 2.6, seed = 1) +
    labs(title = sel$label[i], x = "shrunken log2FC (treated / control)", y = "-log10 p")
  mp[[i]] <- ggplot(t, aes(log10(baseMean), log2FC, colour = cls)) + geom_point(size = 0.6, stroke = 0) +
    scale_colour_manual(values = c(up = RED, down = BLUE, ns = "#c9c8c2"), guide = "none") +
    geom_point(data = lab, colour = INK, shape = 21, size = 2.6, stroke = 0.8) + geom_hline(yintercept = 0, colour = INK2, linewidth = 0.3) +
    labs(title = sel$label[i], x = "log10 mean normalized count", y = "log2FC (MLE)")
}
save3(wrap_plots(vp, ncol = 3) + plot_annotation(title = "Volcano plots (contrast with most DE genes per hormone); red/blue = BH FDR < 0.05"),
      fig, "volcano_selected", 13, 8)
save3(wrap_plots(mp, ncol = 3) + plot_annotation(title = "MA plots (same contrasts)"), fig, "ma_selected", 13, 8)

# ---------------- focal gene ----------------
ff <- fread(file.path(tg, "focal_gene_de_40_contrasts.tsv"))
s <- fread("metadata/samples.tsv")
lab_map <- unique(s[is_control == FALSE, .(condition, hormone, dose_value, dose_unit, time_hours)])
ff <- merge(ff[, !c("hormone", "dose_value", "dose_unit", "time_hours"), with = FALSE], lab_map, by = "condition", all.x = TRUE)
ff[, time_lab := factor(tl(time_hours), levels = c("6 h", "12 h", "24 h", "7 d"))]
ff[, dose_lab := factor(paste(hormone, dose_value, sub("uM", "µM", dose_unit)))]
ff[, sig_mark := fifelse(is.na(padj_BH_across_40_focal_tests), "NA",
                 fifelse(padj_BH_across_40_focal_tests < 0.001, "***", fifelse(padj_BH_across_40_focal_tests < 0.01, "**",
                 fifelse(padj_BH_across_40_focal_tests < 0.05, "*", ""))))]
grid <- CJ(dose_lab = levels(ff$dose_lab), time_lab = levels(ff$time_lab))
hm <- merge(grid, ff[, .(dose_lab = as.character(dose_lab), time_lab = as.character(time_lab), log2FC, sig_mark, batch)],
            by = c("dose_lab", "time_lab"), all.x = TRUE)
hm[, time_lab := factor(time_lab, levels = c("6 h", "12 h", "24 h", "7 d"))]
lim <- max(abs(hm$log2FC), na.rm = TRUE)
ph <- ggplot(hm, aes(time_lab, dose_lab)) + geom_tile(aes(fill = log2FC), colour = SURF, linewidth = 1) +
  geom_text(aes(label = ifelse(is.na(log2FC), "not\nestimable", sprintf("%.2f%s", log2FC, sig_mark)),
                colour = !is.na(log2FC) & abs(log2FC) > 0.7 * lim), size = 2.8, lineheight = 0.8, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#ffffff", `FALSE` = INK), guide = "none") +
  scale_fill_gradient2(low = BLUE, mid = NEU, high = RED, midpoint = 0, limits = c(-lim, lim), na.value = "#ffffff",
                       name = "log2FC (treated / batch-matched control)") +
  labs(x = "Exposure time", y = NULL, title = sprintf("Focal gene %s: effect sizes (DESeq2 MLE log2FC)", FOCAL),
       subtitle = "Stars = BH-adjusted across the 40 focal tests (* <0.05, ** <0.01, *** <0.001). Colour shows effect size only, not significance.\n7 d groups (except ABA 5 mg/L) use batch-2 controls; columns are not a time course.") +
  theme(panel.grid = element_blank(), legend.position = "right")
save3(ph, tg, "focal_effect_heatmap", 10, 6.5)
fo <- ff[!is.na(log2FC)][order(hormone, dose_value, time_hours)]
fo[, y := factor(label, levels = rev(unique(label)))]
pf <- ggplot(fo, aes(log2FC, y)) + geom_vline(xintercept = 0, colour = INK2, linewidth = 0.3) +
  geom_errorbar(aes(xmin = ci95_low, xmax = ci95_high), width = 0, orientation = "y", colour = INK2, linewidth = 0.5) +
  geom_point(aes(colour = padj_BH_across_40_focal_tests < FDR), size = 2) +
  scale_colour_manual(values = c(`TRUE` = RED, `FALSE` = "#8f8e88"), labels = c(`TRUE` = "BH (40 focal tests) < 0.05", `FALSE` = "not significant"), name = NULL) +
  labs(x = "log2FC treated / control (MLE, 95% Wald CI)", y = NULL, title = sprintf("Focal gene %s across 40 contrasts", FOCAL))
save3(pf, tg, "focal_forest", 8, 10)
# replicate-level dot/box plot with DESeq2-normalized counts (per-batch size factors)
nb <- rbindlist(lapply(c("1", "2"), function(b) {
  d <- readRDS(file.path(de, sprintf("dds_within_batch%s.rds", b)))
  if (!FOCAL %in% rownames(d)) return(NULL)
  data.table(sample_uid = colnames(d), norm = counts(d, normalized = TRUE)[FOCAL, ], batch = b)
}))
nb <- merge(nb, s[, .(sample_uid, condition, hormone, dose_value, dose_unit, time_hours, is_control, vendor_sample_name)], by = "sample_uid")
nb[, grp := ifelse(is_control, paste0("control (batch ", batch, ")"), paste0(hormone, " ", dose_value, " ", sub("uM", "µM", dose_unit), ", ", tl(as.numeric(time_hours))))]
ord <- c("control (batch 1)", unique(nb[is_control == FALSE][order(hormone, dose_value, as.numeric(time_hours))]$grp), "control (batch 2)")
nb[, grp := factor(grp, levels = ord)]
fwrite(nb, file.path(tg, "focal_normalized_expression_per_sample.tsv"), sep = "\t")
fwrite(nb[, .(n = .N, mean_norm = mean(norm), median_norm = median(norm), sd_norm = sd(norm), min = min(norm), max = max(norm)), by = .(grp, batch)],
       file.path(tg, "focal_group_summary.tsv"), sep = "\t")
pd <- ggplot(nb, aes(grp, norm)) + geom_boxplot(outlier.shape = NA, width = 0.55, colour = INK2, linewidth = 0.3, fill = NA) +
  geom_point(aes(colour = factor(batch)), size = 1.8, position = position_jitter(width = 0.12, seed = 1)) +
  scale_colour_manual(values = CAT[1:2], name = "Batch") + scale_y_log10() +
  labs(x = NULL, y = "DESeq2-normalized molecules (per-batch size factors, log scale)",
       title = sprintf("Focal gene %s: individual replicates", FOCAL)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
save3(pd, tg, "focal_replicates_dotplot", 14, 6)
cat("DE plots done\n")
