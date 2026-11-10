#!/usr/bin/env Rscript
# Descriptive QC figures + normalized matrices for the primary v5.1 counts (no hypothesis tests).
#  - normalized matrices: DESeq2 median-ratio normalized counts, log2(CPM+1) (visualization only), VST
#  - PCA (VST, blind, top 2000 variable genes): by batch, by time (ordinal ramp), small multiples by hormone
#  - sample-sample Pearson correlation heatmap (VST)
#  - per-sample QC bars: depth, dedup retention, assignment, strand, organelles
# Palette: validated reference instance (dataviz skill references/palette.md).
# Usage: Rscript scripts/qc_plots.R config/config.yaml
suppressPackageStartupMessages({
  library(yaml); library(data.table); library(DESeq2); library(ggplot2); library(ggrepel)
  library(ComplexHeatmap); library(patchwork); library(edgeR)
})
args <- commandArgs(trailingOnly = TRUE)
cfg <- yaml::read_yaml(ifelse(length(args), args[1], "config/config.yaml")); setwd(cfg$project_root)
odir <- "results/pn40024_v5.1"; fig <- file.path(odir, "figures"); dir.create(fig, recursive = TRUE, showWarnings = FALSE)
cdir <- file.path(odir, "counts"); dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e0"; SURF <- "#fcfcfb"
CAT <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
SEQ <- c("6" = "#9ec5f4", "12" = "#5598e7", "24" = "#256abf", "168" = "#0d366b")
theme_set(theme_minimal(base_size = 10) + theme(text = element_text(colour = INK), axis.text = element_text(colour = INK2),
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
  plot.background = element_rect(fill = SURF, colour = NA), legend.position = "bottom"))
save3 <- function(p, name, w, h) {
  for (ext in c("pdf", "svg")) ggsave(file.path(fig, paste0(name, ".", ext)), p, width = w, height = h, bg = SURF)
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = SURF)
}

cm <- fread("results/pn_v51/counts/pn_v51.primary.counts.tsv")
mat <- as.matrix(cm[, grep("^PPPJ8P_", names(cm)), with = FALSE]); rownames(mat) <- cm$gene_id
s <- fread("metadata/samples.tsv"); setkey(s, sample_uid); s <- s[colnames(mat)]
s[, time_lab := ifelse(is_control, "control", ifelse(time_hours == "168", "7 d", paste0(time_hours, " h")))]
s[, short := sub("PPPJ8P_", "", sample_uid)]
nuc <- !grepl("^MT:", rownames(mat))
dds <- DESeqDataSetFromMatrix(mat[nuc, ], data.frame(row.names = s$sample_uid, batch = factor(s$batch)), ~ 1)
dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized = TRUE)
fwrite(data.table(gene_id = rownames(norm), round(norm, 3)), file.path(cdir, "pn_v51.normalized_counts_deseq2.tsv"), sep = "\t")
lcpm <- cpm(DGEList(mat[nuc, ]), log = TRUE, prior.count = 1)
fwrite(data.table(gene_id = rownames(lcpm), round(lcpm, 4)), file.path(cdir, "pn_v51.log2cpm_visualization_only.tsv"), sep = "\t")
fwrite(data.table(sample_uid = s$sample_uid, size_factor_all_samples = sizeFactors(dds)), file.path(cdir, "pn_v51.size_factors_all_samples.tsv"), sep = "\t")
keep <- rowSums(mat[nuc, ] >= 10) >= 3
vsd <- vst(dds[keep, ], blind = TRUE)
fwrite(data.table(gene_id = rownames(vsd), round(assay(vsd), 4)), file.path(cdir, "pn_v51.vst_blind_visualization_only.tsv"), sep = "\t")

# ---- PCA ----
v <- assay(vsd); top <- head(order(apply(v, 1, var), decreasing = TRUE), 2000)
pc <- prcomp(t(v[top, ])); ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
d <- cbind(s, PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC3 = pc$x[, 3])
fwrite(d[, .(sample_uid, vendor_sample_name, batch, hormone, time_lab, PC1, PC2, PC3)], file.path(fig, "pca_coordinates.tsv"), sep = "\t")
xl <- sprintf("PC1 (%.1f%%)", ve[1]); yl <- sprintf("PC2 (%.1f%%)", ve[2])
pb <- ggplot(d, aes(PC1, PC2, colour = factor(batch))) + geom_point(size = 2.2, stroke = 0) +
  geom_point(data = d[is_control == TRUE], shape = 21, size = 3.6, colour = INK, stroke = 0.5) +
  geom_text_repel(data = d[is_control == TRUE], aes(label = paste0(vendor_sample_name, " (", short, ")")), colour = INK2, size = 2.6, seed = 1) +
  scale_colour_manual(values = CAT[1:2], name = "Batch (= sequencing lane)") + labs(x = xl, y = yl, title = "PCA by batch",
  subtitle = "VST, top 2,000 variable genes; ringed points = untreated controls")
tl <- c(SEQ, control = INK)
pt <- ggplot(d, aes(PC1, PC2, colour = ifelse(is_control, "control", time_hours))) + geom_point(size = 2.2, stroke = 0) +
  scale_colour_manual(values = tl, breaks = c("6", "12", "24", "168", "control"), labels = c("6 h", "12 h", "24 h", "7 d", "control"), name = "Exposure time") +
  labs(x = xl, y = yl, title = "PCA by exposure time", subtitle = "Time is confounded with batch for 7 d (see audit)")
save3(pb + pt, "pca_batch_time", 11, 5)
hs <- c("ABA", "ACC", "GA3", "MeJA", "NAA")
ph <- lapply(seq_along(hs), function(i) {
  ggplot(d, aes(PC1, PC2)) + geom_point(colour = GRID, size = 1.6) +
    geom_point(data = d[is_control == TRUE], colour = INK, size = 1.8, shape = 17) +
    geom_point(data = d[hormone == hs[i]], colour = CAT[1], size = 2) +
    geom_text_repel(data = d[hormone == hs[i]], aes(label = short), size = 2, colour = INK2, max.overlaps = 30, seed = 1) +
    labs(title = hs[i], x = xl, y = yl) + theme(plot.title = element_text(size = 10))
})
save3(wrap_plots(ph, ncol = 3) + plot_annotation(title = "PCA by hormone (blue = hormone; triangles = controls; labels = OSID)"),
      "pca_by_hormone_small_multiples", 12, 8)

# ---- correlation heatmap ----
r <- cor(v)
lab <- paste0(s$vendor_sample_name, " (", s$short, ")")
ha <- HeatmapAnnotation(batch = factor(s$batch), time = s$time_lab, hormone = s$hormone,
  col = list(batch = c("1" = CAT[1], "2" = CAT[2]), time = c("6 h" = SEQ[[1]], "12 h" = SEQ[[2]], "24 h" = SEQ[[3]], "7 d" = SEQ[[4]], control = INK),
             hormone = c(ABA = CAT[1], ACC = CAT[2], GA3 = CAT[3], MeJA = CAT[4], NAA = CAT[5], none = INK)))
col_fun <- circlize::colorRamp2(c(min(r), (min(r) + 1) / 2, 1), c("#f0efec", "#5598e7", "#0d366b"))
for (ext in c("pdf", "png", "svg")) {
  f <- file.path(fig, paste0("sample_correlation_heatmap.", ext))
  if (ext == "png") png(f, 14, 13, units = "in", res = 300) else if (ext == "pdf") pdf(f, 14, 13) else svglite::svglite(f, 14, 13)
  draw(Heatmap(r, name = "Pearson r\n(VST)", col = col_fun, top_annotation = ha, row_labels = lab, column_labels = lab,
               row_names_gp = grid::gpar(fontsize = 5), column_names_gp = grid::gpar(fontsize = 5),
               clustering_method_rows = "average", clustering_method_columns = "average"))
  dev.off()
}

# ---- QC bars ----
q <- fread("results/pn_v51/qc/pn_v51.per_sample_qc.tsv"); q <- merge(q, s[, .(sample_uid, vendor_sample_name, batch, short)], by = "sample_uid")
q[, o := as.integer(short)]; setorder(q, o)
q[, sample := factor(paste0(short, " ", vendor_sample_name), levels = paste0(short, " ", vendor_sample_name))]
long <- melt(q[, .(sample, batch, `Raw reads (M)` = raw_reads / 1e6, `Unique molecules assigned (M)` = fc_sense_assigned / 1e6,
                   `Dedup retained (fraction)` = dedup_retained_frac, `Uniquely mapped (fraction)` = star_unique_frac,
                   `Sense strand (fraction)` = strand_sense_frac, `Organellar (fraction)` = frac_MT + frac_chloroplast,
                   `Reads past annotated 3' end (fraction)` = frac_downstream_of_annotated_end)],
             id.vars = c("sample", "batch"))
pq <- ggplot(long, aes(sample, value, fill = factor(batch))) + geom_col(width = 0.75) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") + scale_fill_manual(values = CAT[1:2], name = "Batch") +
  labs(x = NULL, y = NULL, title = "Per-sample QC (PN40024 T2T v5.1)") +
  theme(axis.text.x = element_text(angle = 90, size = 4.5, hjust = 1, vjust = 0.5))
save3(pq, "per_sample_qc", 16, 14)
cat("qc plots done\n")
