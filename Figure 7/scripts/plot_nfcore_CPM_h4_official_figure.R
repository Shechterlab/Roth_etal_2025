#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd_args[grep("^--file=", cmd_args)][1]
script_path <- if (!is.na(file_arg)) sub("^--file=", "", file_arg) else NA_character_
root_dir <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE) else getwd()

out_dir <- file.path(root_dir, "results", "nfcore_CPM_h4_official_figure")
data_dir <- file.path(out_dir, "data")
fig_dir <- file.path(out_dir, "figures")
stats_dir <- file.path(out_dir, "stats")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
font_cache_dir <- file.path(out_dir, "fontcache")
dir.create(font_cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache_dir)

font_family <- "Arial"
paper_gray <- "#949494"
gsk_green <- "#02783C"
lly_purple <- "#7570b3"
red <- "#d7191c"
hist1_start <- 25723396 / 1e6
hist1_end <- 27893891 / 1e6
pc <- 1e-3

save_plot <- function(plot, base, width, height) {
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 450)
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(paste0(base, ".svg"), width = width, height = height)
  } else {
    grDevices::svg(paste0(base, ".svg"), width = width, height = height, family = font_family)
  }
  print(plot)
  grDevices::dev.off()
}

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(color = "black", family = font_family),
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.32),
      axis.ticks = element_line(linewidth = 0.28),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = base_size + 1.1),
      plot.subtitle = element_text(size = base_size),
      plot.tag = element_text(face = "bold", size = base_size + 4)
    )
}

nice_cap <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(1)
  raw <- max(x)
  if (raw <= 5) return(ceiling(raw))
  if (raw <= 20) return(ceiling(raw / 5) * 5)
  if (raw <= 100) return(ceiling(raw / 10) * 10)
  ceiling(raw / 100) * 100
}

wilcox_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || length(unique(round(x, 12))) < 2) return(NA_real_)
  wilcox.test(x, mu = 0, exact = FALSE)$p.value
}

all_fc <- fread(file.path(data_dir, "nfcore_CPM_H4_PolIIS5ph_drug_24h_48h_expressed_histone_log2FC.csv"))
all_fc[, time_label := factor(gsub("h", " h", time), levels = c("24 h", "48 h"))]
all_fc[, drug := factor(drug, levels = c("GSK591", "LLY283"))]
all_fc[, antibody := factor(antibody, levels = c("H4", "PolIIS5ph"))]

summary_all <- all_fc[is.finite(log2FC_drug_vs_DMSO), .(
  n = .N,
  median_log2FC = median(log2FC_drug_vs_DMSO),
  q1 = as.numeric(quantile(log2FC_drug_vs_DMSO, 0.25)),
  q3 = as.numeric(quantile(log2FC_drug_vs_DMSO, 0.75)),
  wilcoxon_p = wilcox_p(log2FC_drug_vs_DMSO)
), by = .(antibody, drug, time, time_label)]
summary_all[, wilcoxon_fdr := p.adjust(wilcoxon_p, method = "BH"), by = antibody]
fwrite(all_fc, file.path(stats_dir, "nfcore_CPM_H4_PolIIS5ph_drug_24h_48h_expressed_histone_log2FC.csv"))
fwrite(summary_all, file.path(stats_dir, "nfcore_CPM_H4_PolIIS5ph_drug_24h_48h_summary.csv"))

make_gene_panel <- function() {
  genes <- fread(file.path(root_dir, "results", "h4_no_occupancy_change", "data", "hist1_chr6_gene_models.csv"))
  genes[, `:=`(start_mb = start / 1e6, end_mb = end / 1e6, is_label = gene %in% c("H4C1", "H4C2", "H1-3", "H3C7"))]
  plot_genes <- genes[start_mb < 26.285 & end_mb > 26.0]
  label_genes <- plot_genes[is_label == TRUE]
  label_genes[, `:=`(label_x = (start_mb + end_mb) / 2, label_y = -0.34)]
  label_genes[gene == "H4C1", `:=`(label_x = 26.018, label_y = -0.36)]
  label_genes[gene == "H4C2", `:=`(label_x = 26.040, label_y = -0.62)]
  label_genes[gene == "H1-3", `:=`(label_x = 26.224, label_y = -0.36)]
  label_genes[gene == "H3C7", `:=`(label_x = 26.255, label_y = -0.62)]

  ggplot(plot_genes, aes(xmin = start_mb, xmax = end_mb, ymin = -0.16, ymax = 0.16, fill = histone_class)) +
    geom_rect(color = NA, alpha = 0.96) +
    geom_rect(data = plot_genes[is_label == TRUE], aes(xmin = start_mb, xmax = end_mb, ymin = -0.22, ymax = 0.22), inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.22) +
    geom_segment(data = label_genes, aes(x = (start_mb + end_mb) / 2, xend = label_x, y = -0.19, yend = label_y + 0.05), inherit.aes = FALSE, linewidth = 0.17, color = "grey35") +
    geom_text(data = label_genes, aes(x = label_x, y = label_y, label = gene), inherit.aes = FALSE, size = 1.85, family = font_family) +
    scale_fill_manual(values = c(H1 = "#4d4d4d", H2A = "#5aae61", H2B = "#4393c3", H3 = "#fdae61", H4 = red, other = "#bdbdbd"), breaks = c("H1", "H2A", "H2B", "H3", "H4"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(26, 26.285)) +
    scale_y_continuous(limits = c(-0.76, 0.24)) +
    theme_void(base_size = 7, base_family = font_family) +
    theme(legend.position = "bottom", legend.text = element_text(size = 5.4), legend.key.size = unit(0.12, "in"))
}

make_h4_figure <- function(drug, file_base, supplement = FALSE) {
  drug_name <- drug
  drug_col <- if (drug == "GSK591") gsk_green else lly_purple
  fc <- copy(all_fc[antibody == "H4" & as.character(drug) == drug_name])
  set.seed(if (drug == "GSK591") 1801 else 1802)
  fc[, x := as.numeric(time_label) + runif(.N, -0.09, 0.09), by = time_label]
  summary <- copy(summary_all[antibody == "H4" & as.character(drug) == drug_name])
  summary[, x := as.numeric(time_label)]
  summary[, stat_label := sprintf("median %.2f", median_log2FC)]

  # Keep legacy stat filenames for the main GSK591 figure.
  if (drug == "GSK591") {
    old_fc <- fc[, .(
      symbol, region_id, chr, start, end, strand, time,
      dmso_mean_signal,
      gsk591_mean_signal = drug_mean_signal,
      log2FC_GSK591_vs_DMSO = log2FC_drug_vs_DMSO
    )]
    old_summary <- summary[, .(
      time, time_label, n,
      median_log2FC,
      ci_low = q1,
      ci_high = q3,
      wilcoxon_p,
      x,
      stat_label,
      wilcoxon_fdr
    )]
    fwrite(old_fc, file.path(stats_dir, "nfcore_CPM_H4_GSK591_24h_48h_expressed_histone_log2FC.csv"))
    fwrite(old_summary, file.path(stats_dir, "nfcore_CPM_H4_GSK591_24h_48h_summary.csv"))
  }

  p_fc <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.36, color = "grey35") +
    geom_point(data = fc[is.finite(log2FC_drug_vs_DMSO)], aes(x = x, y = log2FC_drug_vs_DMSO),
               shape = 21, size = 1.58, stroke = 0.28, color = paper_gray, fill = paper_gray, alpha = 0.45) +
    geom_errorbar(data = summary, aes(x = x, ymin = q1, ymax = q3), width = 0.13, linewidth = 0.42, color = "black") +
    geom_point(data = summary, aes(x = x, y = median_log2FC), size = 2.15, color = "black") +
    geom_text(data = summary, aes(x = x, y = 4.55, label = stat_label), size = 1.85, family = font_family) +
    annotate("point", x = 0.48, y = -4.45, shape = 21, size = 1.5, stroke = 0.25, color = paper_gray, fill = paper_gray, alpha = 0.45) +
    annotate("text", x = 0.54, y = -4.45, label = "expressed histone genes", hjust = 0, size = 1.65, family = font_family) +
    scale_x_continuous(breaks = c(1, 2), labels = c("24 h", "48 h"), limits = c(0.4, 2.6)) +
    scale_y_continuous(breaks = c(-5, 0, 5)) +
    coord_cartesian(ylim = c(-5.5, 5.5)) +
    labs(
      tag = "A",
      title = "H4 occupancy is maintained",
      subtitle = paste("Expressed histone genes,", paste0(drug, "/DMSO")),
      x = NULL,
      y = "H4 log2FC"
    ) +
    theme_pub(7)

  tracks_wide <- fread(file.path(data_dir, paste0("hist1_chr6p_nfcore_CPM_tracks_", drug, ".csv")))
  tracks_zoom <- fread(file.path(data_dir, paste0("hist1_chr6_zoom_nfcore_CPM_tracks_", drug, ".csv")))
  track_levels <- c("H3 DMSO 24h", "H4 DMSO 24h", paste("H4", drug, "24h"))
  tracks_wide[, `:=`(mb = mid / 1e6, track = factor(track, levels = track_levels))]
  tracks_zoom[, `:=`(mb = mid / 1e6, track = factor(track, levels = track_levels))]
  wide_cap <- nice_cap(tracks_wide[track != "H3 DMSO 24h"]$signal)
  zoom_cap <- nice_cap(tracks_zoom[track != "H3 DMSO 24h"]$signal)
  tracks_wide[, signal_clip := pmin(signal, wide_cap)]
  tracks_zoom[, signal_clip := pmin(signal, zoom_cap)]
  track_cols <- c("H3 DMSO 24h" = "#bdbdbd", "H4 DMSO 24h" = "#111111")
  track_cols[paste("H4", drug, "24h")] <- drug_col
  track_labels <- c("H3 DMSO 24h" = "H3 DMSO", "H4 DMSO 24h" = "H4 DMSO")
  track_labels[paste("H4", drug, "24h")] <- paste("H4", drug)

  p_chr6 <- ggplot() +
    annotate("segment", x = 0, xend = 58, y = 0, yend = 0, linewidth = 1.7, color = "#6b6b6b", lineend = "round") +
    annotate("segment", x = 62.5, xend = 170.8, y = 0, yend = 0, linewidth = 1.7, color = "#6b6b6b", lineend = "round") +
    annotate("rect", xmin = hist1_start, xmax = hist1_end, ymin = -0.14, ymax = 0.14, fill = red, color = NA) +
    annotate("text", x = mean(c(hist1_start, hist1_end)), y = 0.34, label = "HIST1", size = 2.0, fontface = "bold", color = red, family = font_family) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0, 40, 80, 120, 160), limits = c(0, 170.8)) +
    scale_y_continuous(limits = c(-0.44, 0.48)) +
    labs(x = NULL, y = "chr6") +
    theme_classic(base_size = 6, base_family = font_family) +
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_text(angle = 90, size = 6.5))

  rect <- data.table(track = factor(track_levels, levels = track_levels), xmin = hist1_start, xmax = hist1_end, ymin = -Inf, ymax = Inf)
  p_wide <- ggplot(tracks_wide, aes(mb, signal_clip, fill = track, color = track)) +
    geom_rect(data = rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = red, alpha = 0.035, color = NA) +
    geom_area(linewidth = 0.12, alpha = 0.78) +
    facet_grid(track ~ ., switch = "y", labeller = as_labeller(track_labels)) +
    scale_fill_manual(values = track_cols) +
    scale_color_manual(values = track_cols) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0, 20, 40, 58), limits = c(0, 58)) +
    scale_y_continuous(breaks = NULL, limits = c(0, wide_cap), expand = expansion(mult = c(0, 0.03))) +
    annotate("text", x = 0.2, y = wide_cap, label = wide_cap, hjust = 0, vjust = 1.1, size = 1.55, family = font_family) +
    labs(tag = "B", title = paste("H4 CUT&Tag at the chr6p HIST1 locus, 24 h", drug), x = "hg38 chr6 position (Mb)", y = "CUT&Tag signal") +
    theme_pub(7) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, size = 6.2, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  p_zoom <- ggplot(tracks_zoom, aes(mb, signal_clip, fill = track, color = track)) +
    geom_area(linewidth = 0.12, alpha = 0.78) +
    facet_grid(track ~ ., switch = "y", labeller = as_labeller(track_labels)) +
    scale_fill_manual(values = track_cols) +
    scale_color_manual(values = track_cols) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = seq(26.05, 26.25, 0.05), limits = c(26, 26.285)) +
    scale_y_continuous(breaks = NULL, limits = c(0, zoom_cap), expand = expansion(mult = c(0, 0.03))) +
    annotate("text", x = 26.003, y = zoom_cap, label = zoom_cap, hjust = 0, vjust = 1.1, size = 1.45, family = font_family) +
    labs(title = "Representative histone genes", x = "hg38 chr6 position (Mb)", y = NULL) +
    theme_pub(6.5) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, size = 5.8, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  right <- p_chr6 / p_wide / p_zoom / make_gene_panel() + plot_layout(heights = c(0.55, 2.45, 2.55, 0.72))
  final <- p_fc | right
  final <- final + plot_layout(widths = c(0.82, 1.55))
  save_plot(final, file.path(fig_dir, file_base), 8.2, 4.75)
  final
}

make_polii_test_figure <- function() {
  plot_dt <- copy(all_fc[antibody == "PolIIS5ph"])
  set.seed(1810)
  plot_dt[, x := as.numeric(time_label) + ifelse(drug == "GSK591", -0.12, 0.12) + runif(.N, -0.045, 0.045), by = .(drug, time_label)]
  sum_dt <- copy(summary_all[antibody == "PolIIS5ph"])
  sum_dt[, x := as.numeric(time_label) + ifelse(drug == "GSK591", -0.12, 0.12)]
  sum_dt[, label := sprintf("%.2f", median_log2FC)]
  cols <- c(GSK591 = gsk_green, LLY283 = lly_purple)
  p <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.34, color = "grey35") +
    geom_point(data = plot_dt[is.finite(log2FC_drug_vs_DMSO)], aes(x = x, y = log2FC_drug_vs_DMSO, color = drug, fill = drug),
               shape = 21, size = 1.45, stroke = 0.24, alpha = 0.42) +
    geom_errorbar(data = sum_dt, aes(x = x, ymin = q1, ymax = q3), width = 0.08, linewidth = 0.38, color = "black") +
    geom_point(data = sum_dt, aes(x = x, y = median_log2FC), size = 1.9, color = "black") +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = c(1, 2), labels = c("24 h", "48 h"), limits = c(0.55, 2.45)) +
    scale_y_continuous(breaks = c(-2, 0, 2)) +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    labs(title = "PolII S5ph over expressed histone genes", subtitle = "nf-core CPM gene-body signal, drug/DMSO", x = NULL, y = "PolII S5ph log2FC") +
    theme_pub(7) +
    theme(legend.position = "top", legend.title = element_blank(), legend.key.size = unit(0.13, "in"), plot.title = element_text(face = "bold", size = 8.5))
  save_plot(p, file.path(fig_dir, "nfcore_CPM_PolIIS5ph_expressed_histone_24h_48h_drug_vs_DMSO"), 3.7, 3.0)
}

make_h4_figure("GSK591", "nfcore_CPM_h4_occupancy_maintained_HIST1_hg38")
make_h4_figure("GSK591", "nfcore_CPM_h4_occupancy_maintained_HIST1_hg38_GSK591")
make_h4_figure("LLY283", "nfcore_CPM_h4_occupancy_maintained_HIST1_hg38_LLY283_supplement", supplement = TRUE)
make_polii_test_figure()

cat(
  "# nf-core CPM H4 Official Figure Set\n\n",
  "These figures remake the H4 occupancy/HIST1 analysis using averaged nf-core/cutandrun 3.2.2 CPM-normalized refgenie hg38 bigWigs from `../bw/averaged/*_mean.bigWig`.\n\n",
  "Main figure: `figures/nfcore_CPM_h4_occupancy_maintained_HIST1_hg38.pdf`.\n\n",
  "Supplemental LLY283 figure: `figures/nfcore_CPM_h4_occupancy_maintained_HIST1_hg38_LLY283_supplement.pdf`.\n\n",
  "PolII S5ph test figure: `figures/nfcore_CPM_PolIIS5ph_expressed_histone_24h_48h_drug_vs_DMSO.pdf`.\n\n",
  "Summary medians:\n\n",
  paste(sprintf(
    "- %s %s %s: n=%d, median log2FC %.3f, IQR %.3f to %.3f, Wilcoxon FDR %.3g",
    summary_all$antibody,
    summary_all$drug,
    summary_all$time,
    summary_all$n,
    summary_all$median_log2FC,
    summary_all$q1,
    summary_all$q3,
    summary_all$wilcoxon_fdr
  ), collapse = "\n"),
  "\n\nCoordinate note: this version uses the new nf-core/refgenie hg38 tracks and should supersede figures made from older Ahmad bigWigs.\n",
  sep = "",
  file = file.path(out_dir, "README_nfcore_CPM_h4_official_figure.md")
)

message("Done: ", out_dir)
