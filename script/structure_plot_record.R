# The per-lineage structure-plot record. NOT a figure: this script publishes
# nothing.
#
# It builds the seven per-lineage structure plots -- one per lineage, CD8, CD4,
# Treg, gdT, CD8aa, Tz, DN, each over every GP that reaches AUC > 0.9 for at
# least one cluster of that lineage -- and writes down which GPs, colours and
# clusters each row uses, into output/structure_plot_record/.
#
# Two things consume that record:
#   script/Figure4.R                     panels 4b, 4c and 4d are these GPs and
#                                        these cells in another form, and the
#                                        script stops if they disagree with what
#                                        is recorded here.
#   script/verify_structure_plot_gps.R   re-derives every row's GP set from the
#                                        AUCs published as Extended Data Table 6
#                                        and diffs it against the record.
#
# The row map, the AUC rule, the display filters and the palette all come from
# code/R/structure_plot_panels.R, which Figure4.R uses too.
#
# --- internal ---
# Was script/FigureS4.R, which drew the seven rows stacked into one giant
# Extended Data figure. Ziang removed that figure on 2026-09-10; the script
# survives because Figure 4 checks itself against what it records. The rows are
# still built rather than short-circuited, so `n_cells_drawn` stays "what a row
# draws" rather than "what the rule says" -- the assembly and the ggsave are all
# that were taken out, and restoring the figure means putting them back.
# --- end internal ---
#
# --- internal ---
# Ported from experiments/giant_structure_plot_by_lineage/save_separate_pdfs.R
# and the "Giant structure plot" section of experiments/assess_structure_plot.R,
# whose stacked layout (16 in wide, 3 in per row, plot_grid(align = "v")) this
# figure follows. Two changes from the exploratory rows: DP is dropped (seven
# rows, not eight), and each row is coloured on its own rather than from one
# GP -> colour map shared by all rows, which is what those PDFs and this figure
# did until 2026-09-02 -- see the palette comment in
# code/R/structure_plot_panels.R and the trials in
# experiments/structure_plot_recolor/, whose per_lineage_glasbey variant this
# figure now reproduces exactly (RMSE 0; script/README.md records the check).
# The cells, clusters, GP sets and geometry have never changed.
#
# Took the Extended Data Figure 5 slot on 2026-08-27, when it moved there from
# Figure S8: the protein-program heatmap that had been Figure S5 became Figure
# S6, and the CD69/gating figure became Figure S7. Renumbered to Extended Data
# Figure 4 on 2026-09-10, when Extended Data Figure 3 was folded into Extended
# Data Figure 2 as its panel d and everything after it moved up one.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md's "Data provenance" table for
# the full picture:
#   igt1_96_..._ADTonly.Rds                           [primary input Seurat object]
#   L_pm_filtered.rds                                 [code/pipeline/01b_filter_cells.R]
#   level_2_AUC_list_figure_no_thymocytes_healthy.rds [code/pipeline/02_compute_auc.R]

# --- doc:setup ---
library(ggplot2)
library(ggrastr)
library(cowplot)
library(fastTopics) # structure_plot()

if (!file.exists("code/R/structure_plot_panels.R")) {
  stop("Run this script from the immgenT-GP-analysis repository root.")
}
source("code/R/structure_plot_panels.R")

data_path <- "data/"
record_path <- "output/structure_plot_record/"
dir.create(record_path, recursive = TRUE, showWarnings = FALSE)

# Record when this run started, to assert at the end that the figure is newer.
# --- internal ---
# A figure script here was once seen to exit 0 with a complete log and write
# nothing at all -- see script/README.md, "A re-run can silently not write".
# --- end internal ---
run_started_at <- Sys.time()

# ============================================================
# Load data
# ============================================================
seurat_meta <- readRDS(paste0(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
colnames(L_pm_filtered) <- gsub("^K", "GP", colnames(L_pm_filtered))

level_2_AUC_list <- readRDS(paste0(
  data_path,
  "level_2_AUC_list_figure_no_thymocytes_healthy.rds"
))
auc_level2 <- level_2_AUC_list$auc
colnames(auc_level2) <- gsub("^K", "GP", colnames(auc_level2))

# The AUC matrix and the loading matrix have to be talking about the same GPs
# in the same order: the panels index one by names taken from the other.
if (!identical(colnames(auc_level2), colnames(L_pm_filtered))) {
  stop("The AUC matrix and L_pm_filtered do not carry the same GP columns in the same order.")
}

level1_all <- seurat_meta_filtered$annotation_level1
level2_all <- seurat_meta_filtered$annotation_level2
healthy_non_thymocyte <- which(
  seurat_meta_filtered$condition_broad == "healthy" &
    level1_all != "thymocyte"
)

# ============================================================
# GP selection
# ============================================================
# Which lineage each cluster belongs to, from the metadata.
cluster_lineage <- cluster_lineage_map(level2_all, level1_all)

# --- internal ---
# verify_structure_plot_gps.R has only the published table to work from, so it
# maps clusters to lineages by their name prefix (CD8.A -> CD8) instead. Check
# that shortcut here, where the metadata-derived map is available, so the check
# cannot be re-deriving a different grouping than the figure drew.
# --- end internal ---
auc_clusters <- rownames(auc_level2)
prefix_lineage <- sub("[.].*$", "", auc_clusters)
if (!identical(unname(cluster_lineage[auc_clusters]), prefix_lineage)) {
  disagree <- auc_clusters[unname(cluster_lineage[auc_clusters]) != prefix_lineage]
  stop(sprintf(
    "cluster name prefixes disagree with annotation_level1 for: %s",
    paste(disagree, collapse = ", ")
  ))
}

# One GP set per row. Each row's palette is then assigned inside the loop
# below, independently of the other rows -- see structure_plot_row_colors().
panel_gps <- gps_above_auc_by_lineage(auc_level2, cluster_lineage)

message(sprintf(
  "%d GPs over %d rows (AUC > %.1f): %s",
  length(unique(unlist(panel_gps, use.names = FALSE))), length(panel_gps),
  structure_plot_auc_threshold,
  paste(sprintf("%s %d", names(panel_gps), lengths(panel_gps)), collapse = ", ")
))

# ============================================================
# rows: one structure plot per lineage, built to record what it uses
# ============================================================
# The loop iterates over the names of the row map, so a lineage cannot be drawn
# under another lineage's letter. Each row is kept as a ggplot and the rows are
# assembled below, rather than saved one file per lineage.
# --- internal ---
# The figure is the stack, and assembling it here means no hand layout step can
# fall behind a re-run.
# --- end internal ---
lineage_plots <- list()
cluster_records <- list()
gp_records <- list()

for (lineage in names(structure_plot_panels)) {
  panel <- structure_plot_panels[[lineage]]
  gps_lineage <- panel_gps[[lineage]]

  lineage_cells <- healthy_non_thymocyte[level1_all[healthy_non_thymocyte] == lineage]
  cluster_size <- table(droplevels(factor(level2_all[lineage_cells])))
  small_clusters <- names(cluster_size)[cluster_size < structure_plot_min_cluster_cells]
  lineage_cells <- lineage_cells[!level2_all[lineage_cells] %in% small_clusters]

  # Cap each cluster's width so one large cluster cannot crowd out the rest.
  set.seed(1234)
  keep <- unlist(lapply(
    split(seq_along(lineage_cells), level2_all[lineage_cells]),
    function(idx) {
      if (length(idx) > structure_plot_max_cells_per_cluster) {
        sample(idx, structure_plot_max_cells_per_cluster)
      } else {
        idx
      }
    }
  ))
  lineage_cells <- lineage_cells[keep]

  fit_lineage <- L_pm_filtered[lineage_cells, gps_lineage, drop = FALSE]
  grouping_lineage <- factor(level2_all[lineage_cells])

  # This row's colors, assigned from the top of the palette without reference to
  # any other row. structure_plot_row_colors() returns them in ascending GP
  # order, which is the column order here.
  # --- internal ---
  # structure_plot() renames the colors it is given positionally, by the columns
  # of the matrix, so a palette in any other order would mislabel every bar.
  # --- end internal ---
  colors_lineage <- structure_plot_row_colors(gps_lineage)
  if (!identical(names(colors_lineage), colnames(fit_lineage))) {
    stop(sprintf("panel %s: palette order does not match its GP columns.", panel))
  }

  set.seed(1234)
  p <- structure_plot(
    fit_lineage,
    topics = gps_lineage,
    gap = 40,
    n = 10000,
    colors = colors_lineage,
    grouping = grouping_lineage,
    ggplot_call = rasterized_structure_plot_call
  ) +
    labs(
      y = "membership",
      color = "",
      fill = "",
      title = sprintf(
        "%s (%d GPs, AUC > %.1f)",
        lineage,
        length(gps_lineage),
        structure_plot_auc_threshold
      )
    ) +
    guides(
      fill = guide_legend(ncol = 2),
      color = guide_legend(ncol = 2)
    ) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 10, face = "bold"),
      legend.position = "right",
      legend.key.size = unit(0.25, "cm"),
      legend.text = element_text(size = 5),
      legend.spacing.y = unit(0.02, "cm")
    )

  lineage_plots[[lineage]] <- p

  # What this row used, for the alignment check. Every cluster of the lineage
  # is listed, drawn or not: the GP selection above uses all of them.
  lineage_clusters <- sort(auc_clusters[prefix_lineage == lineage])
  drawn_size <- table(droplevels(factor(level2_all[lineage_cells])))
  cluster_records[[lineage]] <- data.frame(
    panel = panel,
    lineage = lineage,
    cluster = lineage_clusters,
    n_cells_healthy = as.integer(cluster_size[lineage_clusters]),
    n_cells_drawn = as.integer(ifelse(
      lineage_clusters %in% names(drawn_size),
      drawn_size[lineage_clusters],
      0L
    )),
    stringsAsFactors = FALSE
  )
  gp_records[[lineage]] <- data.frame(
    panel = panel,
    lineage = lineage,
    gp = gps_lineage,
    color = unname(colors_lineage),
    max_auc_in_lineage = unname(apply(
      auc_level2[lineage_clusters, gps_lineage, drop = FALSE], 2,
      max, na.rm = TRUE
    )),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# What each panel drew, for the alignment check
# ============================================================
# --- internal ---
# script/verify_structure_plot_gps.R reads these two files, so that it can
# re-derive the panels' GP sets from the published Extended Data Table 6 without
# reloading the 1 GB loading matrix, and check the clusters drawn and omitted
# against the display filters above.
# --- end internal ---
cluster_record <- do.call(rbind, cluster_records)
cluster_record$n_cells_healthy[is.na(cluster_record$n_cells_healthy)] <- 0L
write.csv(
  cluster_record,
  paste0(record_path, "panel_clusters.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, gp_records),
  paste0(record_path, "panel_gps.csv"),
  row.names = FALSE
)

# ============================================================
# Did this run actually write the record?
# ============================================================
# A write that fails leaves the previous CSV in place and the script still exits
# 0, so check the files rather than the exit code -- Figure 4 would then be
# checking itself against a stale record.
expected_records <- paste0(record_path, c("panel_clusters.csv", "panel_gps.csv"))
stale <- expected_records[!file.exists(expected_records) |
                            file.mtime(expected_records) < run_started_at |
                            file.size(expected_records) == 0]
if (length(stale) > 0L) {
  stop("not written by this run: ", paste(stale, collapse = ", "))
}
message(sprintf(
  "wrote the %d-row structure-plot record to %s",
  length(structure_plot_panels), record_path
))
