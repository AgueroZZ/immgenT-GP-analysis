# Figure S6. Cluster-level GP membership within each T cell lineage.
#
# One stacked figure (s6.pdf): seven rows, one per lineage, labelled
# a = CD8, b = CD4, c = Treg, d = gdT, e = CD8aa, f = Tz, g = DN. Each row is a
# structure plot of the healthy non-thymocyte cells of that lineage, grouped by
# their annotation_level2 cluster, over every GP that reaches AUC > 0.9 for at
# least one cluster of that lineage. It is the per-cluster, all-GP counterpart
# of Figure 3B, which shows six hand-picked lineage-defining GPs grouped by
# lineage.
#
# The row map, the AUC rule, the display filters, the palette and the assembled
# geometry all come from code/R/structure_plot_panels.R. This script records what
# it drew into output/FigureS6/, and script/verify_structure_plot_gps.R checks
# that record against the AUCs published as Extended Data Table 6 and against the
# caption on analysis/FigureS6.Rmd.
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
# S6, and the CD69/gating figure became Figure S7. Moved on again to Extended
# Data Figure 6 on 2026-09-09, when a new tissue figure was inserted at 5; only
# numbers and paths changed, nothing was re-rendered.
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
figure_path <- "figures/final-selected/Figure S6/"
record_path <- "output/FigureS6/"
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)
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
# s6: per-lineage rows, stacked into one figure
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

# Rows are stacked in the map's order and labelled a-g. align = "v" equalises
# everything outside the plotting panel -- y-axis labels and the per-row legends,
# which differ in width because the rows show 9 to 44 GPs -- so the cluster
# blocks line up down the figure instead of each row starting at its own x.
p_s6 <- cowplot::plot_grid(
  plotlist = lineage_plots[names(structure_plot_panels)],
  nrow = length(structure_plot_panels),
  align = "v",
  labels = "auto",
  label_size = 14
)
ggsave(
  filename = paste0(figure_path, "s6.pdf"),
  plot = p_s6,
  width = structure_plot_width,
  height = structure_plot_row_height * length(structure_plot_panels),
  dpi = 300,
  limitsize = FALSE
)

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
  paste0(record_path, "s6_panel_clusters.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, gp_records),
  paste0(record_path, "s6_panel_gps.csv"),
  row.names = FALSE
)

# ============================================================
# Did this run actually write the figure?
# ============================================================
expected_panel <- paste0(figure_path, "s6.pdf")
if (!file.exists(expected_panel) ||
      file.mtime(expected_panel) < run_started_at ||
      file.size(expected_panel) == 0) {
  stop(sprintf(
    "%s was not written by this run (missing, empty, or older than the run)",
    expected_panel
  ))
}
message(sprintf(
  "wrote %s (%d rows, %.0f x %.0f in)",
  expected_panel, length(structure_plot_panels), structure_plot_width,
  structure_plot_row_height * length(structure_plot_panels)
))
