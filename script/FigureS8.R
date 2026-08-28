# Figure S8. Cluster-level GP membership within each T cell lineage.
#
# Seven panels, one per lineage (s8a = CD8, s8b = CD4, s8c = Treg, s8d = gdT,
# s8e = CD8aa, s8f = Tz, s8g = DN). Each is a structure plot of the healthy
# non-thymocyte cells of that lineage, grouped by their annotation_level2
# cluster, over every GP that reaches AUC > 0.9 for at least one cluster of
# that lineage. It is the per-cluster, all-GP counterpart of Figure 3B, which
# shows six hand-picked lineage-defining GPs grouped by lineage.
#
# The panel map, the AUC rule, the display filters and the palette all come
# from code/R/structure_plot_panels.R. This script records what it drew into
# output/FigureS8/, and script/verify_structure_plot_gps.R checks that record
# against the AUCs published as Extended Data Table 6 and against the caption
# on analysis/FigureS8.Rmd.
#
# --- internal ---
# Ported from experiments/giant_structure_plot_by_lineage/save_separate_pdfs.R,
# with two changes: the DP panel is dropped (seven panels, not eight), and the
# palette is assigned over the 69 GPs these seven panels show rather than the 72
# that pass in any lineage -- GP9, GP14 and GP48 pass only in DP clusters. Both
# changes move every color, so these panels are not pixel-comparable with those
# exploratory PDFs; the cells, clusters, GP sets and geometry are unchanged.
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
figure_path <- "figures/final-selected/Figure S8/"
record_path <- "output/FigureS8/"
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)
dir.create(record_path, recursive = TRUE, showWarnings = FALSE)

# A figure script here was once seen to exit 0 with a complete log and write
# nothing at all (see script/README.md, "A re-run can silently not write"), so
# record when this run started and assert at the end that the panels are newer.
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
# GP selection and the shared palette
# ============================================================
# Which lineage each cluster belongs to, from the metadata.
cluster_lineage <- cluster_lineage_map(level2_all, level1_all)

# verify_structure_plot_gps.R has only the published table to work from, so it
# maps clusters to lineages by their name prefix (CD8.A -> CD8) instead. Check
# that shortcut here, where the metadata-derived map is available, so the check
# cannot be re-deriving a different grouping than the figure drew.
auc_clusters <- rownames(auc_level2)
prefix_lineage <- sub("[.].*$", "", auc_clusters)
if (!identical(unname(cluster_lineage[auc_clusters]), prefix_lineage)) {
  disagree <- auc_clusters[unname(cluster_lineage[auc_clusters]) != prefix_lineage]
  stop(sprintf(
    "cluster name prefixes disagree with annotation_level1 for: %s",
    paste(disagree, collapse = ", ")
  ))
}

# One GP set per panel, and one color per GP over their union -- so a GP that
# marks clusters in two lineages keeps its color in both panels.
panel_gps <- gps_above_auc_by_lineage(auc_level2, cluster_lineage)
gp_colors <- structure_plot_gp_colors(unlist(panel_gps, use.names = FALSE))

message(sprintf(
  "%d GPs over %d panels (AUC > %.1f): %s",
  length(gp_colors), length(panel_gps), structure_plot_auc_threshold,
  paste(sprintf("%s %d", names(panel_gps), lengths(panel_gps)), collapse = ", ")
))

# ============================================================
# s8a-s8g: per-lineage structure plots
# ============================================================
# The loop iterates over the names of the panel map, so a lineage cannot be
# drawn under another lineage's letter.
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

  # structure_plot() renames the colors it is given positionally, by the columns
  # of the matrix, so a palette in any other order would mislabel every bar.
  colors_lineage <- gp_colors[gps_lineage]
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

  ggsave(
    filename = paste0(figure_path, panel, ".pdf"),
    plot = p,
    width = 11,
    height = 4,
    dpi = 300,
    limitsize = FALSE
  )

  # What this panel used, for the alignment check. Every cluster of the lineage
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
    color = unname(gp_colors[gps_lineage]),
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
# script/verify_structure_plot_gps.R reads these two files, so that it can
# re-derive the panels' GP sets from the published Extended Data Table 6 without
# reloading the 1 GB loading matrix, and check the clusters drawn and omitted
# against the display filters above.
cluster_record <- do.call(rbind, cluster_records)
cluster_record$n_cells_healthy[is.na(cluster_record$n_cells_healthy)] <- 0L
write.csv(
  cluster_record,
  paste0(record_path, "s8_panel_clusters.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, gp_records),
  paste0(record_path, "s8_panel_gps.csv"),
  row.names = FALSE
)

# ============================================================
# Did this run actually write the panels?
# ============================================================
expected_panels <- paste0(figure_path, structure_plot_panels, ".pdf")
stale <- expected_panels[!file.exists(expected_panels) |
                           file.mtime(expected_panels) < run_started_at |
                           file.size(expected_panels) == 0]
if (length(stale)) {
  stop(sprintf(
    "these panels were not written by this run (missing, empty, or older than the run): %s",
    paste(basename(stale), collapse = ", ")
  ))
}
message(sprintf("wrote %d panels to %s", length(expected_panels), figure_path))
