# Figure S5 (data step): EBMF cluster-mean matrix and column/palette metadata for
# the EBMF vs RQVI gene-program comparison across annotation_level2 clusters.
#
#   * Cells: L_pm_filtered gene-program loadings (cells passing the iterative
#     total-loading filter), restricted to non-thymocytes (annotation_level1 !=
#     "thymocyte", all conditions) and to cells that also carry an RQVI loading
#     ("common cells"). The RQVI cell set is taken from the cell_id index of the
#     RQVI loading table.
#   * EBMF matrix: flashier loadings L_pm_filtered (GP1..GP200), averaged within
#     annotation_level2 on the common cells.
#   * Columns (level2 clusters) are ordered by the Figure-1 level1 lineage order,
#     alphabetically within each lineage.
#
# This script writes the EBMF cluster-mean matrix, a cell->annotation table, and
# the column order + lineage palette. The RQVI cluster means and the one-to-one
# EBMF-RQVI matching are computed by script/FigureS5_rematch.py; row ordering,
# per-program [0,1] scaling, and the heatmaps by script/FigureS5_plot.py.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(ZemmourLib)
})

if (!file.exists("code/R/setup_data.R")) {
  stop("Run this script from the immgenT-GP-analysis repository root.")
}
source("code/R/setup_data.R")

outdir <- "figures/generated/Figure S5"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

PKG <- "data/rqvi_loading/RQVI_EBMF_heatmap_data_v1/data"
rqvi_loading_path <- file.path(PKG, "rqvi_matched_200_cell_loadings.parquet")
level1_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN", "DP")

## ---- EBMF loadings + metadata ----
gp <- load_gp_data()
L  <- gp$L_pm_filtered
colnames(L) <- paste0("GP", seq_len(ncol(L)))
meta <- gp$seurat_meta_filtered
stopifnot(identical(rownames(L), rownames(meta)), ncol(L) == 200L)

nonthy <- meta$annotation_level1 != "thymocyte"
if (anyNA(nonthy)) stop("annotation_level1 has missing values.")

## ---- common cells (non-thymocyte, in L_pm, and with an RQVI loading) ----
rqvi_cell_ids <- as.data.frame(arrow::read_parquet(rqvi_loading_path, col_select = "cell_id"))$cell_id
our_nonthy_ids <- rownames(L)[nonthy]
common_ids <- our_nonthy_ids[our_nonthy_ids %in% rqvi_cell_ids]
grp <- droplevels(factor(as.character(meta$annotation_level2[match(common_ids, rownames(meta))])))
if (anyNA(grp) || any(as.character(grp) == "")) stop("Missing level2 labels on common cells.")

message(sprintf("common non-thymocyte cells: %d (healthy %d / non-healthy %d); level2 clusters: %d",
                length(common_ids),
                sum(meta$condition_broad[match(common_ids, rownames(meta))] == "healthy"),
                sum(meta$condition_broad[match(common_ids, rownames(meta))] != "healthy"),
                nlevels(grp)))

## ---- EBMF cluster means (columns F1..F200) ----
Lc <- L[common_ids, , drop = FALSE]
colnames(Lc) <- paste0("F", seq_len(ncol(Lc)))
esum <- rowsum(Lc, grp)
ecnt <- as.integer(table(grp)[rownames(esum)])
ebmf_means <- sweep(esum, 1L, ecnt, "/")                 # K x 200
stopifnot(nrow(ebmf_means) == 107L, ncol(ebmf_means) == 200L)

## ---- level2 column order (Figure-1 level1 order, alphabetical within) ----
l2  <- rownames(ebmf_means)
lmap <- unique(data.frame(level2 = as.character(meta$annotation_level2),
                          level1 = as.character(meta$annotation_level1),
                          stringsAsFactors = FALSE))
if (anyDuplicated(lmap$level2)) stop("A level2 label maps to multiple level1 labels.")
lin <- lmap$level1[match(l2, lmap$level2)]
if (anyNA(lin) || !all(lin %in% level1_order)) stop("Unexpected level1 lineage among clusters.")
ord <- order(match(lin, level1_order), l2)
cluster_order <- data.frame(
  level2_cluster = l2[ord],
  level1         = lin[ord],
  display_column = seq_along(ord) - 1L,
  n_cells        = ecnt[ord],
  stringsAsFactors = FALSE
)

## ---- level1 palette (canonical lineage colors) as hex ----
l1pal <- ZemmourLib::immgent_colors$level1
hex <- vapply(l1pal, function(cc) {
  v <- grDevices::col2rgb(cc)
  grDevices::rgb(v[1], v[2], v[3], maxColorValue = 255)
}, character(1))
pal_df <- data.frame(level1 = names(hex), color = unname(hex), stringsAsFactors = FALSE)

## ---- write outputs ----
# cell -> level1/level2 for all L_pm_filtered cells; consumed by FigureS5_rematch.py
# to define common cells and clusters (avoids any machine-specific path).
fwrite(data.frame(
  cellID            = rownames(L),
  annotation_level1 = as.character(meta$annotation_level1),
  annotation_level2 = as.character(meta$annotation_level2),
  stringsAsFactors  = FALSE
), file.path(outdir, "S5_cell_metadata.csv.gz"))

fwrite(data.frame(level2_cluster = rownames(ebmf_means), ebmf_means, check.names = FALSE),
       file.path(outdir, "S5_ebmf_raw_means_level2.csv"))
fwrite(cluster_order, file.path(outdir, "S5_cluster_order.csv"))
fwrite(pal_df,        file.path(outdir, "S5_level1_palette.csv"))
fwrite(data.frame(
  metric = c("common_cells", "healthy_cells", "nonhealthy_cells", "level2_clusters", "ebmf_factors"),
  value  = c(length(common_ids),
             sum(meta$condition_broad[match(common_ids, rownames(meta))] == "healthy"),
             sum(meta$condition_broad[match(common_ids, rownames(meta))] != "healthy"),
             nlevels(grp), 200L)
), file.path(outdir, "S5_build_summary.csv"))

message("Wrote Fig S5 data inputs to ", normalizePath(outdir))
