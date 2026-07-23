# Figure S5 (data step). EBMF vs matched-RQVI level2-cluster comparison.
#
# Design (version B, Tianze-style plotting done in script/FigureS5_plot.py):
#   * Cells: OUR L_pm_filtered cells (passed the iterative total-loading filter),
#     restricted to non-thymocytes (annotation_level1 != "thymocyte"), ALL
#     conditions (no healthy-only restriction), intersected with the cells that
#     have RQVI loadings -> "common cells".
#   * EBMF matrix: our flashier loadings L_pm_filtered (GP1..GP200), averaged
#     within annotation_level2 on the common cells.
#   * RQVI matrix: Tianze's 200 matched RQVI programs (raw cell loadings),
#     averaged within annotation_level2 on the SAME common cells. Each matched
#     program is placed under the column of its paired EBMF factor. The pairing
#     is Tianze's one-to-one match table; F_k == our GP_k (verified at cell level,
#     cell-level Pearson r = 1.0 for all 200).
#   * Columns (level2 clusters) are ordered by the Figure-1 level1 lineage order
#     and alphabetically within each lineage.
#
# This script writes raw cluster-mean matrices + column/palette metadata.
# Row ordering (hierarchical clustering of EBMF), per-factor 0-1 scaling, and
# the heatmaps are produced by script/FigureS5_plot.py in Tianze's visual style.

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
parquet_path <- file.path(PKG, "rqvi_matched_200_cell_loadings.parquet")
matches_path <- file.path(PKG, "ebmf_rqvi_multiseed_level2_one_to_one_matches.csv")
level1_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN", "DP")

## ---- our EBMF loadings + metadata ----
gp <- load_gp_data()
L  <- gp$L_pm_filtered
colnames(L) <- paste0("GP", seq_len(ncol(L)))
meta <- gp$seurat_meta_filtered
stopifnot(identical(rownames(L), rownames(meta)), ncol(L) == 200L)

nonthy <- meta$annotation_level1 != "thymocyte"
if (anyNA(nonthy)) stop("annotation_level1 has missing values.")

## ---- matched RQVI cell loadings (200 programs) ----
pq   <- as.data.frame(arrow::read_parquet(parquet_path))
prog <- setdiff(colnames(pq), "cell_id")
stopifnot(length(prog) == 200L)
matches <- fread(matches_path)
f_for <- matches$ebmf_factor[match(prog, matches$rqvi_candidate)]   # "F<k>" per program
if (anyNA(f_for) || length(unique(f_for)) != 200L) {
  stop("Could not map every RQVI program column to a unique EBMF factor.")
}

## ---- common cells (non-thymocyte, in L_pm, and with RQVI loading) ----
our_nonthy_ids <- rownames(L)[nonthy]
common_ids <- our_nonthy_ids[our_nonthy_ids %in% pq$cell_id]
grp <- factor(as.character(meta$annotation_level2[match(common_ids, rownames(meta))]))
grp <- droplevels(grp)
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

## ---- matched-RQVI cluster means (columns F1..F200, same pairing) ----
pqc <- as.matrix(pq[match(common_ids, pq$cell_id), prog, drop = FALSE])
colnames(pqc) <- f_for
pqc <- pqc[, paste0("F", seq_len(200)), drop = FALSE]     # reorder to F1..F200
rsum <- rowsum(pqc, grp)
rcnt <- as.integer(table(grp)[rownames(rsum)])
rqvi_means <- sweep(rsum, 1L, rcnt, "/")

stopifnot(identical(rownames(ebmf_means), rownames(rqvi_means)),
          identical(ecnt, rcnt),                           # same cells -> same counts
          nrow(ebmf_means) == 107L, ncol(ebmf_means) == 200L)

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

## ---- level1 palette (our canonical colors) as hex ----
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
fwrite(data.frame(level2_cluster = rownames(rqvi_means), rqvi_means, check.names = FALSE),
       file.path(outdir, "S5_rqvi_matched_raw_means_level2.csv"))
fwrite(cluster_order, file.path(outdir, "S5_cluster_order.csv"))
fwrite(pal_df,        file.path(outdir, "S5_level1_palette.csv"))
fwrite(data.frame(
  metric = c("common_cells", "healthy_cells", "nonhealthy_cells", "level2_clusters",
             "ebmf_factors", "rqvi_matched_programs"),
  value  = c(length(common_ids),
             sum(meta$condition_broad[match(common_ids, rownames(meta))] == "healthy"),
             sum(meta$condition_broad[match(common_ids, rownames(meta))] != "healthy"),
             nlevels(grp), 200L, 200L)
), file.path(outdir, "S5_build_summary.csv"))

message("Wrote Fig S5 data inputs to ", normalizePath(outdir))
