# Shared setup for the Figure 3 / Figure S3 scripts (T-cell activation GPs).
#
# Figure3.R and FigureS3.R both build on the same curated set of "activation"
# GPs, the same CD4/CD8 activated-vs-resting cell groupings, and the same
# semantic color grouping (CD4-only / CD8-only / both-up / both-down). This
# was duplicated near-verbatim across panels in the original
# Figure_Activation.R; factored here so both figure scripts source it
# once instead of repeating ~150 lines of setup each.
#
# Requires L_pm_filtered, F_pm_filtered, seurat_meta_filtered to already be
# loaded (see code/R/setup_data.R::load_gp_data()).

# DGE of one group vs another, per GP, on the loading matrix L1 (using F1's
# genes for AveExpr only -- the "expression" side of a limma-style topTable).
FlashierDGE_corrected <- function(F1, L1, group1, group2, title_plot = "") {
  mean_group1 <- colMeans(L1[group1, , drop = FALSE])
  mean_group2 <- colMeans(L1[group2, , drop = FALSE])
  mean_change_loadings <- mean_group1 - mean_group2
  ave_expr <- colMeans(L1[c(group1, group2), , drop = FALSE])
  diff_factors <- data.frame(
    SYMBOL = colnames(L1),
    mean_change_loadings = mean_change_loadings,
    AveExpr = ave_expr,
    row.names = NULL
  )
  list(diff_factors = diff_factors, title = title_plot)
}

# Standardized mean difference (activated vs resting) per GP.
# d_k = (mean_act - mean_rest) / sd_pooled
std_mean_diff <- function(mat, idx_a, idx_b, sd_ref) {
  d <- (colMeans(mat[idx_a, , drop = FALSE]) -
    colMeans(mat[idx_b, , drop = FALSE])) /
    sd_ref
  d[!is.finite(d)] <- 0
  d
}

# ============================================================
# CD4/CD8 resting vs activated cell groups + per-GP DGE
# ============================================================
CD4_cells <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 == "CD4"]
CD4_resting_cells <- CD4_cells[
  seurat_meta_filtered$annotation_level2_group[match(CD4_cells, seurat_meta_filtered$cellID)] == "resting"
]
CD4_activated_cells <- CD4_cells[
  seurat_meta_filtered$annotation_level2_group[match(CD4_cells, seurat_meta_filtered$cellID)] == "activated"
]
CD4_DGE <- FlashierDGE_corrected(F1 = F_pm_filtered, L1 = L_pm_filtered, group1 = CD4_activated_cells, group2 = CD4_resting_cells)

CD8_cells <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 == "CD8"]
CD8_resting_cells <- CD8_cells[
  seurat_meta_filtered$annotation_level2_group[match(CD8_cells, seurat_meta_filtered$cellID)] == "resting"
]
CD8_activated_cells <- CD8_cells[
  seurat_meta_filtered$annotation_level2_group[match(CD8_cells, seurat_meta_filtered$cellID)] == "activated"
]
CD8_DGE <- FlashierDGE_corrected(F1 = F_pm_filtered, L1 = L_pm_filtered, group1 = CD8_activated_cells, group2 = CD8_resting_cells)

diff_factors_CD4 <- CD4_DGE$diff_factors %>% dplyr::rename(mean_change_loadings_CD4 = mean_change_loadings, AveExpr_CD4 = AveExpr)
diff_factors_CD8 <- CD8_DGE$diff_factors %>% dplyr::rename(mean_change_loadings_CD8 = mean_change_loadings, AveExpr_CD8 = AveExpr)
diff_factors_merged <- diff_factors_CD4 %>% dplyr::inner_join(diff_factors_CD8, by = "SYMBOL")

# ============================================================
# Curated activation-GP set + semantic color grouping
#   Blue = CD4 only, Orange = CD8 only, Darkred = both up, Darkgreen = both down
# ============================================================
GPs_of_interest <- paste0("GP", c(
  25, 26, 10, 12, 58, 171, 9, 79, 35, 177, 161, 152, 162, 36, 56, 181, 32, 80,
  49, 41, 57, 176, 11, 13, 159
))
highlight_colors <- c(
  "GP56" = "blue", "GP162" = "blue", "GP36" = "blue", "GP152" = "blue",
  "GP161" = "blue", "GP177" = "blue", "GP79" = "blue", "GP12" = "blue",
  "GP13" = "blue", "GP159" = "blue",
  "GP10" = "darkorange2", "GP58" = "darkorange2", "GP181" = "darkorange2", "GP176" = "darkorange2",
  "GP25" = "darkred", "GP26" = "darkred", "GP35" = "darkred", "GP32" = "darkred",
  "GP80" = "darkred", "GP57" = "darkred",
  "GP9" = "darkgreen", "GP171" = "darkgreen", "GP49" = "darkgreen", "GP41" = "darkgreen", "GP11" = "darkgreen"
)

# Shared per-GP denominator: SD over all activated+resting CD4 and CD8 cells.
sd_pooled_act_rest <- apply(
  L_pm_filtered[c(CD4_activated_cells, CD4_resting_cells, CD8_activated_cells, CD8_resting_cells), , drop = FALSE],
  2, sd
)
d_CD4 <- std_mean_diff(L_pm_filtered, CD4_activated_cells, CD4_resting_cells, sd_pooled_act_rest)
d_CD8 <- std_mean_diff(L_pm_filtered, CD8_activated_cells, CD8_resting_cells, sd_pooled_act_rest)
d_factors_merged <- data.frame(SYMBOL = colnames(L_pm_filtered), d_CD4 = d_CD4, d_CD8 = d_CD8)

# Explicit GP display order used by the heatmaps (S3F, Fig3d/e): grouped by
# semantic color (blue / orange / darkred / darkgreen blocks).
ordered_GPs <- c(
  "GP56", "GP162", "GP36", "GP152", "GP161", "GP177", "GP79", "GP12", "GP13", "GP159", # "Blue" group
  "GP10", "GP58", "GP181", "GP176", # "Orange" group
  "GP25", "GP26", "GP35", "GP32", "GP80", "GP57", # "Darkred" group (both up)
  "GP9", "GP171", "GP49", "GP41", "GP11" # "Green" group
)

# ============================================================
# Normalized gene-score matrix (max |score| = 1 per GP) -- used by the
# GP-gene network (Fig3b), the TF-GP networks (Fig3c, S3E), and reused
# wherever a per-GP-normalized score is needed.
# ============================================================
F_pm_filtered_norm <- scale_cols(F_pm_filtered, 1 / apply(abs(F_pm_filtered), 2, max))
colnames(F_pm_filtered_norm) <- paste0("GP", seq_len(ncol(F_pm_filtered_norm)))

# ============================================================
# GP group colors / membership (derived from highlight_colors), and the
# curated-GP loading submatrix -- used by Fig3d/e and S3F.
# ============================================================
group_colors <- c("CD4 only" = "blue", "CD8 only" = "darkorange2", "both up" = "darkred", "both down" = "darkgreen")
color_to_group <- setNames(names(group_colors), unname(group_colors))
gp_to_group <- setNames(color_to_group[unname(highlight_colors)], names(highlight_colors))
gp_groups <- split(names(gp_to_group), factor(gp_to_group, levels = names(group_colors)))
gp_row_order <- unlist(gp_groups, use.names = FALSE)
lineages <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN")
L_subset <- L_pm_filtered[, GPs_of_interest, drop = FALSE]
