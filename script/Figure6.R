# Figure 6. Linking GPs to proteins.
#
# Panels produced (see figures/Figure6_CITEseq/Figure6_caption.md for the
# full caption text):
#   6a  Schematic of the projection (EBMF on scRNA, then EBMF on CITE-seq
#       with cell loadings fixed). Hand-drawn in other software -- NOT
#       code-generated, no source to port. Not produced by this script.
#   6b  Heatmap of the re-estimated protein matrix U (sparse: |score| < 0.5
#       shown white), after removing isotype/low-quality proteins and 4
#       contamination GPs (GP40/50/55/188).
#   6c-6f  Protein-gate vs. GP-loading comparison on the MDE embedding, for
#       the 4 curated GPs GP171/GP23/GP12/GP80.
#   6g,6h  KLRG1 modulation: CD8 vs CD4 (g) and CD8 vs Treg (h).
#   6i  Heatmap of top-5 gene scores per GP for the 10 GPs most associated
#       with CD69, with a CD69-correlation strip.
#   6j,6k  Mean loading of those 10 GPs per tissue (j) and per lineage (k).
#
# Source: ported from Figure_CITEseq.R (panels b, g, h, i, j, k) and
# gated_protein_loading_plot.R (panels c-f, using
# plot_gated_gp_vs_protein() from code/R/gated_protein_helpers.R,
# shared with FigureS6.R).
#
# Required inputs (data/), read via code/R/citeseq_shared_setup.R below --
# see code/README.md's "Data provenance" table for the full picture:
#   L_pm_filtered.rds, F_pm_filtered.rds        [code/pipeline/01b_filter_cells.R]
#   igt1_96_..._ADTonly.Rds                     [primary input Seurat object]
#   protein_mat_normalized_lognorm.rds          [gap, no producer script here]
#   umap_result.rds                             [gap, no producer script here]
#   protein_flash_selected_summary_lognorm_backfit200.rds
#     [gap -- code/pipeline/04_protein_projection.R produces the
#     non-backfit200 variant, not this exact file]
#   TableS4_citeseq_qc_20250513.csv             [external: manuscript's own Table S4]
#   Thresholds_Selected_Proteins.csv            [code/pipeline/03_protein_thresholds.R]
#   CITEseq_markers_full.rds                    [code/pipeline/04_protein_projection.R, using the
#     non-backfit200 protein summary -- see caveat above]

library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)
library(pheatmap)
library(tidyr)
library(Matrix) # protein matrices are dgCMatrix; must be attached for `[` to dispatch

data_path <- "data/"
figure_path <- "figures/generated/Figure 6/"
source("code/R/gated_protein_helpers.R")

# 6a: hand-drawn schematic -- not code-generated, no output here.

# ============================================================
# Load data (shared with FigureS6.R)
# ============================================================
source("code/R/citeseq_shared_setup.R")

# Re-derive the normalized protein matrix used only by this script's panel
# 6b (Protein_F_pm_raw filtered/scaled -- FigureS6.R doesn't need it).
Protein_F_pm <- Protein_F_pm_raw[!rownames(Protein_F_pm_raw) %in% isotype_proteins, ]
Protein_F_pm <- Protein_F_pm[rownames(Protein_F_pm) %in% good_proteins, ]
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% exclude_proteins, ]
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% thy11_proteins, ]
D_lognorm <- diag(1 / apply(Protein_F_pm, 2, function(x) max(abs(x))))
Protein_F_pm <- Protein_F_pm %*% D_lognorm
colnames(Protein_F_pm) <- paste0("GP", 1:ncol(Protein_F_pm))
Protein_F_pm[is.na(Protein_F_pm)] <- 0

# ============================================================
# 6b: sparse protein-program heatmap, contamination GPs removed
# ============================================================
threshold_simplified <- 0
keep_rows_simplified <- apply(Protein_F_pm, 1, function(v) any(abs(v) > threshold_simplified, na.rm = TRUE))
Protein_F_pm_simplified <- Protein_F_pm[keep_rows_simplified, , drop = FALSE]
keep_cols_simplified <- apply(Protein_F_pm_simplified, 2, function(v) any(abs(v) > threshold_simplified, na.rm = TRUE))
Protein_F_pm_simplified <- Protein_F_pm_simplified[, keep_cols_simplified, drop = FALSE]

GP_contamination <- c("GP40", "GP50", "GP55", "GP188")
Protein_F_pm_simplified_no_contamination <- Protein_F_pm_simplified[, !colnames(Protein_F_pm_simplified) %in% GP_contamination, drop = FALSE]

sparse_cutoff <- 0.5
bk_sparse <- unique(c(seq(-1, -sparse_cutoff, length.out = 26), seq(-sparse_cutoff, sparse_cutoff, length.out = 51), seq(sparse_cutoff, 1, length.out = 26)))
cols_sparse <- c(colorRampPalette(c("#4575B4", "white"))(25), rep("white", 50), colorRampPalette(c("white", "#D73027"))(25))

pdf(paste0(figure_path, "6b.pdf"), width = 20, height = 40)
pheatmap::pheatmap(
  t(Protein_F_pm_simplified_no_contamination),
  main = sprintf("Protein programs (sparse, no contamination GPs): %d proteins x %d GPs", nrow(Protein_F_pm_simplified_no_contamination), ncol(Protein_F_pm_simplified_no_contamination)),
  color = cols_sparse,
  breaks = bk_sparse,
  border_color = "black"
)
dev.off()

# ============================================================
# 6c-6f: protein-gate vs. GP-loading comparison for the 4 curated main-figure GPs
# (df_markers2, thymocyte/proliferating/miniverse_cells, L_pm_for_gating,
# select_proteins, threshold_results_subset_manual all come from
# citeseq_shared_setup.R above)
# ============================================================
GPs_fig6 <- c("GP171", "GP23", "GP12", "GP80")
fig6_letter <- c("GP171" = "6c", "GP23" = "6d", "GP12" = "6e", "GP80" = "6f")
enlarge_gps <- c("GP8", "GP30", "GP170", "GP107")
for (gp in GPs_fig6) {
  k_name <- paste0("K", sub("^GP", "", gp))
  plot_gated_gp_vs_protein(
    gp_name = k_name,
    df_markers = df_markers2,
    protein_mat = protein_mat_normalized_lognorm,
    loading_mat = L_pm_for_gating,
    mde_emb = mde_result,
    missing_threshold_action = "skip",
    threshold_df = threshold_results_subset_manual,
    exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells),
    selected_proteins = select_proteins,
    loading_q = NULL,
    min_pointsize = if (gp %in% enlarge_gps) 3L else 0L,
    save_path = paste0(figure_path, fig6_letter[gp], ".pdf")
  )
}

# ============================================================
# 6g/6h: KLRG1 modulation (CD8 vs CD4, CD8 vs Treg)
# ============================================================
FlashierDGE_corrected <- function(F1, L1, group1, group2, title_plot = "") {
  loadings_group1 <- colMeans(L1[group1, ])
  loadings_group2 <- colMeans(L1[group2, ])
  mean_change_loadings <- loadings_group1 - loadings_group2
  vplot <- data.frame(SYMBOL = names(mean_change_loadings), mean_change_loadings = mean_change_loadings, AveExpr = colMeans(L1[c(group1, group2), ]))
  list(diff_factors = vplot)
}
get_klrg1_split <- function(cell_type_label, meta, protein_data, threshold) {
  cells <- meta$cellID[meta$annotation_level1 == cell_type_label]
  cells <- intersect(cells, rownames(protein_data))
  list(pos = cells[protein_data[cells, "KLRG1"] >= threshold], neg = cells[protein_data[cells, "KLRG1"] < threshold])
}
run_checked_dge <- function(group_list, F_mat, L_mat, label) {
  if (length(group_list$pos) < 3 || length(group_list$neg) < 3) stop(paste("Insufficient data:", label))
  df <- FlashierDGE_corrected(F1 = F_mat, L1 = L_mat, group1 = group_list$pos, group2 = group_list$neg)$diff_factors
  if (!"SYMBOL" %in% colnames(df)) df$SYMBOL <- rownames(df)
  df
}
plot_target_gps <- function(df, x_var, y_var, label_var, target_gps, highlight_color = "darkorange", background_color = "black",
                             x_limits = c(-0.5, 0.5), y_limits = c(-0.5, 0.5), background_alpha = 0.5,
                             xlab = "Difference in Mean Loading", ylab = "Difference in Mean Loading", title = "Comparison of Specific GP Loadings") {
  highlight_df <- df %>% filter({{ label_var }} %in% target_gps) %>% mutate(.label_display = as.character({{ label_var }}))
  ggplot(df, aes(x = {{ x_var }}, y = {{ y_var }})) +
    geom_point(color = background_color, alpha = background_alpha) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
    geom_point(data = highlight_df, aes(color = {{ label_var }}), size = 2) +
    ggrepel::geom_text_repel(data = highlight_df, aes(label = .label_display, color = {{ label_var }}), max.overlaps = Inf, size = 3.5, box.padding = 0.35, point.padding = 0.5, segment.color = "grey50", show.legend = FALSE) +
    scale_color_manual(values = highlight_color, guide = "none") +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal()
}

klrg1_threshold <- threshold_results_subset_manual$Threshold[threshold_results_subset_manual$Protein == "KLRG1"]
cd8_split <- get_klrg1_split("CD8", seurat_meta_filtered, protein_mat_normalized_lognorm, klrg1_threshold)
CD4_split <- get_klrg1_split("CD4", seurat_meta_filtered, protein_mat_normalized_lognorm, klrg1_threshold)
treg_split <- get_klrg1_split("Treg", seurat_meta_filtered, protein_mat_normalized_lognorm, klrg1_threshold)

diff_CD8 <- run_checked_dge(cd8_split, F_pm_filtered, L_pm_filtered, "CD8") %>% rename(mean_change_CD8 = mean_change_loadings, AveExpr_CD8 = AveExpr)
diff_CD4 <- run_checked_dge(CD4_split, F_pm_filtered, L_pm_filtered, "CD4") %>% rename(mean_change_CD4 = mean_change_loadings, AveExpr_CD4 = AveExpr)
diff_Treg <- run_checked_dge(treg_split, F_pm_filtered, L_pm_filtered, "Treg") %>% rename(mean_change_Treg = mean_change_loadings, AveExpr_Treg = AveExpr)

# 6g: CD8 vs CD4
merged_cd4 <- inner_join(diff_CD4, diff_CD8, by = "SYMBOL")
p_6g <- plot_target_gps(
  df = merged_cd4, x_var = mean_change_CD8, y_var = mean_change_CD4, label_var = SYMBOL,
  target_gps = c("GP10", "GP58", "GP25", "GP26", "GP43"), background_alpha = 0.8, x_limits = c(-0.2, 0.4), y_limits = c(-0.2, 0.4),
  highlight_color = c("GP10" = "darkorange2", "GP25" = "blue", "GP43" = "blue", "GP26" = "blue", "GP58" = "darkorange2"),
  title = "KLRG1 Modulation: CD8 vs CD4", xlab = "Effect Size in CD8 (KLRG1+ - KLRG1-)", ylab = "Effect Size in CD4 (KLRG1+ - KLRG1-)"
) + theme_bw()
ggsave(paste0(figure_path, "6g.pdf"), p_6g, width = 7, height = 6)

# 6h: CD8 vs Treg
merged_treg <- inner_join(diff_Treg, diff_CD8, by = "SYMBOL")
p_6h <- plot_target_gps(
  df = merged_treg, x_var = mean_change_CD8, y_var = mean_change_Treg, label_var = SYMBOL,
  target_gps = c("GP6", "GP10", "GP12", "GP27", "GP68", "GP58"), background_alpha = 0.8, x_limits = c(-0.2, 0.4), y_limits = c(-0.2, 0.4),
  highlight_color = c("GP10" = "darkorange2", "GP27" = "deeppink", "GP6" = "deeppink", "GP68" = "deeppink", "GP12" = "deeppink", "GP58" = "darkorange2"),
  title = "KLRG1 Modulation: CD8 vs Treg", xlab = "Effect Size in CD8 (KLRG1+ - KLRG1-)", ylab = "Effect Size in Treg (KLRG1+ - KLRG1-)"
) + theme_bw()
ggsave(paste0(figure_path, "6h.pdf"), p_6h, width = 7, height = 6)

# ============================================================
# 6i/6j/6k: the 10 GPs most associated with CD69
# ============================================================
D_scale6 <- diag(1 / apply(F_pm_filtered, 2, function(x) max(abs(x), na.rm = TRUE)))
F_pm_filtered_scaled <- F_pm_filtered %*% D_scale6
colnames(F_pm_filtered_scaled) <- paste0("GP", 1:ncol(F_pm_filtered_scaled))

cd69_top_gps_subset <- c("GP35", "GP6", "GP170", "GP26", "GP58", "GP171", "GP63", "GP62", "GP3", "GP29")
shared_cells_cd69 <- intersect(rownames(L_pm_filtered), rownames(protein_mat_normalized_lognorm))
cd69_expr_vec <- protein_mat_normalized_lognorm[shared_cells_cd69, "CD69"]
cd69_corr <- sapply(cd69_top_gps_subset, function(gp) cor(L_pm_filtered[shared_cells_cd69, gp], cd69_expr_vec, method = "spearman"))
cd69_top_gps_sorted <- names(sort(cd69_corr, decreasing = FALSE)) # most-correlated GP ends up at top of y-axis

plot_factor_heatmap <- function(F_matrix, gp_vector, n_top = 5, min_abs_loading = 0.5, transpose = FALSE,
                                 title = "Factor loadings", low_color = "steelblue", mid_color = "white", high_color = "firebrick", font_size = 9) {
  F_sub <- F_matrix[, gp_vector, drop = FALSE]
  selected_genes <- lapply(gp_vector, function(gp) {
    vals <- F_sub[, gp]
    top_pos <- names(sort(vals, decreasing = TRUE))[seq_len(min(n_top, sum(vals > 0)))]
    top_neg <- names(sort(vals, decreasing = FALSE))[seq_len(min(n_top, sum(vals < 0)))]
    c(top_pos, top_neg)
  })
  selected_genes <- unique(unlist(selected_genes))
  if (min_abs_loading > 0) {
    max_abs <- apply(F_sub[selected_genes, , drop = FALSE], 1, function(x) max(abs(x), na.rm = TRUE))
    selected_genes <- names(max_abs[max_abs >= min_abs_loading])
  }
  hc_genes <- hclust(dist(F_sub[selected_genes, , drop = FALSE]))
  gene_order <- rownames(F_sub[selected_genes, , drop = FALSE])[hc_genes$order]
  plot_df <- F_sub[selected_genes, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene") %>%
    tidyr::pivot_longer(cols = -Gene, names_to = "GP", values_to = "Loading") %>%
    mutate(GP = factor(GP, levels = gp_vector), Gene = factor(Gene, levels = gene_order))
  limit <- max(abs(plot_df$Loading), na.rm = TRUE)
  x_aes <- if (transpose) "Gene" else "GP"
  y_aes <- if (transpose) "GP" else "Gene"
  ggplot(plot_df, aes(x = .data[[x_aes]], y = .data[[y_aes]], fill = Loading)) +
    geom_tile() +
    scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, midpoint = 0, limits = c(-limit, limit), name = "Loading") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = font_size) +
    theme(axis.text.x = element_text(angle = if (transpose) 90 else 45, hjust = 1, size = font_size), axis.text.y = element_text(size = font_size), panel.grid = element_blank())
}

p_heatmap <- plot_factor_heatmap(F_matrix = F_pm_filtered_scaled, gp_vector = cd69_top_gps_sorted, n_top = 5, font_size = 9, transpose = TRUE, low_color = "#4DAF4A", mid_color = "white", high_color = "#984EA3")

corr_strip_df <- data.frame(GP = factor(cd69_top_gps_sorted, levels = cd69_top_gps_sorted), Correlation = cd69_corr[cd69_top_gps_sorted], x = "Corr")
corr_limit <- max(abs(corr_strip_df$Correlation))
p_corr_strip <- ggplot(corr_strip_df, aes(x = x, y = GP, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "tomato", midpoint = 0, limits = c(-corr_limit, corr_limit), name = "Corr\n(CD69)") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())

p_6i <- p_corr_strip + p_heatmap + patchwork::plot_layout(widths = c(0.06, 1), guides = "collect")
ggsave(paste0(figure_path, "6i.pdf"), p_6i, width = 11, height = 5)

# 6j/6k: mean loading of these GPs per tissue (j) and per lineage (k)
cells_for_heatmap <- intersect(rownames(L_pm_filtered), rownames(seurat_meta_filtered))
L_cd69_sub <- L_pm_filtered[cells_for_heatmap, cd69_top_gps_sorted, drop = FALSE]
meta_hm <- seurat_meta_filtered[cells_for_heatmap, c("annotation_level1", "organ_simplified")]

mean_loading_long <- function(L_mat, group_vec, gp_levels) {
  as.data.frame(L_mat) %>%
    mutate(group = group_vec) %>%
    tidyr::pivot_longer(cols = -group, names_to = "GP", values_to = "Loading") %>%
    group_by(group, GP) %>%
    summarise(mean_loading = mean(Loading, na.rm = TRUE), .groups = "drop") %>%
    mutate(GP = factor(GP, levels = gp_levels))
}
make_mean_loading_heatmap <- function(df, title) {
  fill_max <- max(df$mean_loading, na.rm = TRUE)
  ggplot(df, aes(x = group, y = GP, fill = mean_loading)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, fill_max), name = "Mean\nloading") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), axis.text.y = element_text(size = 9), panel.grid = element_blank())
}

df_organ <- mean_loading_long(L_cd69_sub, meta_hm$organ_simplified, cd69_top_gps_sorted)
p_6j <- make_mean_loading_heatmap(df_organ, "Mean GP loading by tissue (organ_simplified)")
ggsave(paste0(figure_path, "6j.pdf"), p_6j, width = 9, height = 5)

df_level1 <- mean_loading_long(L_cd69_sub, meta_hm$annotation_level1, cd69_top_gps_sorted)
p_6k <- make_mean_loading_heatmap(df_level1, "Mean GP loading by cell type (level1)")
ggsave(paste0(figure_path, "6k.pdf"), p_6k, width = 7, height = 5)
