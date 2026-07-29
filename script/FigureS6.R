# Figure S6. Protein programs, CD69-associated GPs, and further surface-protein
# gating strategies.
#
# Panels produced (see analysis/FigureS6.Rmd for the caption text):
#   s6a  Heatmap of the re-estimated protein matrix U (sparse: |score| < 0.5
#        shown white), after removing isotype/low-quality proteins and 4
#        contamination GPs (GP40/50/55/188). GPs are columns ordered from dense
#        to sparse; proteins are rows ordered by their rightmost visible GP to
#        expose a triangular support boundary.
#   s6b,s6c  Mean loading of the 10 curated CD69-associated GPs
#        (cd69_top_gps_subset, defined in code/R/citeseq_shared_setup.R) per
#        tissue (b) and per lineage (c). Figure 6d shows the same GPs' genes,
#        in the same order, which is why the subset lives in the shared setup.
#   s6d-s6g  Protein-gate vs. GP-loading comparison on the MDE embedding, for
#        the 4 curated GPs GP29/GP58/GP22/GP68.
#
# Reworked 2026-07-28. This page used to be a two-page gallery (s6-1, s6-2) of
# all 23 well-aligned GPs not shown in Figure 6; that gallery is retired. For
# anyone diffing against figures/Previous/bits/:
#
#   now   was                    what
#   s6a   Figure 6/6b            protein-program heatmap
#   s6b   Figure 6/6j            CD69 GPs, mean activity per tissue
#   s6c   Figure 6/6k            CD69 GPs, mean activity per lineage
#   s6d   inside Figure S6/s6-1  gating, GP29
#   s6e   inside Figure S6/s6-1  gating, GP58
#   s6f   inside Figure S6/s6-1  gating, GP22
#   s6g   inside Figure S6/s6-2  gating, GP68
#   --    Figure S6/s6-1, s6-2   the 2-page gallery, retired
#
# Source: ported from Figure_CITEseq.R (panels a, b, c) and
# gated_protein_loading_plot.R (panels d-g, using plot_gated_gp_vs_protein()
# from code/R/gated_protein_helpers.R, shared with Figure6.R).
#
# Required inputs (data/), read via code/R/citeseq_shared_setup.R below --
# see code/README.md's "Data provenance" table for the full picture:
#   L_pm_filtered.rds, F_pm_filtered.rds        [code/pipeline/01b_filter_cells.R]
#   igt1_96_..._ADTonly.Rds                     [primary input Seurat object]
#   protein_mat_normalized_lognorm.rds          [code/other/prepare_citeseq_protein_matrices_20260206.R]
#   umap_result.rds                             [gap, no producer script here]
#   protein_flash_selected_summary_lognorm_backfit200.rds
#     [code/other/fit_citeseq_fixed_loading_ebmf_20260206.R]
#   TableS4_citeseq_qc_20250513.csv             [external: manuscript's own Table S4]
#   Thresholds_Selected_Proteins.csv            [curated input, hand-revised; NOT regenerated
#     by code/pipeline/03_protein_thresholds.R -- see that script's header]
#   CITEseq_markers_full.rds                    [code/pipeline/04_protein_projection.R, using the
#     non-backfit200 protein summary -- see caveat above]

# --- doc:setup ---
library(ggplot2)
library(dplyr)
library(patchwork)
library(pheatmap)
library(tidyr)
library(Matrix) # protein matrices are dgCMatrix; must be attached for `[` to dispatch

data_path <- "data/"
figure_path <- "figures/final-selected/Figure S6/"
source("code/R/gated_protein_helpers.R")
source("code/R/citeseq_shared_setup.R")

# A figure script here was once seen to exit 0 with a complete log and write
# nothing at all (see script/README.md, "A re-run can silently not write"), so
# record when this run started and assert at the end that every panel is newer.
run_started_at <- Sys.time()

# ============================================================
# s6a: sparse protein-program heatmap, contamination GPs removed
# ============================================================
# The normalized protein matrix is derived here rather than in
# citeseq_shared_setup.R because this is its only consumer among the figures
# (Extended Data Table 6 rebuilds it the same way for its own export).
Protein_F_pm <- Protein_F_pm_raw[!rownames(Protein_F_pm_raw) %in% isotype_proteins, ]
Protein_F_pm <- Protein_F_pm[rownames(Protein_F_pm) %in% good_proteins, ]
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% exclude_proteins, ]
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% thy11_proteins, ]
D_lognorm <- diag(1 / apply(Protein_F_pm, 2, function(x) max(abs(x))))
Protein_F_pm <- Protein_F_pm %*% D_lognorm
colnames(Protein_F_pm) <- paste0("GP", 1:ncol(Protein_F_pm))
Protein_F_pm[is.na(Protein_F_pm)] <- 0

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

# Display proteins as rows and GPs as columns. Order GP columns from most to
# fewest visible proteins. Order protein rows by their rightmost visible GP, so
# proteins extending into the sparse right side appear first and form a
# triangular boundary. Visible count and a rarity-weighted support score provide
# deterministic secondary ordering.
wide_matrix_s6a <- as.matrix(Protein_F_pm_simplified_no_contamination)
wide_visible_mask_s6a <- abs(wide_matrix_s6a) >= sparse_cutoff
wide_gp_visible_count_s6a <- colSums(wide_visible_mask_s6a)
wide_protein_visible_count_s6a <- rowSums(wide_visible_mask_s6a)
wide_gp_number_s6a <- as.integer(sub("^GP", "", colnames(wide_matrix_s6a)))

wide_gp_order_s6a <- order(-wide_gp_visible_count_s6a, wide_gp_number_s6a)
wide_mask_ordered_cols_s6a <- wide_visible_mask_s6a[
  ,
  wide_gp_order_s6a,
  drop = FALSE
]
wide_rightmost_visible_gp_s6a <- apply(
  wide_mask_ordered_cols_s6a,
  1,
  function(values) max(which(values))
)
wide_rarity_weights_s6a <- seq_len(ncol(wide_mask_ordered_cols_s6a))^2
wide_protein_rarity_score_s6a <- as.numeric(
  wide_mask_ordered_cols_s6a %*% wide_rarity_weights_s6a
)
wide_protein_order_s6a <- order(
  -wide_rightmost_visible_gp_s6a,
  -wide_protein_visible_count_s6a,
  -wide_protein_rarity_score_s6a,
  rownames(wide_matrix_s6a)
)
wide_ordered_matrix_s6a <- wide_matrix_s6a[
  wide_protein_order_s6a,
  wide_gp_order_s6a,
  drop = FALSE
]

pdf(paste0(figure_path, "s6a.pdf"), width = 48, height = 14)
pheatmap::pheatmap(
  wide_ordered_matrix_s6a,
  main = sprintf(
    paste0(
      "Protein programs - GP columns, triangular-first protein rows ",
      "(|score| >= %.1f; protein-row sparsity is not monotone)"
    ),
    sparse_cutoff
  ),
  color = cols_sparse,
  breaks = bk_sparse,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey75",
  fontsize = 16,
  fontsize_row = 16,
  fontsize_col = 16,
  angle_col = 90,
  legend_breaks = c(-1, -sparse_cutoff, 0, sparse_cutoff, 1),
  legend_labels = c("-1", "-0.5", "0 (white)", "0.5", "1")
)
dev.off()

# ============================================================
# s6b/s6c: mean loading of the 10 curated CD69-associated GPs, per tissue (b)
# and per lineage (c). cd69_top_gps_sorted comes from citeseq_shared_setup.R
# and is the same GP order Figure 6d draws.
# ============================================================
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
p_s6b <- make_mean_loading_heatmap(df_organ, "Mean GP loading by tissue (organ_simplified)")
ggsave(paste0(figure_path, "s6b.pdf"), p_s6b, width = 9, height = 5)

df_level1 <- mean_loading_long(L_cd69_sub, meta_hm$annotation_level1, cd69_top_gps_sorted)
p_s6c <- make_mean_loading_heatmap(df_level1, "Mean GP loading by cell type (level1)")
ggsave(paste0(figure_path, "s6c.pdf"), p_s6c, width = 7, height = 5)

# ============================================================
# s6d-s6g: protein-gate vs. GP-loading comparison for the 4 curated
# supplementary GPs. Same helper, same inputs and same panel geometry as
# Figure 6e-6j -- only the GPs differ, and the two sets are disjoint.
# ============================================================
# As in Figure6.R, the loop iterates over the names of the letter map so a GP
# cannot be drawn under another GP's letter.
figs6_gating <- c("GP29" = "s6d", "GP58" = "s6e", "GP22" = "s6f", "GP68" = "s6g")
for (gp in names(figs6_gating)) {
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
    save_path = paste0(figure_path, figs6_gating[gp], ".pdf")
  )
}

# ============================================================
# Did this run actually write the panels?
# ============================================================
expected_panels <- paste0(figure_path, c("s6a", "s6b", "s6c", figs6_gating), ".pdf")
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
