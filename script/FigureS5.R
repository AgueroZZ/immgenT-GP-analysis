# Figure S5. Protein programs across GPs.
#
# One panel (see analysis/FigureS5.Rmd for the caption text):
#   s5   Heatmap of the re-estimated protein matrix U (sparse: |score| < 0.5
#        shown white), after removing isotype/low-quality proteins and 4
#        contamination GPs (GP40/50/55/188). GPs are columns ordered from dense
#        to sparse; proteins are rows ordered by their rightmost visible GP to
#        expose a triangular support boundary.
#
# Split out of Figure S6 on 2026-07-30. It was that figure's panel s6a (and
# before 2026-07-28, the published Figure 6b); it is now a standalone figure and
# takes the Extended Data Figure 5 slot that main Figure 7b vacated. Figure S6's
# remaining panels each dropped one letter (its s6b-s6g are now s6a-s6f). For
# anyone diffing against figures/Previous/bits/:
#
#   now   was                    what
#   s5    Figure 6/6b            protein-program heatmap
#
# Source: ported from Figure_CITEseq.R (its panel a).
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
#
# The shared setup loads more than this one panel needs (it also builds the
# gating inputs Figure 6 / Figure S6 use); it is sourced as-is rather than
# trimmed so the protein filters here cannot drift from the ones those figures
# and Extended Data Table 6 apply.

# --- doc:setup ---
library(dplyr)
library(pheatmap)
library(Matrix) # protein matrices are dgCMatrix; must be attached for `[` to dispatch

data_path <- "data/"
figure_path <- "figures/final-selected/Figure S5/"
source("code/R/citeseq_shared_setup.R")

# A figure script here was once seen to exit 0 with a complete log and write
# nothing at all (see script/README.md, "A re-run can silently not write"), so
# record when this run started and assert at the end that the panel is newer.
run_started_at <- Sys.time()

# ============================================================
# s5: sparse protein-program heatmap, contamination GPs removed
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
wide_matrix_s5 <- as.matrix(Protein_F_pm_simplified_no_contamination)
wide_visible_mask_s5 <- abs(wide_matrix_s5) >= sparse_cutoff
wide_gp_visible_count_s5 <- colSums(wide_visible_mask_s5)
wide_protein_visible_count_s5 <- rowSums(wide_visible_mask_s5)
wide_gp_number_s5 <- as.integer(sub("^GP", "", colnames(wide_matrix_s5)))

wide_gp_order_s5 <- order(-wide_gp_visible_count_s5, wide_gp_number_s5)
wide_mask_ordered_cols_s5 <- wide_visible_mask_s5[
  ,
  wide_gp_order_s5,
  drop = FALSE
]
wide_rightmost_visible_gp_s5 <- apply(
  wide_mask_ordered_cols_s5,
  1,
  function(values) max(which(values))
)
wide_rarity_weights_s5 <- seq_len(ncol(wide_mask_ordered_cols_s5))^2
wide_protein_rarity_score_s5 <- as.numeric(
  wide_mask_ordered_cols_s5 %*% wide_rarity_weights_s5
)
wide_protein_order_s5 <- order(
  -wide_rightmost_visible_gp_s5,
  -wide_protein_visible_count_s5,
  -wide_protein_rarity_score_s5,
  rownames(wide_matrix_s5)
)
wide_ordered_matrix_s5 <- wide_matrix_s5[
  wide_protein_order_s5,
  wide_gp_order_s5,
  drop = FALSE
]

pdf(paste0(figure_path, "s5.pdf"), width = 48, height = 14)
pheatmap::pheatmap(
  wide_ordered_matrix_s5,
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
# Did this run actually write the panel?
# ============================================================
expected_panels <- paste0(figure_path, "s5", ".pdf")
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
