# Figure S1. GP reproducibility across IGTs.
#
# Panels produced (see figures/Previous/bits/Figure S1/FigureS1_caption.md
# for the full caption text):
#   S1A  Cumulative number of GPs validated (cosine >= threshold, thresholds
#        0.2-0.8) as IGTs are added one at a time, in IGT index order.
#   S1B  Number of GPs validated by at least X IGTs, vs X (log-log), for the
#        same thresholds.
#   S1C  Between-IGT variability: each GP's mean-of-per-IGT-mean-loading (x)
#        vs. variance-of-per-IGT-mean-loading (y).
#   S1D  Heatmap of per-IGT mean loading for S1C's 10 highest-variance GPs.
#
#        S1C and S1D share one per-IGT mean-loading matrix, over the
#        standard-spleen subset restricted to the 18 IGTs with >= 500 such
#        cells. Figure_batch.R built it twice on different IGT sets (S1C over
#        all 35, S1D over the 18), which made the published panels disagree --
#        see the note in the "Load data for S1C/S1D" section.
#   S1E  Scatter of the NUMBER of active genes (x) vs. proportion of active
#        cells (y) per GP, using the same hard-threshold definitions as Figure 2
#        (|normalized score| > 0.25 for genes; normalized loading > 0.1 for
#        cells), over non-thymocyte cells -- not the EBMF sparsity prior. One
#        dot per GP.
#
# Source: S1A/S1B ported from Figure_Saturation.R; S1C/S1D from
# Figure_batch.R (panels a/b only -- that script's `plot_gp_loading()`
# helper is defined but never called for a saved output, so it's dropped).
#
# S1A/S1B reuse the per-IGT cosine-matching score matrix
# (data/igt_specific_cosine_scores.csv) rather than recomputing it here --
# recomputing requires Hungarian-matching each of the ~80 per-IGT
# refactorizations in data/igt_specific/*.qs against the full model, which is
# the job of code/pipeline/05_igt_validation.R (run once upstream).
#
# Required inputs (data/) -- see code/README.md's "Data provenance" table
# for the full picture:
#   igt_specific_cosine_scores.csv           [code/pipeline/05_igt_validation.R]
#   L_pm_filtered.rds                        [code/pipeline/01b_filter_cells.R]
#   igt1_96_..._ADTonly.Rds                  [primary input Seurat object]

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

# --- doc:setup-ab ---
data_path <- "data/"
figure_path <- "figures/final-selected/Figure S1/"
gp_label <- function(x) sub("^K(\\d+)$", "GP\\1", x)

# ============================================================
# S1A/S1B: load the cached per-IGT cosine score matrix
# (GPs x IGTs; produced by code/pipeline/05_igt_validation.R)
# ============================================================
score_mat <- as.matrix(read.csv(paste0(data_path, "igt_specific_cosine_scores.csv"), row.names = 1, check.names = FALSE))

# ============================================================
# S1B: number of GPs validated by at least X IGTs, vs X (log-log)
# ============================================================
thresholds <- seq(0.2, 0.8, by = 0.1)
X_grid <- 1:50
plot_df_b <- tidyr::crossing(threshold = thresholds, X = X_grid) %>%
  mutate(n_GP = purrr::map2_int(threshold, X, \(t, x) {
    rowSums(score_mat >= t, na.rm = TRUE) |> (\(v) sum(v >= x))()
  }))

p_S1B <- ggplot(plot_df_b, aes(x = X, y = n_GP, group = factor(threshold))) +
  geom_line() +
  geom_point(size = 1) +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "X (validated by at least X IGTs)", y = "Number of GPs", color = "Threshold") +
  aes(color = factor(threshold)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")
ggsave(paste0(figure_path, "S1B.pdf"), plot = p_S1B, width = 6, height = 4)

# ============================================================
# S1A: cumulative number of GPs validated as IGTs are added, in IGT-index order
# ============================================================
igt_idx <- as.integer(gsub("^IGT", "", colnames(score_mat)))
o <- order(igt_idx)
score_mat_ord <- score_mat[, o, drop = FALSE]

cum_validated_counts <- function(score_mat_ord, threshold) {
  validated <- score_mat_ord >= threshold
  ever_validated <- t(apply(validated, 1, cummax)) # 200 x nIGT logical
  colSums(ever_validated)
}
plot_df_a <- lapply(thresholds, function(t) {
  y <- cum_validated_counts(score_mat_ord, t)
  data.frame(n_IGTs_included = seq_along(y), validated_GPs = y, threshold = factor(t))
}) %>% bind_rows()

p_S1A <- ggplot(plot_df_a, aes(x = n_IGTs_included, y = validated_GPs, color = threshold)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  labs(x = "Number of IGTs included (in IGT index order)", y = "Number of validated GPs (cumulative union)", color = "Threshold") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0, ncol(score_mat_ord), by = 5)) +
  scale_y_continuous(breaks = seq(0, max(plot_df_a$validated_GPs), by = 20))
ggsave(paste0(figure_path, "S1A.pdf"), plot = p_S1A, width = 6, height = 4)

# ============================================================
# Load data for S1C/S1D
# ============================================================
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
seurat_meta_filtered_spleen <- seurat_meta_filtered %>% filter(spleen_standard == TRUE)

# S1C and S1D are both computed from ONE per-IGT mean-loading matrix
# (IGTs x GPs) over the standard-spleen subset, restricted to the IGTs with
# >= 500 such cells -- an IGT mean over a handful of cells is too noisy to
# either rank on or draw. Both panels therefore describe the same 18 IGTs, and
# S1C's ten labelled GPs are exactly S1D's ten rows.
#
# Figure_batch.R built this matrix twice instead: S1C over all 35 spleen IGTs,
# S1D over the >= 500-cell subset. The two rankings disagree -- GP2 is 7th on
# the subset but 11th over all 35, GP25 is 6th over all 35 but 23rd on the
# subset -- so the published S1C labels GP25 while the published S1D shows GP2,
# and "the top ten from (C)" could not be followed across the two panels.
# Computing it once removes the inconsistency and the chance of the two
# drifting apart again.
spleen_cells <- intersect(rownames(L_pm_filtered), rownames(seurat_meta_filtered_spleen))
igt_vec <- seurat_meta_filtered_spleen[spleen_cells, "IGT"]
selected_igts <- names(table(igt_vec))[table(igt_vec) >= 500]
igt_mean_mat <- do.call(rbind, lapply(selected_igts, function(igt) {
  colMeans(L_pm_filtered[spleen_cells[igt_vec == igt], , drop = FALSE])
}))
rownames(igt_mean_mat) <- selected_igts

gp_igt_var <- apply(igt_mean_mat, 2, var)
gp_overall <- colMeans(igt_mean_mat)

# ============================================================
# S1C: GP mean-of-IGT-mean-loading vs. variance-of-IGT-mean-loading (spleen)
# ============================================================
gp_stats <- data.frame(GP = colnames(igt_mean_mat), x = gp_overall, y = gp_igt_var) %>%
  arrange(desc(y)) %>%
  mutate(label = ifelse(row_number() <= 10, gp_label(as.character(GP)), ""))

p_S1C <- ggplot(gp_stats, aes(x = x, y = y, label = label)) +
  geom_point(size = 1.5, alpha = 0.7, color = "steelblue") +
  ggrepel::geom_text_repel(seed = 42, size = 3, box.padding = 0.4, max.overlaps = Inf, segment.color = "grey50") +
  cowplot::theme_cowplot() +
  labs(
    title = "GP Mean of IGT Mean Loading vs. Between-IGT VAR",
    x = "Mean of IGT Mean Loading",
    y = "Variance of IGT Mean Loading"
  )
ggsave(paste0(figure_path, "S1C.pdf"), plot = p_S1C, width = 6, height = 5, dpi = 300)

# ============================================================
# S1D: heatmap of per-IGT mean loading for S1C's ten highest-variance GPs
# ============================================================
# Same `igt_mean_mat` and same `gp_igt_var` as S1C, so these ten rows are the
# ten GPs S1C labels -- see the note where that matrix is built.
top10_var_gps <- names(sort(gp_igt_var, decreasing = TRUE))[1:10]

plot_mat <- t(igt_mean_mat[, top10_var_gps, drop = FALSE])
plot_mat[plot_mat < 0] <- 0
rownames(plot_mat) <- gp_label(rownames(plot_mat))

pdf(paste0(figure_path, "S1D.pdf"), width = 5, height = 5)
pheatmap(
  plot_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Top 10 GPs by Variance of IGT Mean Loading",
  color = colorRampPalette(c("white", "red"))(100),
  border_color = "white",
  fontsize_row = 8,
  angle_col = 45
)
dev.off()

# ============================================================
# S1E: active-gene vs active-cell scatter per GP, using the SAME hard-threshold
# definitions as Figure 2 (per-GP-normalized): number of active genes = count of
# genes with |score| > 0.25 of the GP's max; proportion of active cells = fraction
# of cells with loading > 0.1 of the GP's max. Non-thymocyte cells, matching
# Figure 2. (Replaces the earlier EBMF-sparsity-prior version.)
# ============================================================
non_thymo_s1e <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 != "thymocyte"]
L_s1e <- L_pm_filtered[non_thymo_s1e, ]
L_norm_s1e <- L_s1e / matrix(apply(L_s1e, 2, max), nrow = nrow(L_s1e), ncol = ncol(L_s1e), byrow = TRUE)
prop_cells <- colSums(L_norm_s1e > 0.1) / nrow(L_norm_s1e)   # proportion of active cells per GP

F_s1e <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
F_norm_s1e <- F_s1e / matrix(apply(F_s1e, 2, function(x) max(abs(x))), nrow = nrow(F_s1e), ncol = ncol(F_s1e), byrow = TRUE)
n_genes_act <- colSums(abs(F_norm_s1e) > 0.25)               # number of active genes per GP

scatter_df_s1e <- data.frame(n_genes = n_genes_act, prop_cells = prop_cells)
pct_breaks <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 1)
p_S1E <- ggplot(scatter_df_s1e, aes(x = n_genes, y = prop_cells)) +
  geom_point(size = 2, alpha = 0.7, color = "steelblue") +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(breaks = pct_breaks, labels = function(x) paste0(x * 100, "%")) +
  annotation_logticks(sides = "bl") +
  labs(
    x = "Number of active genes per GP (log scale)",
    y = "Proportion of active cells per GP (log scale)",
    title = "Active genes vs. active-cell proportion per GP"
  ) +
  theme_minimal(base_size = 13)
ggsave(paste0(figure_path, "S1E.pdf"), plot = p_S1E, width = 6, height = 5, dpi = 300)
