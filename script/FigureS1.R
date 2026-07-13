# Figure S1. GP reproducibility across IGTs.
#
# Panels produced (see figures/final-selected/bits/Figure S1/FigureS1_caption.md
# for the full caption text):
#   S1A  Cumulative number of GPs validated (cosine >= threshold, thresholds
#        0.2-0.8) as IGTs are added one at a time, in IGT index order.
#   S1B  Number of GPs validated by at least X IGTs, vs X (log-log), for the
#        same thresholds.
#   S1C  Between-IGT variability: each GP's mean-of-per-IGT-mean-loading (x)
#        vs. variance-of-per-IGT-mean-loading (y), spleen-standard subset.
#   S1D  Heatmap of per-IGT mean loading for the 10 highest-variance GPs from
#        S1C (IGTs with >= 500 spleen-standard cells only).
#   S1E  Scatter of the EBMF-learned expected NUMBER of active genes (x) vs.
#        expected proportion of active cells (y) per GP -- from the (1 - pi0)
#        sparsity parameter of each factor's point-mass prior (times the total
#        gene count for the x-axis), not a hard threshold on the fitted
#        matrices. One dot per GP.
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

data_path <- "data/"
figure_path <- "figures/generated/Figure S1/"
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
all_gps <- paste0("K", 1:200)
selected_igts <- names(table(seurat_meta_filtered_spleen$IGT))[table(seurat_meta_filtered_spleen$IGT) >= 500]

# ============================================================
# S1C: GP mean-of-IGT-mean-loading vs. variance-of-IGT-mean-loading (spleen)
# ============================================================
common_cells_c <- intersect(rownames(L_pm_filtered), rownames(seurat_meta_filtered_spleen))
L_sub_c <- L_pm_filtered[common_cells_c, ]
igt_vec_c <- seurat_meta_filtered_spleen[common_cells_c, "IGT"]
igt_levels_c <- unique(igt_vec_c)
igt_mat_c <- do.call(rbind, lapply(igt_levels_c, function(igt) colMeans(L_sub_c[igt_vec_c == igt, , drop = FALSE])))

gp_igt_var <- apply(igt_mat_c, 2, var)
gp_overall <- colMeans(igt_mat_c)
gp_stats <- data.frame(GP = colnames(L_sub_c), x = gp_overall, y = gp_igt_var) %>%
  arrange(desc(y)) %>%
  mutate(label = ifelse(row_number() <= 10, gp_label(as.character(GP)), ""))

p_S1C <- ggplot(gp_stats, aes(x = x, y = y, label = label)) +
  geom_point(size = 1.5, alpha = 0.7, color = "steelblue") +
  ggrepel::geom_text_repel(size = 3, box.padding = 0.4, max.overlaps = Inf, segment.color = "grey50") +
  cowplot::theme_cowplot() +
  labs(
    title = "GP Mean of IGT Mean Loading vs. Between-IGT VAR",
    x = "Mean of IGT Mean Loading",
    y = "Variance of IGT Mean Loading"
  )
ggsave(paste0(figure_path, "S1C.pdf"), plot = p_S1C, width = 6, height = 5, dpi = 300)

# ============================================================
# S1D: heatmap of per-IGT mean loading for the top-10 variance GPs from S1C
#      (IGTs with >= 500 spleen-standard cells only)
# ============================================================
top10_var_gps <- names(sort(gp_igt_var, decreasing = TRUE))[1:10]

spleen_cells_act <- intersect(rownames(L_pm_filtered), rownames(seurat_meta_filtered_spleen))
igt_vec_act <- seurat_meta_filtered_spleen[spleen_cells_act, "IGT"]
igt_mean_mat_act <- do.call(rbind, lapply(selected_igts, function(igt) {
  cells_i <- spleen_cells_act[igt_vec_act == igt]
  colMeans(L_pm_filtered[cells_i, all_gps, drop = FALSE])
}))
rownames(igt_mean_mat_act) <- selected_igts

plot_mat <- t(igt_mean_mat_act[, top10_var_gps, drop = FALSE])
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
# S1E: EBMF sparsity scatter -- expected proportion of active genes (x)
# vs. expected proportion of active cells (y) per GP, taken directly from
# each factor's point-mass prior weight (p = 1 - pi0).
# (Same quantities as Figure 1's scatter_e_active_cells_vs_genes panel,
# re-drawn here in ggplot with log-log axes.)
# ============================================================
load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
n_genes  <- nrow(readRDS(paste0(data_path, "F_pm_filtered.rds")))  # total genes in the factorization

n_gp     <- length(flashier_snmf_fitted_prior$L_ghat)
l_pi_vec <- sapply(seq_len(n_gp), function(i) flashier_snmf_fitted_prior$L_ghat[[i]]$pi[1])
f_pi_vec <- sapply(seq_len(n_gp), function(i) flashier_snmf_fitted_prior$F_ghat[[i]]$pi[1])
p_cells     <- 1 - l_pi_vec               # expected proportion of active cells per GP
n_genes_act <- (1 - f_pi_vec) * n_genes   # expected NUMBER of active genes per GP

scatter_df_s1e <- data.frame(n_genes = n_genes_act, prop_cells = p_cells)

pct_breaks <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 1)
p_S1E <- ggplot(scatter_df_s1e, aes(x = n_genes, y = prop_cells)) +
  geom_point(size = 2, alpha = 0.7, color = "steelblue") +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(breaks = pct_breaks, labels = function(x) paste0(x * 100, "%")) +
  annotation_logticks(sides = "bl") +
  labs(
    x = "Expected number of active genes per GP (log scale)",
    y = "Expected proportion of active cells per GP (log scale)",
    title = "Expected active genes vs. active-cell proportion per GP\n(EBMF sparsity priors)"
  ) +
  theme_minimal(base_size = 13)
ggsave(paste0(figure_path, "S1E.pdf"), plot = p_S1E, width = 6, height = 5, dpi = 300)
