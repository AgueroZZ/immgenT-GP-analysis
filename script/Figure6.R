# Figure 6. GPs and tissue.
#
# --- internal ---
# NOTE the renumbering: this figure's published counterpart is
# figures/Previous/bits/Figure *4* (4a-4e), not its "Figure 6" -- that one is
# the published CITE-seq figure, which is ours Figure 7. Full caption text:
# ../immgen-t-factors/figures/Figure_Organ/Figure_Organ_caption.md.
# --- end internal ---
# Panels produced:
#   6a  Row-centered mean GP activity across the 18 tissues, for the 32
#       tissue-associated GPs, tissues in organ_simplified factor order.
#   6b  The same 32 GPs and row order, across the 8 lineages.
#   6c  GP37+ rate by lineage, mammary gland vs. the same lineage elsewhere.
#   6d  Organ AUC vs Level-2 (fine-grained sub-lineage/cluster) AUC, with the
#       7 organ-specific GPs (red) and a contrasting cluster-specific set
#       (blue) highlighted.
#   6e  Alluvial diagram: organ of origin -> GP -> Level-2 cell type, for
#       GP+ cells of the 7 organ-specific GPs.
# --- internal ---
# Reworked 2026-09-10 from the experiments/fig_n5/ draft. Two panels were
# dropped: the organ-vs-level-1 AUC scatter that used to be 6a, and the
# marker-gene heatmap that used to be 6c (published 4a and 4c). The two
# heatmaps are new; the surviving three are the old 6b, 6d and 6e, unchanged,
# with 6b relettered to 6c.
# --- end internal ---
#
# --- internal ---
# Source: ported from Figure_Organ.R, which mixed these 5 panels with
# other exploratory analyses (extra AUC scatter variants, per-organ ROC
# curves, a broken/undefined-object "gp_decomposition.pdf" panel) that are
# dropped here since they don't correspond to a final figure panel.
#
# --- end internal ---
# Required inputs (data/) -- see code/README.md's "Data provenance" table
# for the full picture:
#   L_pm_filtered.rds                        [code/pipeline/01b_filter_cells.R]
#   igt1_96_..._ADTonly.Rds                  [primary input Seurat object]
#   level_2_AUC_list_figure_no_thymocytes_healthy.rds,
#   organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds
#     [code/pipeline/02_compute_auc.R]

library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(Matrix)
library(viridis)
library(cowplot)
library(ggalluvial)
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ZemmourLib)
})

data_path <- "data/"
figure_path <- "figures/final-selected/Figure 6/"
run_started_at <- Sys.time()

# Group means, row centering, dominant-group ordering and the heatmap call for
# 6a/6b; shared with Extended Data Figure 2d, which draws the cluster version.
source("code/R/centered_mean_heatmap.R")

# ============================================================
# Load data (healthy, non-thymocyte reference)
# ============================================================
level_2_AUC_list <- readRDS(paste0(
  data_path, "level_2_AUC_list_figure_no_thymocytes_healthy.rds"
))
organ_AUC_list <- readRDS(paste0(
  data_path, "organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds"
))
# Metadata is read from the Seurat object directly.
# --- internal ---
# Not from the cached data/seurat_meta.rds, which is stale -- see
# code/R/setup_data.R for why.
# --- end internal ---
seurat_meta <- readRDS(paste0(
  data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]

# Rename K## to GP## for display consistency
colnames(L_pm_filtered) <- gsub("^K", "GP", colnames(L_pm_filtered))
colnames(level_2_AUC_list$auc) <- gsub("^K", "GP", colnames(level_2_AUC_list$auc))
colnames(level_2_AUC_list$threshold) <- gsub("^K", "GP", colnames(level_2_AUC_list$threshold))
colnames(organ_AUC_list$auc) <- gsub("^K", "GP", colnames(organ_AUC_list$auc))
colnames(organ_AUC_list$threshold) <- gsub("^K", "GP", colnames(organ_AUC_list$threshold))

# Restrict reference to healthy, non-thymocyte cells
seurat_meta_filtered_no_thymocytes_healthy <- seurat_meta_filtered %>%
  filter(annotation_level1 != "thymocyte", condition_broad == "healthy")

# The 7 organ-specific GPs highlighted throughout this figure (caption 6d/6e)
gps_of_interest <- c("GP3", "GP6", "GP11", "GP26", "GP29", "GP37", "GP177")

# Column order for 6a/6b: levels(so_orig$organ_simplified) and the Figure 1
# lineage order. Organ levels with no healthy non-thymocyte cells drop out.
organ_display_levels <- c(
  "blood", "spleen", "LN", "SLO", "bone marrow", "thymus",
  "lung", "liver", "peritoneal cavity",
  "colon epi", "small intestine epi", "colon LP", "small intestine LP",
  "skin", "mammary gland", "uterus", "placenta", "prostate",
  "submandibular gland", "kidney", "pancreas", "CNS", "synovial fluid"
)
level1_display_levels <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN", "DP")
organ_color_limit <- 0.2
level1_color_limit <- 0.1

# Labels a highlighted point with its top categories above `threshold` AUC.
top_cats_label <- function(factor_name, auc_matrix, positive_mask, threshold = 0.85, n = 3) {
  vals <- auc_matrix[, factor_name]
  vals <- vals[positive_mask[, factor_name]]
  vals <- sort(vals[vals > threshold], decreasing = TRUE)
  cats <- names(vals)[seq_len(min(n, length(vals)))]
  if (length(cats) == 0) return(factor_name)
  paste0(factor_name, ":\n", paste(cats, collapse = "\n"))
}

# ============================================================
# Shared AUC setup: the organ AUC matrix and the healthy reference means
# ============================================================
# Organs with fewer than 100 healthy non-thymocyte cells are dropped, matching
# the filter every AUC panel in this figure uses.
organ_AUC <- organ_AUC_list$auc
organ_small_count <- table(seurat_meta_filtered_no_thymocytes_healthy$organ_simplified)
organ_small_count <- names(organ_small_count[organ_small_count < 100])
organ_AUC <- organ_AUC[!rownames(organ_AUC) %in% organ_small_count, ]

healthy_cells <- rownames(seurat_meta_filtered_no_thymocytes_healthy)
L_healthy <- L_pm_filtered[healthy_cells, ]
overall_mean <- colMeans(L_healthy, na.rm = TRUE)
stopifnot(nrow(L_healthy) == nrow(seurat_meta_filtered_no_thymocytes_healthy), !anyNA(L_healthy))

# ============================================================
# 6a/6b selection: the 32 tissue-associated GPs
# ============================================================
# The union of the two criteria the manuscript uses for "tissue-associated":
#
#   A  raw (uncentered) mean loading >= 0.1 in at least one of the 18 tissues
#      -- 31 GPs, the "programs active across tissues";
#   B  one-vs-rest tissue AUC > 0.9 under the positivity mask (mean loading in
#      the tissue above the GP's overall mean) -- 7 GPs, the organ-specific set
#      this figure highlights in 6d and 6e.
#
# Six GPs satisfy both, so 32 are shown. GP37 is the only member of B alone:
# its mammary-gland AUC is 0.9918, but it is active in 2.5% of mammary gland
# cells, which dilutes its mean loading there to 0.011 -- an order of magnitude
# under A's cutoff. Its row is therefore near-white in 6a and 6b by
# construction, which is why it is parked on the last row rather than left as a
# hole inside the mammary gland block; 6c is the panel that shows what it does.
#
# Rule B is re-derived rather than reusing gps_of_interest, and then asserted
# equal to it, so the selection and the highlighting cannot drift apart.
raw_mean_cutoff <- 0.1
auc_cutoff <- 0.9
last_row_gp <- "GP37"

organ_raw <- mean_loading_by_group(L_healthy, seurat_meta_filtered_no_thymocytes_healthy$organ_simplified)$matrix
level1_raw <- mean_loading_by_group(L_healthy, seurat_meta_filtered_no_thymocytes_healthy$annotation_level1)$matrix

tissue_active <- rowSums(organ_raw >= raw_mean_cutoff) > 0L

organ_AUC_positive <- sweep(
  t(sapply(rownames(organ_AUC), function(cat) {
    colMeans(L_healthy[seurat_meta_filtered_no_thymocytes_healthy$organ_simplified == cat, , drop = FALSE], na.rm = TRUE)
  })),
  2, overall_mean, "-"
) > 0
tissue_specific <- rownames(organ_raw) %in%
  colnames(organ_AUC)[colSums((organ_AUC > auc_cutoff) & organ_AUC_positive) > 0L]

selected <- tissue_active | tissue_specific
tissue_associated_gps <- rownames(organ_raw)[selected]
expected_gps <- paste0("GP", c(
  1, 3, 4, 6, 8, 9, 11, 22, 23, 25, 26, 29, 30, 32, 35, 37, 41, 43, 49, 51, 58,
  62, 63, 72, 80, 93, 100, 166, 170, 171, 174, 177
))
stopifnot(
  setequal(rownames(organ_raw)[tissue_specific], gps_of_interest),
  sum(tissue_active) == 31L,
  sum(tissue_specific) == 7L,
  identical(tissue_associated_gps, expected_gps),
  last_row_gp %in% expected_gps
)

organ_selected_raw <- organ_raw[selected, , drop = FALSE]
level1_selected_raw <- level1_raw[selected, , drop = FALSE]
organ_centered <- center_by_gp_mean(organ_selected_raw)
level1_centered <- center_by_gp_mean(level1_selected_raw)

# Columns follow levels(so_orig$organ_simplified) and levels(annotation_level1)
# rather than a dominant-group ordering. Five organ levels (SLO, thymus,
# prostate, pancreas, synovial fluid) have no healthy non-thymocyte cells and
# drop out; the other 18 keep their relative order. Rows are dominant-group
# blocks computed against that fixed column order, with GP37 appended last.
observed_organs <- intersect(organ_display_levels, colnames(organ_centered))
if (!setequal(observed_organs, colnames(organ_centered))) {
  stop("A tissue in the reference is missing from organ_display_levels: ",
       paste(setdiff(colnames(organ_centered), organ_display_levels), collapse = ", "))
}
organ_column_order <- match(observed_organs, colnames(organ_centered))

block_rows <- setdiff(rownames(organ_selected_raw), last_row_gp)
block_order <- dominant_group_order(
  organ_selected_raw[block_rows, observed_organs, drop = FALSE],
  fixed_column_order = seq_along(observed_organs)
)
organ_row_order <- match(
  c(block_rows[block_order$row_order], last_row_gp), rownames(organ_centered)
)
level1_column_order <- match(level1_display_levels, colnames(level1_centered))
level1_row_order <- match(rownames(organ_centered)[organ_row_order], rownames(level1_centered))

stopifnot(
  nrow(organ_centered) == 32L, ncol(organ_centered) == 18L, ncol(level1_centered) == 8L,
  identical(rownames(organ_centered), rownames(level1_centered)),
  !anyNA(organ_column_order), !anyNA(level1_column_order), !anyNA(organ_row_order),
  identical(rownames(organ_centered)[organ_row_order[32]], last_row_gp),
  max(abs(rowMeans(organ_centered))) < 1e-12,
  max(abs(rowMeans(level1_centered))) < 1e-12
)

# ============================================================
# 6a: the 32 tissue-associated GPs across the 18 tissues
# ============================================================
render_centered_heatmap(
  organ_centered,
  palette_for_groups(colnames(organ_centered), ZemmourLib::immgent_colors$organ_simplified, "organ_simplified"),
  "tissue (organ_simplified)",
  paste0(figure_path, "6a.pdf"),
  organ_row_order, organ_column_order, organ_color_limit,
  "32 tissue-associated GPs; organ_simplified order",
  palette = heatmap_palettes$blue_red
)

# ============================================================
# 6b: the same 32 GPs across the 8 lineages
# ============================================================
# Green-white-purple at +/-0.1 rather than blue-white-red at +/-0.2: spreading a
# program over 8 lineages instead of 18 tissues gives correspondingly smaller
# deviations, and 6a's scale renders this panel almost blank. Purple is the
# positive end.
render_centered_heatmap(
  level1_centered,
  palette_for_groups(colnames(level1_centered), ZemmourLib::immgent_colors$level1, "annotation_level1"),
  "lineage (annotation_level1)",
  paste0(figure_path, "6b.pdf"),
  level1_row_order, level1_column_order, level1_color_limit,
  "the same 32 GPs, row order from (a)",
  palette = heatmap_palettes$green_purple
)

# ============================================================
# 6d prep: Max AUC Organ vs Level-2
# ============================================================
level_2_AUC <- level_2_AUC_list$auc
level_2_small_count <- table(seurat_meta_filtered_no_thymocytes_healthy$annotation_level2)
level_2_small_count <- names(level_2_small_count[level_2_small_count < 100])
level_2_AUC <- level_2_AUC[!rownames(level_2_AUC) %in% level_2_small_count, ]

level_2_cat_mean <- t(sapply(rownames(level_2_AUC), function(cat) {
  idx <- seurat_meta_filtered_no_thymocytes_healthy$annotation_level2 == cat
  colMeans(L_healthy[idx, , drop = FALSE], na.rm = TRUE)
}))
level_2_AUC_positive <- sweep(level_2_cat_mean, 2, overall_mean, "-") > 0

# 6d/6e reuse `organ_AUC_max_name`, but recomputed against the Level-2
# category-count filter.
organ_AUC_masked_l2 <- organ_AUC
organ_AUC_positive_l2 <- sweep(
  t(sapply(rownames(organ_AUC), function(cat) {
    idx <- seurat_meta_filtered_no_thymocytes_healthy$organ_simplified == cat
    colMeans(L_healthy[idx, , drop = FALSE], na.rm = TRUE)
  })),
  2, overall_mean, "-"
) > 0
organ_AUC_masked_l2[!organ_AUC_positive_l2] <- NA
organ_AUC_max <- apply(organ_AUC_masked_l2, 2, max, na.rm = TRUE)
organ_AUC_max_name <- apply(organ_AUC_masked_l2, 2, function(x) rownames(organ_AUC_masked_l2)[which.max(x)])

# 6d's own max-AUC table, over the Level-2 categories.
# --- internal ---
# It must not reuse the organ-vs-level-1 `df` that the retired 6a built, whose
# Max_AUC_Level1 column holds the Level-1 maxima -- doing so silently plots
# Level-1 AUC on this panel's "Max AUC (Level-2)" axis. That panel is gone as
# of 2026-09-10, but the trap survives in Figure_Organ.R.
#
# The original Figure_Organ.R rebuilds `max_AUC_df`/`df` at this point from
# `table_level_2_AUC`, storing the Level-2 maxima in a column it still calls
# `Max_AUC_Level1` -- a misleading name we drop here in favour of
# `Max_AUC_Level2`.
# --- end internal ---
level_2_AUC_masked <- level_2_AUC
level_2_AUC_masked[!level_2_AUC_positive] <- NA
level_2_AUC_max <- apply(level_2_AUC_masked, 2, max, na.rm = TRUE)
level_2_AUC_max_name <- apply(level_2_AUC_masked, 2, function(x) {
  rownames(level_2_AUC_masked)[which.max(x)]
})
o_l2 <- order(level_2_AUC_max, decreasing = TRUE)
table_level_2_AUC <- data.frame(
  Factor = colnames(level_2_AUC)[o_l2],
  Max_AUC = level_2_AUC_max[o_l2],
  Annotation = level_2_AUC_max_name[o_l2]
)
df_l2 <- data.frame(
  Factor = table_level_2_AUC$Factor,
  annotation_Level2 = table_level_2_AUC$Annotation,
  annotation_Organ = organ_AUC_max_name[match(table_level_2_AUC$Factor, names(organ_AUC_max))],
  Max_AUC_Organ = organ_AUC_max[match(table_level_2_AUC$Factor, names(organ_AUC_max))],
  Max_AUC_Level2 = table_level_2_AUC$Max_AUC
) %>%
  mutate(residual = Max_AUC_Level2 - Max_AUC_Organ, abs_res = abs(residual))

# ============================================================
# 6d: Max AUC Organ vs Level-2, 7 organ-specific GPs (red) vs.
#     contrasting cluster-specific GPs (blue) highlighted
# ============================================================
seven_gp_df <- df_l2 |>
  dplyr::filter(Factor %in% gps_of_interest) |>
  dplyr::mutate(label_text = sapply(Factor, top_cats_label, auc_matrix = level_2_AUC, positive_mask = level_2_AUC_positive, threshold = 0.9, n = 3))

top_left_gps <- c("GP14", "GP36", "GP16", "GP151", "GP21", "GP122", "GP2", "GP171", "GP5", "GP13")
top_left_df <- df_l2 |>
  dplyr::filter(Factor %in% top_left_gps) |>
  dplyr::mutate(label_text = sapply(Factor, top_cats_label, auc_matrix = level_2_AUC, positive_mask = level_2_AUC_positive, threshold = 0.9, n = 3))

p_6d <- ggplot(df_l2, aes(Max_AUC_Organ, Max_AUC_Level2)) +
  geom_point(alpha = 0.2, size = 1.5, color = "grey60") +
  geom_point(data = top_left_df, color = "#1f78b4", size = 2.2, alpha = 0.9) +
  geom_text_repel(
    seed = 42,
    data = top_left_df, aes(label = label_text), color = "#1f78b4", lineheight = 0.85, size = 2.5,
    direction = "y", nudge_x = -0.1, segment.color = "#1f78b4",
    arrow = arrow(length = unit(0.008, "npc"), type = "closed", angle = 20),
    force = 4, force_pull = 0.05, box.padding = 0.5, point.padding = 0.15,
    max.time = 10, max.iter = 2e4, max.overlaps = 30, min.segment.length = 0.01, segment.alpha = 0.7
  ) +
  geom_point(data = seven_gp_df, color = "#e31a1c", size = 2.2, alpha = 0.9) +
  geom_text_repel(
    seed = 42,
    data = seven_gp_df, aes(label = label_text), color = "#e31a1c", size = 2.5, lineheight = 0.85,
    direction = "y", nudge_x = 0.18, xlim = c(1.0, NA), segment.color = "#e31a1c",
    arrow = arrow(length = unit(0.008, "npc"), type = "closed", angle = 20),
    force = 6, force_pull = 0.02, box.padding = 0.6, point.padding = 0.15,
    max.time = 10, max.iter = 2e4, max.overlaps = 30, min.segment.length = 0.01, segment.alpha = 0.7
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(0.46, 1.04), ylim = c(0.5, 1.02), expand = FALSE, clip = "off") +
  labs(x = "Max AUC (Organ Simplified)", y = "Max AUC (Level-2)", title = "Max AUC: Organ vs Level-2 - organ-specific GPs") +
  theme_minimal(base_size = 13) +
  theme(plot.margin = margin(10, 80, 10, 80))
ggsave(filename = paste0(figure_path, "6d.pdf"), plot = p_6d, width = 8, height = 8, dpi = 300)

# ============================================================
# 6c: GP37+ rate by lineage, mammary gland vs. elsewhere
# ============================================================
plot_gp_threshold_group_activation_rate <- function(
  gp, organ, threshold, loading_mat, organ_info, group_info,
  group_label = "Level-2", base_size = 13, min_in_organ = 10,
  group_colors = ZemmourLib::immgent_colors$level2, fallback_group_color = "grey60",
  reference = c("not_in_group", "not_in_organ")
) {
  reference <- match.arg(reference)
  if (!gp %in% colnames(loading_mat)) stop(sprintf("GP '%s' not found in loading matrix.", gp))
  if (!organ %in% organ_info) stop(sprintf("Organ '%s' not found in organ_info.", organ))

  loading <- loading_mat[, gp]
  keep <- !(is.na(loading) | is.na(organ_info) | is.na(group_info))
  loading <- loading[keep]
  organ_info <- organ_info[keep]
  group_info <- as.character(group_info[keep])

  in_organ <- organ_info == organ
  positive <- loading > threshold
  group_levels <- sort(unique(group_info))

  rate_df <- data.frame(
    group = group_levels,
    n_in_organ = vapply(group_levels, function(l) sum(group_info == l & in_organ), integer(1)),
    n_pos_in_organ = vapply(group_levels, function(l) sum(group_info == l & in_organ & positive), integer(1))
  )

  if (reference == "not_in_group") {
    rate_df$n_ref <- vapply(group_levels, function(l) sum(group_info != l & in_organ), integer(1))
    rate_df$n_pos_ref <- vapply(group_levels, function(l) sum(group_info != l & in_organ & positive), integer(1))
    ref_label <- "Not in group (same organ)"
    title_vs <- sprintf("%s vs. same-organ non-group", organ)
  } else {
    rate_df$n_ref <- vapply(group_levels, function(l) sum(group_info == l & !in_organ), integer(1))
    rate_df$n_pos_ref <- vapply(group_levels, function(l) sum(group_info == l & !in_organ & positive), integer(1))
    ref_label <- "Not in organ (same group)"
    title_vs <- sprintf("%s vs. same-group non-organ", organ)
  }

  rate_df$rate_in_organ <- rate_df$n_pos_in_organ / rate_df$n_in_organ
  rate_df$rate_ref <- rate_df$n_pos_ref / rate_df$n_ref
  rate_df <- rate_df[rate_df$n_in_organ >= min_in_organ, , drop = FALSE]
  if (nrow(rate_df) == 0) stop(sprintf("No %s type has >= %d cells in '%s'.", group_label, min_in_organ, organ))

  long_df <- data.frame(
    group = rep(rate_df$group, 2),
    type = factor(rep(c("In organ", ref_label), each = nrow(rate_df)), levels = c("In organ", ref_label)),
    rate = c(rate_df$rate_in_organ, rate_df$rate_ref)
  )
  level_order <- rate_df$group[order(rate_df$rate_in_organ, decreasing = TRUE)]
  long_df$group <- factor(long_df$group, levels = level_order)

  fill_values <- group_colors[as.character(level_order)]
  missing_colors <- is.na(fill_values)
  if (any(missing_colors)) {
    fill_values[missing_colors] <- fallback_group_color
    warning(sprintf(
      "%s annotations missing from group_colors and colored %s: %s",
      group_label, fallback_group_color, paste(level_order[missing_colors], collapse = ", ")
    ))
  }
  alpha_vals <- c(1, 0.35)
  names(alpha_vals) <- c("In organ", ref_label)

  ggplot(long_df, aes(x = group, y = rate, fill = group, alpha = type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75, color = "grey35", linewidth = 0.15) +
    scale_fill_manual(values = fill_values, guide = "none") +
    scale_alpha_manual(values = alpha_vals, guide = guide_legend(override.aes = list(fill = "grey40"))) +
    labs(
      x = sprintf("%s annotation", group_label),
      y = sprintf("Proportion of cells with %s > %.3g", gp, threshold),
      alpha = NULL,
      title = sprintf("%s+ rate by %s: %s", gp, group_label, title_vs),
      subtitle = sprintf("threshold = %.3g; %s types with < %d cells in %s dropped", threshold, group_label, min_in_organ, organ)
    ) +
    theme_minimal(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
}

p_6c <- plot_gp_threshold_group_activation_rate(
  gp = "GP37",
  organ = "mammary gland",
  threshold = organ_AUC_list$threshold["mammary gland", "GP37"],
  min_in_organ = 100,
  loading_mat = L_pm_filtered[rownames(seurat_meta_filtered_no_thymocytes_healthy), ],
  organ_info = seurat_meta_filtered_no_thymocytes_healthy$organ_simplified,
  group_info = seurat_meta_filtered_no_thymocytes_healthy$annotation_level1,
  group_label = "Level-1",
  group_colors = ZemmourLib::immgent_colors$level1,
  reference = "not_in_organ"
)
ggsave(filename = paste0(figure_path, "6c.pdf"), plot = p_6c, width = 8, height = 5, dpi = 300)

# ============================================================
# 6e: alluvial, organ -> GP -> Level-2, for GP+ cells of the
#     7 organ-specific GPs
# ============================================================
best_organ_per_gp <- organ_AUC_max_name[gps_of_interest]
gp_thresholds <- mapply(function(gp, organ) organ_AUC_list$threshold[organ, gp], gps_of_interest, best_organ_per_gp)
names(gp_thresholds) <- gps_of_interest

n_cap_gp <- 300
set.seed(42)
alluvial_rows <- lapply(gps_of_interest, function(gp) {
  positive_idx <- L_healthy[, gp] > gp_thresholds[gp]
  meta_pos <- seurat_meta_filtered_no_thymocytes_healthy[positive_idx, ]
  d <- data.frame(gp_program = gp, organ = meta_pos$organ_simplified, level2 = meta_pos$annotation_level2, stringsAsFactors = FALSE)
  if (nrow(d) > n_cap_gp) d <- dplyr::slice_sample(d, n = n_cap_gp)
  d
})

count_df <- do.call(rbind, alluvial_rows) |>
  dplyr::count(organ, gp_program, level2, name = "n") |>
  dplyr::filter(!is.na(organ), !is.na(level2), n >= 5)

organ_order <- count_df |> dplyr::summarise(total = sum(n), .by = organ) |> dplyr::arrange(dplyr::desc(total)) |> dplyr::pull(organ)
level2_order <- count_df |> dplyr::summarise(total = sum(n), .by = level2) |> dplyr::arrange(dplyr::desc(total)) |> dplyr::pull(level2)

count_df <- count_df |>
  dplyr::mutate(
    organ = factor(organ, levels = rev(organ_order)),
    gp_program = factor(gp_program, levels = rev(gps_of_interest)),
    level2 = factor(level2, levels = rev(level2_order))
  )

gp_colors <- ZemmourLib::immgent_colors$organ_simplified[unname(best_organ_per_gp)]
gp_colors[is.na(gp_colors)] <- "grey60"
names(gp_colors) <- gps_of_interest

p_6e <- ggplot(count_df, aes(axis1 = organ, axis2 = gp_program, axis3 = level2, y = n)) +
  ggalluvial::geom_alluvium(aes(fill = gp_program), width = 1 / 4, alpha = 0.6, knot.pos = 0.4) +
  ggalluvial::geom_stratum(width = 1 / 4, fill = "grey92", color = "grey50", linewidth = 0.3) +
  ggplot2::geom_text(stat = ggalluvial::StatStratum, aes(label = after_stat(stratum)), size = 3, angle = 90) +
  scale_fill_manual(values = gp_colors, guide = "none") +
  scale_x_discrete(limits = c("Organ", "GP", "Level-2"), expand = c(0.12, 0.12)) +
  labs(y = "Number of GP+ cells", title = "GP+ cells: organ origin and cell type") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  coord_flip()
ggsave(filename = paste0(figure_path, "6e.pdf"), plot = p_6e, width = 20, height = 10, dpi = 300)

# ============================================================
# Every panel this run should have written
# ============================================================
# A panel that fails to write leaves the previous PDF in place and the script
# still exits 0, so check the files rather than the exit code.
expected_panels <- paste0(figure_path, c("6a", "6b", "6c", "6d", "6e"), ".pdf")
stale <- expected_panels[!file.exists(expected_panels) |
                           file.mtime(expected_panels) < run_started_at |
                           file.size(expected_panels) == 0]
if (length(stale) > 0L) {
  stop("These panels were not written by this run: ", paste(stale, collapse = ", "))
}
message(sprintf("wrote %d panels to %s", length(expected_panels), figure_path))
