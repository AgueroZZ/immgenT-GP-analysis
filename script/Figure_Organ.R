library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(fastTopics)
library(dplyr)
library(pheatmap)
library(scales)
library(scattermore)
library(ggpubr)

# ============================================================
# Setup: paths and helper functions
# ============================================================
data_path <- "data/"
code_path <- "code/"
figure_path <- "figures/Figure_Organ/"
source(paste0(code_path, "ROC.R"))
source(paste0(code_path, "filtering_membership.R"))

# ============================================================
# Load data (healthy, non-thymocyte reference)
# ============================================================
level_1_AUC_list <- readRDS(
  file = paste0(data_path, "level_1_AUC_list_figure_no_thymocytes_healthy.rds")
)
level_2_AUC_list <- readRDS(
  file = paste0(data_path, "level_2_AUC_list_figure_no_thymocytes_healthy.rds")
)
organ_AUC_list <- readRDS(paste0(
  data_path,
  "organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds"
))
seurat_meta <- readRDS(paste0(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]

# Rename K## to GP## for display consistency
colnames(L_pm_filtered) <- gsub("^K", "GP", colnames(L_pm_filtered))
colnames(level_1_AUC_list$auc) <- gsub(
  "^K",
  "GP",
  colnames(level_1_AUC_list$auc)
)
colnames(level_2_AUC_list$auc) <- gsub(
  "^K",
  "GP",
  colnames(level_2_AUC_list$auc)
)
colnames(level_2_AUC_list$threshold) <- gsub(
  "^K",
  "GP",
  colnames(level_2_AUC_list$threshold)
)
colnames(organ_AUC_list$auc) <- gsub("^K", "GP", colnames(organ_AUC_list$auc))
colnames(organ_AUC_list$threshold) <- gsub(
  "^K",
  "GP",
  colnames(organ_AUC_list$threshold)
)

saveRDS(
  seurat_meta_filtered,
  file = paste0(data_path, "seurat_meta_filtered.rds")
)
qs::qsave(L_pm_filtered, file = paste0(data_path, "L_pm_filtered.qs"))

# Restrict reference to healthy, non-thymocyte cells
seurat_meta_filtered_no_thymocytes_healthy <- seurat_meta_filtered %>%
  filter(annotation_level1 != "thymocyte", condition_broad == "healthy")


# ============================================================
# Scatterplot: Max AUC Organ vs Level-1
#   -> AUC_Organ_vs_Level1.pdf
# ============================================================
top_n <- 100
# remove categories with small count
level_1_small_count <- table(
  seurat_meta_filtered_no_thymocytes_healthy$annotation_level1
)
level_1_small_count <- names(level_1_small_count[level_1_small_count < 1000])
level_1_AUC <- level_1_AUC_list$auc
level_1_AUC <- level_1_AUC[!rownames(level_1_AUC) %in% level_1_small_count, ]

organ_AUC <- organ_AUC_list$auc
# remove categories with small count
organ_small_count <- table(
  seurat_meta_filtered_no_thymocytes_healthy$organ_simplified
)
organ_small_count <- names(organ_small_count[organ_small_count < 100])
organ_AUC <- organ_AUC[!rownames(organ_AUC) %in% organ_small_count, ]

# Positivity masks: mean loading in category > overall mean → high loading predicts membership
# Computed before max AUC so the axes also reflect positive-direction predictions only
healthy_cells <- rownames(seurat_meta_filtered_no_thymocytes_healthy)
L_healthy <- L_pm_filtered[healthy_cells, ]
overall_mean <- colMeans(L_healthy, na.rm = TRUE)

level_1_cat_mean <- t(sapply(rownames(level_1_AUC), function(cat) {
  idx <- seurat_meta_filtered_no_thymocytes_healthy$annotation_level1 == cat
  colMeans(L_healthy[idx, , drop = FALSE], na.rm = TRUE)
}))
level_1_AUC_positive <- sweep(level_1_cat_mean, 2, overall_mean, "-") > 0

organ_cat_mean <- t(sapply(rownames(organ_AUC), function(cat) {
  idx <- seurat_meta_filtered_no_thymocytes_healthy$organ_simplified == cat
  colMeans(L_healthy[idx, , drop = FALSE], na.rm = TRUE)
}))
organ_AUC_positive <- sweep(organ_cat_mean, 2, overall_mean, "-") > 0

# Max AUC restricted to positively-predicted categories
level_1_AUC_masked <- level_1_AUC
level_1_AUC_masked[!level_1_AUC_positive] <- NA
level_1_AUC_max <- apply(level_1_AUC_masked, 2, max, na.rm = TRUE)
level_1_AUC_max_name <- apply(level_1_AUC_masked, 2, function(x) {
  rownames(level_1_AUC_masked)[which.max(x)]
})
o <- order(level_1_AUC_max, decreasing = TRUE)
table_level_1_AUC <- data.frame(
  Factor = colnames(level_1_AUC)[o],
  Max_AUC = level_1_AUC_max[o],
  Annotation = level_1_AUC_max_name[o]
)

organ_AUC_masked <- organ_AUC
organ_AUC_masked[!organ_AUC_positive] <- NA
organ_AUC_max <- apply(organ_AUC_masked, 2, max, na.rm = TRUE)
organ_AUC_max_name <- apply(organ_AUC_masked, 2, function(x) {
  rownames(organ_AUC_masked)[which.max(x)]
})
o <- order(organ_AUC_max, decreasing = TRUE)
table_organ_AUC <- data.frame(
  Max_AUC = organ_AUC_max[o],
  Annotation = organ_AUC_max_name[o]
)

# merge the two AUC tables by factor
max_AUC_df <- data.frame(
  Factor = table_level_1_AUC$Factor,
  annotation_Level1 = table_level_1_AUC$Annotation,
  annotation_Organ = organ_AUC_max_name[match(
    table_level_1_AUC$Factor,
    names(organ_AUC_max)
  )],
  Max_AUC_Organ = organ_AUC_max[match(
    table_level_1_AUC$Factor,
    names(organ_AUC_max)
  )],
  Max_AUC_Level1 = table_level_1_AUC$Max_AUC
)
df <- max_AUC_df %>%
  mutate(
    residual = Max_AUC_Level1 - Max_AUC_Organ,
    abs_res = abs(residual)
  )

top_cats_label <- function(
  factor_name,
  auc_matrix,
  positive_mask,
  threshold = 0.85,
  n = 3
) {
  vals <- auc_matrix[, factor_name]
  vals <- vals[positive_mask[, factor_name]]
  vals <- sort(vals[vals > threshold], decreasing = TRUE)
  cats <- names(vals)[seq_len(min(n, length(vals)))]
  if (length(cats) == 0) {
    return(factor_name)
  }
  paste0(factor_name, ":\n", paste(cats, collapse = "\n"))
}

# Factors to highlight: AUC > 0.9 in at least one axis (organ or level-1).
# This set is also reused in the Level-2 scatterplot below.
highlighted_factors <- df %>%
  filter(is.finite(residual), Max_AUC_Organ > 0.9 | Max_AUC_Level1 > 0.9) %>%
  pull(Factor)

label_above <- df %>%
  filter(Factor %in% highlighted_factors, residual > 0) %>%
  mutate(
    nudge_x = -0.035,
    label_text = sapply(
      Factor,
      top_cats_label,
      auc_matrix = level_1_AUC,
      positive_mask = level_1_AUC_positive
    )
  )
label_below <- df %>%
  filter(Factor %in% highlighted_factors, residual <= 0) %>%
  mutate(
    nudge_x = 0.035,
    label_text = sapply(
      Factor,
      top_cats_label,
      auc_matrix = organ_AUC,
      positive_mask = organ_AUC_positive
    )
  )
base_plot <- ggplot(df, aes(Max_AUC_Organ, Max_AUC_Level1)) +
  geom_point(alpha = 0.3, size = 1.8) +
  geom_point(data = label_above, color = "#1f78b4", alpha = 0.8, size = 1.8) +
  geom_point(data = label_below, color = "#e31a1c", alpha = 0.8, size = 1.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(0.46, 1.04), ylim = c(0.5, 1.02), expand = FALSE) +
  labs(
    x = "Max AUC (Organ Simplified)",
    y = "Max AUC (Level-1)",
    title = "Max AUC: Organ vs Level-1"
  ) +
  theme_minimal(base_size = 13)

p_all <- base_plot +
  geom_text_repel(
    data = label_above,
    aes(label = label_text),
    color = "#1f78b4",
    size = 2.5,
    lineheight = 0.85,
    direction = "y",
    nudge_x = label_above$nudge_x,
    segment.color = "#1f78b4",
    arrow = arrow(length = unit(0.008, "npc"), type = "closed", angle = 20),
    force = 3,
    force_pull = 0.1,
    box.padding = 0.4,
    point.padding = 0.15,
    max.time = 10,
    max.iter = 2e4,
    max.overlaps = 20,
    min.segment.length = 0.01,
    segment.alpha = 0.7
  ) +
  geom_text_repel(
    data = label_below,
    aes(label = label_text),
    color = "#e31a1c",
    size = 2.5,
    lineheight = 0.85,
    direction = "y",
    nudge_x = label_below$nudge_x,
    segment.color = "#e31a1c",
    arrow = arrow(length = unit(0.008, "npc"), type = "closed", angle = 20),
    force = 3,
    force_pull = 0.1,
    box.padding = 0.4,
    point.padding = 0.15,
    max.time = 10,
    max.iter = 2e4,
    max.overlaps = 20,
    min.segment.length = 0.01,
    segment.alpha = 0.7
  )
p_all
ggsave(
  filename = paste0(figure_path, "AUC_Organ_vs_Level1.pdf"),
  plot = p_all,
  width = 8,
  height = 8,
  dpi = 300
)


# ============================================================
# Scatterplot: Max AUC Organ vs Level-2
#   -> AUC_Organ_vs_Level2.pdf
# ============================================================
top_n <- 10
level_2_AUC <- level_2_AUC_list$auc
# remove categories with small count
level_2_small_count <- table(
  seurat_meta_filtered_no_thymocytes_healthy$annotation_level2
)
level_2_small_count <- names(level_2_small_count[level_2_small_count < 100])
level_2_AUC <- level_2_AUC[!rownames(level_2_AUC) %in% level_2_small_count, ]

organ_AUC <- organ_AUC_list$auc
# remove categories with small count
organ_small_count <- table(
  seurat_meta_filtered_no_thymocytes_healthy$organ_simplified
)
organ_small_count <- names(organ_small_count[organ_small_count < 100])
organ_AUC <- organ_AUC[!rownames(organ_AUC) %in% organ_small_count, ]

# Positivity masks for Level-2 plot (L_healthy and overall_mean reused from Level-1 section)
level_2_cat_mean <- t(sapply(rownames(level_2_AUC), function(cat) {
  idx <- seurat_meta_filtered_no_thymocytes_healthy$annotation_level2 == cat
  colMeans(L_healthy[idx, , drop = FALSE], na.rm = TRUE)
}))
level_2_AUC_positive <- sweep(level_2_cat_mean, 2, overall_mean, "-") > 0

organ_cat_mean_l2 <- t(sapply(rownames(organ_AUC), function(cat) {
  idx <- seurat_meta_filtered_no_thymocytes_healthy$organ_simplified == cat
  colMeans(L_healthy[idx, , drop = FALSE], na.rm = TRUE)
}))
organ_AUC_positive_l2 <- sweep(organ_cat_mean_l2, 2, overall_mean, "-") > 0

# Max AUC restricted to positively-predicted categories
level_2_AUC_masked <- level_2_AUC
level_2_AUC_masked[!level_2_AUC_positive] <- NA
level_2_AUC_max <- apply(level_2_AUC_masked, 2, max, na.rm = TRUE)
level_2_AUC_max_name <- apply(level_2_AUC_masked, 2, function(x) {
  rownames(level_2_AUC_masked)[which.max(x)]
})
o <- order(level_2_AUC_max, decreasing = TRUE)
table_level_2_AUC <- data.frame(
  Factor = colnames(level_2_AUC)[o],
  Max_AUC = level_2_AUC_max[o],
  Annotation = level_2_AUC_max_name[o]
)

organ_AUC_masked_l2 <- organ_AUC
organ_AUC_masked_l2[!organ_AUC_positive_l2] <- NA
organ_AUC_max <- apply(organ_AUC_masked_l2, 2, max, na.rm = TRUE)
organ_AUC_max_name <- apply(organ_AUC_masked_l2, 2, function(x) {
  rownames(organ_AUC_masked_l2)[which.max(x)]
})
o <- order(organ_AUC_max, decreasing = TRUE)
table_organ_AUC <- data.frame(
  Max_AUC = organ_AUC_max[o],
  Annotation = organ_AUC_max_name[o]
)

# merge the two AUC tables by factor
max_AUC_df <- data.frame(
  Factor = table_level_2_AUC$Factor,
  annotation_Level2 = table_level_2_AUC$Annotation,
  annotation_Organ = organ_AUC_max_name[match(
    table_level_2_AUC$Factor,
    names(organ_AUC_max)
  )],
  Max_AUC_Organ = organ_AUC_max[match(
    table_level_2_AUC$Factor,
    names(organ_AUC_max)
  )],
  Max_AUC_Level1 = table_level_2_AUC$Max_AUC
)
df <- max_AUC_df %>%
  mutate(
    residual = Max_AUC_Level1 - Max_AUC_Organ,
    abs_res = abs(residual)
  )

# Manual labels: "GP: main_cat" header + up to 2 additional eligible categories below
manual_cats_label <- function(
  factor_name,
  main_cat,
  auc_matrix,
  positive_mask,
  threshold = 0.85,
  n_extra = 2
) {
  vals <- auc_matrix[, factor_name]
  vals <- vals[positive_mask[, factor_name]]
  vals <- sort(vals[vals > threshold], decreasing = TRUE)
  extra <- names(vals)[!names(vals) %in% main_cat]
  extra <- extra[seq_len(min(n_extra, length(extra)))]
  base <- paste0(factor_name, ": ", main_cat)
  if (length(extra) == 0) {
    return(base)
  }
  paste0(base, "\n", paste(extra, collapse = "\n"))
}

# Highlight the same GP programs as selected in the Level-1 scatterplot above.
label_above <- df %>%
  filter(Factor %in% highlighted_factors, residual > 0) %>%
  mutate(
    nudge_x = -0.035,
    label_text = sapply(
      Factor,
      top_cats_label,
      auc_matrix = level_2_AUC,
      positive_mask = level_2_AUC_positive
    )
  )
label_below <- df %>%
  filter(Factor %in% highlighted_factors, residual <= 0) %>%
  mutate(
    nudge_x = 0.035,
    label_text = sapply(
      Factor,
      top_cats_label,
      auc_matrix = organ_AUC,
      positive_mask = organ_AUC_positive_l2
    )
  )

# Manually emphasized programs for the Level-2 comparison plot.
# Left-top programs have high Level-2 AUC but much lower organ AUC.
# Right-edge programs have the highest organ AUC in the current plot.
manual_left_top_level2 <- c("GP14", "GP16", "GP36", "GP122", "GP151", "GP21")
manual_right_edge_level2 <- c(
  "GP37",
  "GP6",
  "GP23",
  "GP29",
  "GP26",
  "GP177",
  "GP11"
)
manual_level2_programs <- c(manual_left_top_level2, manual_right_edge_level2)

label_above_auto <- label_above %>%
  filter(!Factor %in% manual_level2_programs)
label_below_auto <- label_below %>%
  filter(!Factor %in% manual_level2_programs)

manual_level2_df <- df %>%
  filter(Factor %in% manual_level2_programs) %>%
  mutate(
    label_above_cand = mapply(
      manual_cats_label,
      Factor,
      annotation_Level2,
      MoreArgs = list(
        auc_matrix = level_2_AUC,
        positive_mask = level_2_AUC_positive
      )
    ),
    label_below_cand = mapply(
      manual_cats_label,
      Factor,
      annotation_Organ,
      MoreArgs = list(
        auc_matrix = organ_AUC,
        positive_mask = organ_AUC_positive_l2
      )
    ),
    manual_label = ifelse(residual >= 0, label_above_cand, label_below_cand)
  )
manual_above_df <- manual_level2_df %>%
  filter(residual >= 0)
manual_below_df <- manual_level2_df %>%
  filter(residual < 0)

base_plot <- ggplot(df, aes(Max_AUC_Organ, Max_AUC_Level1)) +
  geom_point(alpha = 0.3, size = 1.8) +
  geom_point(
    data = label_above_auto,
    color = "#1f78b4",
    alpha = 0.8,
    size = 1.8
  ) +
  geom_point(
    data = label_below_auto,
    color = "#e31a1c",
    alpha = 0.8,
    size = 1.8
  ) +
  geom_point(
    data = manual_above_df,
    color = "#1f78b4",
    alpha = 0.8,
    size = 1.8
  ) +
  geom_point(
    data = manual_below_df,
    color = "#e31a1c",
    alpha = 0.8,
    size = 1.8
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(0.46, 1.04), ylim = c(0.5, 1.02), expand = FALSE) +
  labs(
    x = "Max AUC (Organ Simplified)",
    y = "Max AUC (Level-2)",
    title = "Max AUC: Organ vs Level-2"
  ) +
  theme_minimal(base_size = 13)

p_all <- base_plot +
  geom_text_repel(
    data = label_above_auto,
    aes(label = label_text),
    color = "#1f78b4",
    size = 2.5,
    lineheight = 0.85,
    direction = "y",
    nudge_x = label_above_auto$nudge_x,
    segment.color = "#1f78b4",
    arrow = arrow(length = unit(0.008, "npc"), type = "closed", angle = 20),
    force = 3,
    force_pull = 0.1,
    box.padding = 0.4,
    point.padding = 0.15,
    max.time = 10,
    max.iter = 2e4,
    max.overlaps = 20,
    min.segment.length = 0.01,
    segment.alpha = 0.7
  ) +
  geom_text_repel(
    data = label_below_auto,
    aes(label = label_text),
    color = "#e31a1c",
    size = 2.5,
    lineheight = 0.85,
    direction = "y",
    nudge_x = label_below_auto$nudge_x,
    segment.color = "#e31a1c",
    arrow = arrow(length = unit(0.008, "npc"), type = "closed", angle = 20),
    force = 3,
    force_pull = 0.1,
    box.padding = 0.4,
    point.padding = 0.15,
    max.time = 10,
    max.iter = 2e4,
    max.overlaps = 20,
    min.segment.length = 0.01,
    segment.alpha = 0.7
  ) +
  geom_text_repel(
    data = manual_above_df,
    aes(label = manual_label),
    color = "#1f78b4",
    size = 2.5,
    lineheight = 0.85,
    direction = "y",
    nudge_x = -0.035,
    segment.color = "#1f78b4",
    arrow = arrow(length = unit(0.008, "npc"), type = "closed", angle = 20),
    force = 3,
    force_pull = 0.1,
    box.padding = 0.4,
    point.padding = 0.15,
    max.time = 10,
    max.iter = 2e4,
    max.overlaps = 20,
    min.segment.length = 0.01,
    segment.alpha = 0.7
  ) +
  geom_text_repel(
    data = manual_below_df,
    aes(label = manual_label),
    color = "#e31a1c",
    size = 2.5,
    lineheight = 0.85,
    direction = "y",
    nudge_x = 0.035,
    segment.color = "#e31a1c",
    arrow = arrow(length = unit(0.008, "npc"), type = "closed", angle = 20),
    force = 3,
    force_pull = 0.1,
    box.padding = 0.4,
    point.padding = 0.15,
    max.time = 10,
    max.iter = 2e4,
    max.overlaps = 20,
    min.segment.length = 0.01,
    segment.alpha = 0.7
  )
p_all
ggsave(
  filename = paste0(figure_path, "AUC_Organ_vs_Level2.pdf"),
  plot = p_all,
  width = 8,
  height = 8,
  dpi = 300
)

# New version: only the 7 organ-specific GPs highlighted in red,
# labeled with their best Level-2 categories (AUC > 0.9 only).
# Organ labels omitted — already shown in the Level-1 figure.
# Top-left GPs (high Level-2 AUC, lower organ AUC) labeled in blue.
seven_gp_df <- df |>
  dplyr::filter(Factor %in% gps_of_interest) |>
  dplyr::mutate(
    label_text = sapply(
      Factor,
      top_cats_label,
      auc_matrix = level_2_AUC,
      positive_mask = level_2_AUC_positive,
      threshold = 0.9,
      n = 3
    )
  )

top_left_gps <- c("GP14", "GP36", "GP16", "GP151", "GP21", "GP122",
                   "GP2", "GP171", "GP5", "GP13")
top_left_df <- df |>
  dplyr::filter(Factor %in% top_left_gps) |>
  dplyr::mutate(
    label_text = sapply(
      Factor,
      top_cats_label,
      auc_matrix = level_2_AUC,
      positive_mask = level_2_AUC_positive,
      threshold = 0.9,
      n = 3
    )
  )

p_seven_gps <- ggplot(df, aes(Max_AUC_Organ, Max_AUC_Level1)) +
  geom_point(alpha = 0.2, size = 1.5, color = "grey60") +
  geom_point(data = top_left_df, color = "#1f78b4", size = 2.2, alpha = 0.9) +
  geom_text_repel(
    data = top_left_df,
    aes(label = label_text),
    color = "#1f78b4",
    lineheight = 0.85,
    size = 2.5,
    direction = "y",
    nudge_x = -0.1,
    segment.color = "#1f78b4",
    arrow = arrow(length = unit(0.008, "npc"), type = "closed", angle = 20),
    force = 4,
    force_pull = 0.05,
    box.padding = 0.5,
    point.padding = 0.15,
    max.time = 10,
    max.iter = 2e4,
    max.overlaps = 30,
    min.segment.length = 0.01,
    segment.alpha = 0.7
  ) +
  geom_point(data = seven_gp_df, color = "#e31a1c", size = 2.2, alpha = 0.9) +
  geom_text_repel(
    data = seven_gp_df,
    aes(label = label_text),
    color = "#e31a1c",
    size = 2.5,
    lineheight = 0.85,
    direction = "y",
    nudge_x = 0.18,
    xlim = c(1.0, NA),
    segment.color = "#e31a1c",
    arrow = arrow(length = unit(0.008, "npc"), type = "closed", angle = 20),
    force = 6,
    force_pull = 0.02,
    box.padding = 0.6,
    point.padding = 0.15,
    max.time = 10,
    max.iter = 2e4,
    max.overlaps = 30,
    min.segment.length = 0.01,
    segment.alpha = 0.7
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(0.46, 1.04), ylim = c(0.5, 1.02), expand = FALSE,
                  clip = "off") +
  labs(
    x = "Max AUC (Organ Simplified)",
    y = "Max AUC (Level-2)",
    title = "Max AUC: Organ vs Level-2 — organ-specific GPs"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.margin = margin(10, 80, 10, 80))

p_seven_gps
ggsave(
  filename = paste0(figure_path, "AUC_Organ_vs_Level2_seven_gps.pdf"),
  plot = p_seven_gps,
  width = 8,
  height = 8,
  dpi = 300
)


# ============================================================
# ROC curves: a single GP predicting a single organ
#   (CD4/CD8 healthy cells only)
#   -> ROC_GP37_mammary_gland_cd4cd8.pdf
#   -> ROC_GP6_skin_cd4cd8.pdf
#   -> ROC_GP26_submandibular_gland_cd4cd8.pdf
#   -> ROC_GP23_skin_cd4cd8.pdf
# ============================================================
plot_gp_roc <- function(
  gp,
  group,
  loading_mat,
  group_info,
  threshold = NULL, # NULL => use optimal threshold (best_method below)
  best_method = c("closest.topleft", "youden"),
  direction = "auto",
  base_size = 13,
  line_color = "#1f78b4",
  highlight_color = "#e31a1c"
) {
  if (!gp %in% colnames(loading_mat)) {
    stop(sprintf("GP '%s' not found in loading matrix.", gp))
  }
  if (!group %in% group_info) {
    stop(sprintf("Group '%s' not found in group_info.", group))
  }
  best_method <- match.arg(best_method)

  loading <- loading_mat[, gp]
  keep <- !(is.na(loading) | is.na(group_info))
  loading <- loading[keep]
  group_info <- group_info[keep]
  labels <- as.numeric(group_info == group)

  roc_obj <- pROC::roc(
    response = labels,
    predictor = loading,
    quiet = TRUE,
    direction = direction,
    levels = c(0, 1)
  )
  auc_val <- as.numeric(pROC::auc(roc_obj))

  roc_df <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  )
  roc_df <- roc_df[order(roc_df$FPR, roc_df$TPR), ]

  # Compute the (FPR, TPR) for the highlight threshold
  if (is.null(threshold)) {
    coords_pt <- pROC::coords(
      roc_obj,
      x = "best",
      best.method = best_method,
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    )
    label_prefix <- "Optimal threshold"
  } else {
    coords_pt <- pROC::coords(
      roc_obj,
      x = threshold,
      input = "threshold",
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    )
    label_prefix <- "Threshold"
  }
  thr_used <- as.numeric(coords_pt$threshold)[1]
  hl_df <- data.frame(
    FPR = 1 - as.numeric(coords_pt$specificity)[1],
    TPR = as.numeric(coords_pt$sensitivity)[1]
  )
  hl_label <- sprintf(
    "%s = %.3g\nTPR = %.2f, FPR = %.2f",
    label_prefix,
    thr_used,
    hl_df$TPR,
    hl_df$FPR
  )

  ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_line(color = line_color, linewidth = 0.9) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray60"
    ) +
    geom_point(
      data = hl_df,
      aes(x = FPR, y = TPR),
      color = highlight_color,
      size = 3,
      inherit.aes = FALSE
    ) +
    ggrepel::geom_label_repel(
      data = hl_df,
      aes(x = FPR, y = TPR, label = hl_label),
      color = highlight_color,
      fill = alpha("white", 0.8),
      label.size = 0,
      size = base_size / 3.5,
      nudge_x = 0.15,
      nudge_y = -0.1,
      segment.color = highlight_color,
      inherit.aes = FALSE
    ) +
    annotate(
      "text",
      x = 0.98,
      y = 0.04,
      label = sprintf("AUC = %.3f", auc_val),
      hjust = 1,
      vjust = 0,
      size = base_size / 3
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(
      x = "False Positive Rate",
      y = "True Positive Rate",
      title = sprintf("ROC: %s predicting %s", gp, group)
    ) +
    theme_minimal(base_size = base_size)
}

# CD4/CD8 T cells from the healthy reference
CD4_CD8 <- seurat_meta_filtered_no_thymocytes_healthy %>%
  filter(annotation_level1 %in% c("CD4", "CD8")) %>%
  rownames()

# ---- GP37 / mammary gland (CD4/CD8) ----
p_gp37_mammary_cd4cd8 <- plot_gp_roc(
  gp = "GP37",
  group = "mammary gland",
  loading_mat = L_pm_filtered[CD4_CD8, ],
  group_info = seurat_meta_filtered_no_thymocytes_healthy[
    CD4_CD8,
  ]$organ_simplified
)
p_gp37_mammary_cd4cd8
ggsave(
  filename = paste0(figure_path, "ROC_GP37_mammary_gland_cd4cd8.pdf"),
  plot = p_gp37_mammary_cd4cd8,
  width = 5,
  height = 5,
  dpi = 300
)

# ---- GP6 / skin (CD4/CD8) ----
p_gp6_skin_cd4cd8 <- plot_gp_roc(
  gp = "GP6",
  group = "skin",
  loading_mat = L_pm_filtered[CD4_CD8, ],
  group_info = seurat_meta_filtered_no_thymocytes_healthy[
    CD4_CD8,
  ]$organ_simplified
)
p_gp6_skin_cd4cd8
ggsave(
  filename = paste0(figure_path, "ROC_GP6_skin_cd4cd8.pdf"),
  plot = p_gp6_skin_cd4cd8,
  width = 5,
  height = 5,
  dpi = 300
)

# ---- GP26 / submandibular gland (CD4/CD8) ----
p_gp26_submandibular_gland_cd4cd8 <- plot_gp_roc(
  gp = "GP26",
  group = "submandibular gland",
  loading_mat = L_pm_filtered[CD4_CD8, ],
  group_info = seurat_meta_filtered_no_thymocytes_healthy[
    CD4_CD8,
  ]$organ_simplified
)
p_gp26_submandibular_gland_cd4cd8
ggsave(
  filename = paste0(figure_path, "ROC_GP26_submandibular_gland_cd4cd8.pdf"),
  plot = p_gp26_submandibular_gland_cd4cd8,
  width = 5,
  height = 5,
  dpi = 300
)

# ---- GP23 / skin (CD4/CD8) ----
p_gp23_skin_cd4cd8 <- plot_gp_roc(
  gp = "GP23",
  group = "skin",
  loading_mat = L_pm_filtered[CD4_CD8, ],
  group_info = seurat_meta_filtered_no_thymocytes_healthy[
    CD4_CD8,
  ]$organ_simplified
)
p_gp23_skin_cd4cd8
ggsave(
  filename = paste0(figure_path, "ROC_GP23_skin_cd4cd8.pdf"),
  plot = p_gp23_skin_cd4cd8,
  width = 5,
  height = 5,
  dpi = 300
)


# ============================================================
# Threshold-defined GP+ rates for organ- and level2-driven programs
#   Organ-higher GPs are shown on the level-1 x axis with organ thresholds.
#   Level2-higher GPs are shown on the organ x axis with level-2 thresholds.
# ============================================================
# For each annotation group X (e.g. level-1 or level-2):
#   In organ            = #(group=X & in organ & loading>thr) / #(group=X & in organ)
#   Reference (healthy) = #(group=X & loading>thr)            / #(group=X)
plot_gp_threshold_group_activation_rate <- function(
  gp,
  organ,
  threshold,
  loading_mat,
  organ_info,
  group_info,
  group_label = "Level-2",
  base_size = 13,
  min_in_organ = 10,
  group_colors = ZemmourLib::immgent_colors$level2,
  fallback_group_color = "grey60",
  reference = c("not_in_group", "not_in_organ")
) {
  reference <- match.arg(reference)

  if (!gp %in% colnames(loading_mat)) {
    stop(sprintf("GP '%s' not found in loading matrix.", gp))
  }
  if (!organ %in% organ_info) {
    stop(sprintf("Organ '%s' not found in organ_info.", organ))
  }

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
    n_in_organ = vapply(
      group_levels,
      function(l) sum(group_info == l & in_organ),
      integer(1)
    ),
    n_pos_in_organ = vapply(
      group_levels,
      function(l) sum(group_info == l & in_organ & positive),
      integer(1)
    )
  )

  if (reference == "not_in_group") {
    rate_df$n_ref <- vapply(
      group_levels,
      function(l) sum(group_info != l & in_organ),
      integer(1)
    )
    rate_df$n_pos_ref <- vapply(
      group_levels,
      function(l) sum(group_info != l & in_organ & positive),
      integer(1)
    )
    ref_label <- "Not in group (same organ)"
    title_vs <- sprintf("%s vs. same-organ non-group", organ)
  } else {
    rate_df$n_ref <- vapply(
      group_levels,
      function(l) sum(group_info == l & !in_organ),
      integer(1)
    )
    rate_df$n_pos_ref <- vapply(
      group_levels,
      function(l) sum(group_info == l & !in_organ & positive),
      integer(1)
    )
    ref_label <- "Not in organ (same group)"
    title_vs <- sprintf("%s vs. same-group non-organ", organ)
  }

  rate_df$rate_in_organ <- rate_df$n_pos_in_organ / rate_df$n_in_organ
  rate_df$rate_ref <- rate_df$n_pos_ref / rate_df$n_ref

  rate_df <- rate_df[rate_df$n_in_organ >= min_in_organ, , drop = FALSE]
  if (nrow(rate_df) == 0) {
    stop(sprintf(
      "No %s type has >= %d cells in '%s'.",
      group_label,
      min_in_organ,
      organ
    ))
  }

  long_df <- data.frame(
    group = rep(rate_df$group, 2),
    type = factor(
      rep(c("In organ", ref_label), each = nrow(rate_df)),
      levels = c("In organ", ref_label)
    ),
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
      group_label,
      fallback_group_color,
      paste(level_order[missing_colors], collapse = ", ")
    ))
  }

  alpha_vals <- c(1, 0.35)
  names(alpha_vals) <- c("In organ", ref_label)

  ggplot(long_df, aes(x = group, y = rate, fill = group, alpha = type)) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.75,
      color = "grey35",
      linewidth = 0.15
    ) +
    scale_fill_manual(values = fill_values, guide = "none") +
    scale_alpha_manual(
      values = alpha_vals,
      guide = guide_legend(override.aes = list(fill = "grey40"))
    ) +
    labs(
      x = sprintf("%s annotation", group_label),
      y = sprintf("Proportion of cells with %s > %.3g", gp, threshold),
      alpha = NULL,
      title = sprintf("%s+ rate by %s: %s", gp, group_label, title_vs),
      subtitle = sprintf(
        "threshold = %.3g; %s types with < %d cells in %s dropped",
        threshold,
        group_label,
        min_in_organ,
        organ
      )
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
}

save_threshold_plot <- function(
  plot,
  filename,
  width = 8,
  height = 5,
  dpi = 300
) {
  ggsave(
    filename = paste0(figure_path, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

# ============================================================
# Level2-higher GP+ rate: target level-2 vs. same-organ background
#   (Organ x axis; level-2-specific threshold)
#   -> GP6_threshold_level2_organ_rate_CD4_Y.pdf
#   -> GP23_threshold_level2_organ_rate_gdT_R.pdf
#   -> GP26_threshold_level2_organ_rate_CD8_Q.pdf
#   -> GP29_threshold_level2_organ_rate_gdT_I.pdf
#   -> GP177_threshold_level2_organ_rate_CD4_X.pdf
# ============================================================
plot_level2_higher_organ_rate <- function(
  gp,
  level2,
  filename,
  threshold = level_2_AUC_list$threshold[level2, gp],
  loading_mat = L_pm_filtered[
    rownames(seurat_meta_filtered_no_thymocytes_healthy),
  ],
  organ_info = seurat_meta_filtered_no_thymocytes_healthy$organ_simplified,
  level2_info = seurat_meta_filtered_no_thymocytes_healthy$annotation_level2,
  base_size = 13,
  min_in_level2_organ = 20,
  min_background_organ = 100,
  organ_colors = ZemmourLib::immgent_colors$organ_simplified,
  fallback_organ_color = "grey60"
) {
  if (!gp %in% colnames(loading_mat)) {
    stop(sprintf("GP '%s' not found in loading matrix.", gp))
  }
  if (!level2 %in% level2_info) {
    stop(sprintf("Level-2 annotation '%s' not found in level2_info.", level2))
  }

  threshold <- as.numeric(threshold)
  if (length(threshold) != 1 || !is.finite(threshold)) {
    stop(sprintf("Invalid level-2 threshold for %s / %s.", gp, level2))
  }

  loading <- loading_mat[, gp]
  keep <- !(is.na(loading) | is.na(organ_info) | is.na(level2_info))
  loading <- loading[keep]
  organ_info <- as.character(organ_info[keep])
  level2_info <- as.character(level2_info[keep])

  in_level2 <- level2_info == level2
  positive <- loading > threshold
  organ_levels <- sort(unique(organ_info))

  rate_df <- data.frame(
    organ = organ_levels,
    n_in_level2 = vapply(
      organ_levels,
      function(l) sum(organ_info == l & in_level2),
      integer(1)
    ),
    n_pos_in_level2 = vapply(
      organ_levels,
      function(l) sum(organ_info == l & in_level2 & positive),
      integer(1)
    ),
    n_background = vapply(
      organ_levels,
      function(l) sum(organ_info == l & !in_level2),
      integer(1)
    ),
    n_pos_background = vapply(
      organ_levels,
      function(l) sum(organ_info == l & !in_level2 & positive),
      integer(1)
    )
  )
  rate_df$rate_in_level2 <- rate_df$n_pos_in_level2 / rate_df$n_in_level2
  rate_df$rate_background <- rate_df$n_pos_background / rate_df$n_background

  rate_df <- rate_df[
    rate_df$n_in_level2 >= min_in_level2_organ &
      rate_df$n_background >= min_background_organ,
    ,
    drop = FALSE
  ]
  if (nrow(rate_df) == 0) {
    stop(sprintf(
      "No organ has >= %d %s cells and >= %d same-organ background cells.",
      min_in_level2_organ,
      level2,
      min_background_organ
    ))
  }

  long_df <- data.frame(
    organ = rep(rate_df$organ, 2),
    type = factor(
      rep(c("In level2", "Not in level2"), each = nrow(rate_df)),
      levels = c("In level2", "Not in level2")
    ),
    rate = c(rate_df$rate_in_level2, rate_df$rate_background)
  )

  organ_order <- rate_df$organ[order(rate_df$rate_in_level2, decreasing = TRUE)]
  long_df$organ <- factor(long_df$organ, levels = organ_order)

  fill_values <- organ_colors[as.character(organ_order)]
  missing_colors <- is.na(fill_values)
  if (any(missing_colors)) {
    fill_values[missing_colors] <- fallback_organ_color
    warning(sprintf(
      "Organ annotations missing from organ_colors and colored %s: %s",
      fallback_organ_color,
      paste(organ_order[missing_colors], collapse = ", ")
    ))
  }

  p <- ggplot(long_df, aes(x = organ, y = rate, fill = organ, alpha = type)) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.75,
      color = "grey35",
      linewidth = 0.15
    ) +
    scale_fill_manual(
      values = fill_values,
      guide = "none"
    ) +
    scale_alpha_manual(
      values = c(
        "In level2" = 1,
        "Not in level2" = 0.35
      ),
      guide = guide_legend(override.aes = list(fill = "grey40"))
    ) +
    labs(
      x = "Organ (simplified)",
      y = sprintf("Proportion of cells with %s > %.3g", gp, threshold),
      alpha = NULL,
      title = sprintf(
        "%s+ rate by organ: %s vs. same-organ non-level2",
        gp,
        level2
      ),
      subtitle = sprintf(
        paste0(
          "level-2 threshold = %.3g; organs with < %d %s cells ",
          "or < %d background cells dropped"
        ),
        threshold,
        min_in_level2_organ,
        level2,
        min_background_organ
      )
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  save_threshold_plot(
    plot = p,
    filename = filename,
    width = 8,
    height = 5
  )
  p
}

plot_best_level2_higher_organ_rate <- function(gp) {
  level2 <- rownames(level_2_AUC_list$auc)[which.max(level_2_AUC_list$auc[,
    gp
  ])]
  level2_filename <- gsub("[^A-Za-z0-9]+", "_", level2)
  plot_level2_higher_organ_rate(
    gp = gp,
    level2 = level2,
    filename = sprintf(
      "%s_threshold_level2_organ_rate_%s.pdf",
      gp,
      level2_filename
    )
  )
}

# clearly level2
p_gp13_level2_organ <- plot_best_level2_higher_organ_rate("GP13")

# hard to see if level2 or organ
p_gp29_level2_organ <- plot_best_level2_higher_organ_rate("GP29")


# ============================================================
# Organ-higher GP+ rate: in selected organ vs. healthy reference
#   (Level-1 x axis; organ-specific threshold)
#   -> GP37_threshold_level1_rate_mammary_gland.pdf
#   -> GP49_threshold_level1_rate_mammary_gland.pdf
#   -> GP11_threshold_level1_rate_placenta.pdf
#   -> GP166_threshold_level1_rate_placenta.pdf
# ============================================================
plot_organ_higher_level1_rate <- function(gp, organ, filename) {
  p <- plot_gp_threshold_group_activation_rate(
    gp = gp,
    organ = organ,
    threshold = organ_AUC_list$threshold[organ, gp],
    min_in_organ = 100,
    loading_mat = L_pm_filtered[
      rownames(seurat_meta_filtered_no_thymocytes_healthy),
    ],
    organ_info = seurat_meta_filtered_no_thymocytes_healthy$organ_simplified,
    group_info = seurat_meta_filtered_no_thymocytes_healthy$annotation_level1,
    group_label = "Level-1",
    group_colors = ZemmourLib::immgent_colors$level1,
    reference = "not_in_organ"
  )
  save_threshold_plot(
    plot = p,
    filename = filename,
    width = 8,
    height = 5
  )
  p
}
p_gp37_level1 <- plot_organ_higher_level1_rate(
  gp = "GP37",
  organ = "mammary gland",
  filename = "GP37_threshold_level1_rate_mammary_gland.pdf"
)
p_gp11_level1 <- plot_organ_higher_level1_rate(
  gp = "GP11",
  organ = "placenta",
  filename = "GP11_threshold_level1_rate_placenta.pdf"
)


# ============================================================
# Swarm plot: per-organ AUC distribution across GPs
#   -> AUC_swarm_organ.pdf
# ============================================================
library(tidyr)

organ_AUC_sw <- organ_AUC_list$auc
organ_small_count_sw <- table(
  seurat_meta_filtered_no_thymocytes_healthy$organ_simplified
)
organ_small_count_sw <- names(organ_small_count_sw[organ_small_count_sw < 100])
organ_AUC_sw <- organ_AUC_sw[
  !rownames(organ_AUC_sw) %in% organ_small_count_sw,
]

# Mean loading per organ and overall mean (to determine prediction direction)
meta_sw <- seurat_meta_filtered_no_thymocytes_healthy
L_sw <- L_pm_filtered[rownames(meta_sw), ]
organ_levels_sw <- rownames(organ_AUC_sw)
mean_loading_by_organ <- t(sapply(organ_levels_sw, function(org) {
  cells_in_org <- which(meta_sw$organ_simplified == org)
  colMeans(L_sw[cells_in_org, , drop = FALSE], na.rm = TRUE)
}))
overall_loading_sw <- colMeans(L_sw, na.rm = TRUE)

# Sort organs by max AUC across GPs (descending)
organ_order <- names(sort(apply(organ_AUC_sw, 1, max), decreasing = TRUE))

organ_AUC_long <- as.data.frame(organ_AUC_sw)
organ_AUC_long$Organ <- rownames(organ_AUC_long)
organ_AUC_long <- tidyr::pivot_longer(
  organ_AUC_long,
  -Organ,
  names_to = "GP",
  values_to = "AUC"
)

mean_loading_long <- as.data.frame(mean_loading_by_organ)
mean_loading_long$Organ <- rownames(mean_loading_long)
mean_loading_long <- tidyr::pivot_longer(
  mean_loading_long,
  -Organ,
  names_to = "GP",
  values_to = "Mean_Loading"
)

organ_AUC_long <- organ_AUC_long %>%
  left_join(mean_loading_long, by = c("Organ", "GP")) %>%
  left_join(
    data.frame(
      GP = names(overall_loading_sw),
      Overall_Loading = unname(overall_loading_sw)
    ),
    by = "GP"
  ) %>%
  mutate(
    is_positive = Mean_Loading >= Overall_Loading,
    is_top = AUC >= 0.9 & is_positive,
    label_GP = ifelse(is_top, as.character(GP), ""),
    plot_shape = ifelse(is_top, 21L, 16L)
  )

# Keep only positively predictive GPs; organs with at least one GP AUC >= 0.9
organs_with_top <- organ_AUC_long %>%
  filter(is_top) %>%
  pull(Organ) %>%
  unique()

organ_AUC_plot <- organ_AUC_long %>%
  filter(is_positive, Organ %in% organs_with_top) %>%
  mutate(
    Organ = factor(
      as.character(Organ),
      levels = intersect(organ_order, organs_with_top)
    )
  )

pos <- position_jitter(width = 0.2, height = 0, seed = 42)

p_swarm_organ <- ggplot(
  organ_AUC_plot,
  aes(x = Organ, y = AUC, color = Organ)
) +
  geom_hline(
    yintercept = 0.9,
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.8
  ) +
  geom_jitter(
    aes(
      shape = I(plot_shape),
      size = I(ifelse(is_top, 2.5, 1.5)),
      alpha = I(ifelse(is_top, 1, 0.5)),
      stroke = I(ifelse(is_top, 1.2, 0)),
      fill = I(ifelse(is_top, "white", "transparent"))
    ),
    position = pos
  ) +
  geom_text_repel(
    aes(label = label_GP),
    position = pos,
    size = 3.5,
    color = "black",
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  labs(
    title = "GP Predictive Performance (AUC) by Organ",
    subtitle = "Positively predictive GPs only; labeled if AUC ≥ 0.9; dashed line = AUC 0.9",
    x = "Organ",
    y = "Area Under the Curve (AUC)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_swarm_organ
ggsave(
  filename = paste0(figure_path, "AUC_swarm_organ.pdf"),
  plot = p_swarm_organ,
  width = 10,
  height = 6,
  dpi = 300
)

# ============================================================
# Mean loading per organ for a single GP, in a chosen cell subset
#   (e.g. focus on healthy CD4/CD8); bars sorted high -> low
# ============================================================
plot_gp_mean_loading_by_organ <- function(
  gp,
  cells,
  loading_mat,
  organ_info, # named vector: organ_simplified per cell (names = cell ids)
  relative = FALSE, # if TRUE, scale bars so the highest organ = 1
  mark_largest_drop = TRUE, # mark the largest drop between consecutive bars
  base_size = 13,
  organ_colors = ZemmourLib::immgent_colors$organ_simplified
) {
  if (!gp %in% colnames(loading_mat)) {
    stop(sprintf("GP '%s' not found in loading matrix.", gp))
  }

  cells <- intersect(cells, rownames(loading_mat))
  cells <- intersect(cells, names(organ_info))

  loading <- loading_mat[cells, gp]
  organ <- organ_info[cells]
  keep <- !(is.na(loading) | is.na(organ))
  loading <- loading[keep]
  organ <- organ[keep]

  df <- data.frame(organ = organ, loading = loading) %>%
    group_by(organ) %>%
    summarise(
      n = dplyr::n(),
      mean_loading = mean(loading),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_loading)) %>%
    mutate(organ = factor(organ, levels = organ))

  if (relative) {
    df$value <- df$mean_loading / max(df$mean_loading)
    y_lab <- sprintf("Relative mean %s loading (max = 1)", gp)
  } else {
    df$value <- df$mean_loading
    y_lab <- sprintf("Mean %s loading", gp)
  }

  # Largest consecutive drop: boundary between bar i and i+1 with the biggest
  # fall in sorted loading -> how cleanly this GP separates its top organ(s).
  drops <- -diff(df$value)
  i_drop <- if (length(drops) > 0) which.max(drops) else integer(0)

  p <- ggplot(df, aes(x = organ, y = value, fill = organ)) +
    geom_col() +
    scale_fill_manual(values = organ_colors, guide = "none") +
    labs(
      x = "Organ (simplified)",
      y = y_lab,
      title = sprintf("Mean %s loading per organ", gp)
    ) +
    theme_minimal(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (mark_largest_drop && length(i_drop) == 1) {
    p <- p +
      geom_vline(
        xintercept = i_drop + 0.5,
        linetype = "dashed",
        color = "grey30"
      ) +
      annotate(
        "text",
        x = i_drop + 0.5,
        y = max(df$value),
        label = sprintf("largest drop = %.2f", drops[i_drop]),
        hjust = -0.05,
        vjust = 1,
        size = base_size / 3.5,
        color = "grey30"
      )
  }
  p
}

# Organ lookup aligned to cell ids (so we can subset by arbitrary cells)
organ_by_cell <- setNames(
  seurat_meta_filtered_no_thymocytes_healthy$organ_simplified,
  rownames(seurat_meta_filtered_no_thymocytes_healthy)
)

plot_gp_mean_loading_by_organ(
  gp = "GP37",
  cells = CD4_CD8,
  loading_mat = L_pm_filtered,
  organ_info = organ_by_cell
)

plot_gp_mean_loading_by_organ(
  gp = "GP11",
  cells = CD4_CD8,
  loading_mat = L_pm_filtered,
  organ_info = organ_by_cell
)

plot_gp_mean_loading_by_organ(
  gp = "GP6",
  cells = CD4_CD8,
  loading_mat = L_pm_filtered,
  organ_info = organ_by_cell
)

plot_gp_mean_loading_by_organ(
  gp = "GP47",
  cells = CD4_CD8,
  relative = FALSE,
  loading_mat = L_pm_filtered,
  organ_info = organ_by_cell
)


# ============================================================
# Side-by-side boxplot: GP loading for CD4.Y vs. non-CD4.Y,
#   split by in-organ vs. not-in-organ
# ============================================================
plot_gp_loading_boxplot <- function(
  gp,
  organ,
  level2, # the focal level-2 group (e.g. "CD4.Y")
  loading_mat,
  organ_info, # named vector: organ per cell
  level2_info, # named vector: level-2 annotation per cell
  base_size = 13,
  level2_colors = ZemmourLib::immgent_colors$level2,
  fallback_color = "grey60",
  max_cells_per_group = 5000 # downsample for speed / overplotting
) {
  cells <- intersect(names(organ_info), rownames(loading_mat))
  cells <- intersect(cells, names(level2_info))

  loading <- loading_mat[cells, gp]
  organ_vec <- organ_info[cells]
  level2_vec <- as.character(level2_info[cells])

  keep <- !(is.na(loading) | is.na(organ_vec) | is.na(level2_vec))
  loading <- loading[keep]
  organ_vec <- organ_vec[keep]
  level2_vec <- level2_vec[keep]

  x_group <- ifelse(level2_vec == level2, level2, paste0("non-", level2))
  organ_grp <- ifelse(organ_vec == organ, organ, paste0("non-", organ))

  plot_df <- data.frame(
    loading = loading,
    x_group = factor(x_group, levels = c(level2, paste0("non-", level2))),
    organ_grp = factor(organ_grp, levels = c(organ, paste0("non-", organ)))
  )

  # downsample each of the 4 combinations
  set.seed(42)
  plot_df <- plot_df |>
    dplyr::group_by(x_group, organ_grp) |>
    dplyr::slice_sample(n = max_cells_per_group, replace = FALSE) |>
    dplyr::ungroup()

  # colors: focal level2 gets its palette color, non- gets grey
  focal_color <- if (
    !is.null(level2_colors) && level2 %in% names(level2_colors)
  ) {
    level2_colors[[level2]]
  } else {
    fallback_color
  }
  fill_vals <- c(focal_color, "grey70")
  names(fill_vals) <- c(organ, paste0("non-", organ))

  ggplot(plot_df, aes(x = x_group, y = loading, fill = organ_grp)) +
    geom_boxplot(
      outlier.size = 0.3,
      outlier.alpha = 0.3,
      linewidth = 0.4,
      width = 0.6,
      position = position_dodge(width = 0.75)
    ) +
    scale_fill_manual(values = fill_vals) +
    ggpubr::stat_compare_means(
      aes(group = organ_grp),
      method = "wilcox.test",
      label = "p.signif",
      label.y.npc = 0.95,
      hide.ns = FALSE
    ) +
    labs(
      x = "Level-2 annotation",
      y = sprintf("%s loading", gp),
      fill = "Organ",
      title = sprintf(
        "%s loading: %s vs. non-%s  |  %s vs. non-%s",
        gp,
        level2,
        level2,
        organ,
        organ
      )
    ) +
    theme_minimal(base_size = base_size) +
    theme(legend.position = "top")
}

# ---- GP6 / skin: CD4.Y vs non-CD4.Y, skin vs non-skin ----
level2_by_cell <- setNames(
  seurat_meta_filtered_no_thymocytes_healthy$annotation_level2,
  rownames(seurat_meta_filtered_no_thymocytes_healthy)
)

p_gp6_skin_cd4y_boxplot <- plot_gp_loading_boxplot(
  gp = "GP6",
  organ = "skin",
  level2 = "CD4.Y",
  loading_mat = L_pm_filtered,
  organ_info = organ_by_cell,
  level2_info = level2_by_cell
)
p_gp6_skin_cd4y_boxplot
ggsave(
  filename = paste0(figure_path, "GP6_skin_CD4Y_boxplot.pdf"),
  plot = p_gp6_skin_cd4y_boxplot,
  width = 6,
  height = 5,
  dpi = 300
)


p_gp177_skin_cd4x_boxplot <- plot_gp_loading_boxplot(
  gp = "GP177",
  organ = "skin",
  level2 = "CD4.X",
  loading_mat = L_pm_filtered,
  organ_info = organ_by_cell,
  level2_info = level2_by_cell
)
p_gp177_skin_cd4x_boxplot
ggsave(
  filename = paste0(figure_path, "GP177_skin_CD4X_boxplot.pdf"),
  plot = p_gp177_skin_cd4x_boxplot,
  width = 6,
  height = 5,
  dpi = 300
)


library(ggalluvial)
# ============================================================
# Alluvial: Organ → Level-2  (GP+ cells only, organ-specific GPs)
#   For each of the 7 GPs, select GP+ cells (loading > organ threshold),
#   subsample to equal size, then show where those cells come from (organ,
#   left) and what cell type they are (level-2, right). Flows colored by GP.
#   A cell can appear more than once if it is GP+ for multiple programs.
#   -> alluvial_gp_organ_level2.pdf
# ============================================================
gps_of_interest <- c("GP3", "GP6", "GP11", "GP26", "GP29", "GP37", "GP177")
# For each GP, use the threshold from its best-predicting organ
best_organ_per_gp <- organ_AUC_max_name[gps_of_interest]
gp_thresholds <- mapply(
  function(gp, organ) {
    organ_AUC_list$threshold[organ, gp]
  },
  gps_of_interest,
  best_organ_per_gp
)
names(gp_thresholds) <- gps_of_interest
print(data.frame(
  gp = gps_of_interest,
  best_organ = unname(best_organ_per_gp),
  threshold = unname(gp_thresholds)
))

# Switch between organ-specific thresholds (default) and a fixed threshold
gp_thresholds_use <- gp_thresholds # organ-specific
# gp_thresholds_use <- setNames(rep(0.1, length(gps_of_interest)), gps_of_interest) # fixed 0.1
# Select GP+ cells for each GP; subsample to n_cap_gp to keep groups balanced
n_cap_gp <- 300
set.seed(42)
alluvial_rows <- lapply(gps_of_interest, function(gp) {
  positive_idx <- L_healthy[, gp] > gp_thresholds_use[gp]
  meta_pos <- seurat_meta_filtered_no_thymocytes_healthy[positive_idx, ]
  df <- data.frame(
    gp_program = gp,
    organ = meta_pos$organ_simplified,
    level2 = meta_pos$annotation_level2,
    stringsAsFactors = FALSE
  )
  if (nrow(df) > n_cap_gp) {
    df <- dplyr::slice_sample(df, n = n_cap_gp)
  }
  df
})

count_df <- do.call(rbind, alluvial_rows) |>
  dplyr::count(organ, gp_program, level2, name = "n") |>
  dplyr::filter(!is.na(organ), !is.na(level2), n >= 5)

organ_order <- count_df |>
  dplyr::summarise(total = sum(n), .by = organ) |>
  dplyr::arrange(dplyr::desc(total)) |>
  dplyr::pull(organ)
level2_order <- count_df |>
  dplyr::summarise(total = sum(n), .by = level2) |>
  dplyr::arrange(dplyr::desc(total)) |>
  dplyr::pull(level2)

count_df <- count_df |>
  dplyr::mutate(
    organ = factor(organ, levels = rev(organ_order)),
    gp_program = factor(gp_program, levels = rev(gps_of_interest)),
    level2 = factor(level2, levels = rev(level2_order))
  )

# Color each GP by its best-predicted organ
gp_colors <- ZemmourLib::immgent_colors$organ_simplified[unname(
  best_organ_per_gp
)]
gp_colors[is.na(gp_colors)] <- "grey60"
names(gp_colors) <- gps_of_interest

p_alluvial_gp_organ <- ggplot(
  count_df,
  aes(axis1 = organ, axis2 = gp_program, axis3 = level2, y = n)
) +
  ggalluvial::geom_alluvium(
    aes(fill = gp_program),
    width = 1 / 4,
    alpha = 0.6,
    knot.pos = 0.4
  ) +
  ggalluvial::geom_stratum(
    width = 1 / 4,
    fill = "grey92",
    color = "grey50",
    linewidth = 0.3
  ) +
  ggplot2::geom_text(
    stat = ggalluvial::StatStratum,
    aes(label = after_stat(stratum)),
    size = 3,
    angle = 90
  ) +
  scale_fill_manual(values = gp_colors, guide = "none") +
  scale_x_discrete(
    limits = c("Organ", "GP", "Level-2"),
    expand = c(0.12, 0.12)
  ) +
  labs(
    y = "Number of GP+ cells",
    title = "GP+ cells: organ origin and cell type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_flip()

p_alluvial_gp_organ

#

ggsave(
  filename = paste0(figure_path, "alluvial_gp_organ_level2.pdf"),
  plot = p_alluvial_gp_organ,
  width = 20,
  height = 10,
  dpi = 300
)


# ============================================================
# GP Decomposition: per-GP alluvial plots
#   For each GP: balanced sample from target organ vs all others,
#   then flow through GP+/GP- status to level-2 cluster.
#   Layout: Origin → GP Status → Level-2
#   -> gp_decomposition.pdf  (7 panels, one per GP)
# ============================================================

plot_gp_alluvial <- function(
  gp,
  n_cap = 3000,
  min_l2_cells = 100,
  min_l2_frac = Inf
) {
  tgt <- best_organ_per_gp[gp]
  thr <- gp_thresholds[gp]
  meta <- seurat_meta_filtered_no_thymocytes_healthy
  is_tgt <- meta$organ_simplified == tgt
  meta_tgt <- meta[is_tgt, ]
  meta_oth <- meta[!is_tgt, ]

  # Use all target organ cells (up to n_cap); match other group to same size
  n_use <- min(n_cap, nrow(meta_tgt))
  if (nrow(meta_tgt) > n_use) {
    meta_tgt <- meta_tgt[sample(nrow(meta_tgt), n_use), ]
  }
  if (nrow(meta_oth) > n_use) {
    meta_oth <- meta_oth[sample(nrow(meta_oth), n_use), ]
  }

  cell_df <- rbind(
    data.frame(
      origin = tgt,
      level2 = as.character(meta_tgt$annotation_level2),
      row.names = rownames(meta_tgt)
    ),
    data.frame(
      origin = "Other",
      level2 = as.character(meta_oth$annotation_level2),
      row.names = rownames(meta_oth)
    )
  )
  cell_df$gp_status <- ifelse(
    L_healthy[rownames(cell_df), gp] > thr,
    "GP+",
    "GP-"
  )

  # Show level-2 groups that exceed both the fraction and absolute count thresholds
  top_l2 <- cell_df |>
    dplyr::filter(!is.na(level2)) |>
    dplyr::count(level2, sort = TRUE) |>
    dplyr::mutate(frac = n / sum(n)) |>
    dplyr::filter(n > min_l2_cells | frac > min_l2_frac) |>
    dplyr::pull(level2)

  level2_pal <- ZemmourLib::immgent_colors$level2
  l2_colors <- c(
    setNames(
      ifelse(is.na(level2_pal[top_l2]), "grey55", level2_pal[top_l2]),
      top_l2
    ),
    "Other" = "grey82"
  )

  # GP+ rates for subtitle
  rate_tgt <- mean(cell_df$gp_status[cell_df$origin == tgt] == "GP+")
  rate_oth <- mean(cell_df$gp_status[cell_df$origin == "Other"] == "GP+")

  count_df <- cell_df |>
    dplyr::filter(!is.na(level2)) |>
    dplyr::mutate(
      level2_grp = factor(
        ifelse(level2 %in% top_l2, level2, "Other"),
        levels = c(top_l2, "Other")
      ),
      origin = factor(origin, levels = c(tgt, "Other")),
      gp_status = factor(gp_status, levels = c("GP+", "GP-"))
    ) |>
    dplyr::count(origin, gp_status, level2_grp, name = "n")

  ggplot(
    count_df,
    aes(axis1 = origin, axis2 = gp_status, axis3 = level2_grp, y = n)
  ) +
    ggalluvial::geom_alluvium(
      aes(fill = level2_grp),
      alpha = 0.55,
      width = 1 / 4,
      knot.pos = 0.4
    ) +
    ggalluvial::geom_stratum(
      width = 1 / 4,
      fill = "grey92",
      color = "grey50",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      stat = ggalluvial::StatStratum,
      aes(label = after_stat(stratum)),
      size = 2.5
    ) +
    scale_fill_manual(values = l2_colors, guide = "none") +
    scale_x_discrete(
      limits = c("Origin", "GP Status", "Level-2"),
      expand = c(0.12, 0.12)
    ) +
    labs(
      title = paste0(gp, "  |  best organ: ", tgt),
      subtitle = sprintf(
        "GP+ rate — target: %.0f%%  |  other: %.0f%%",
        100 * rate_tgt,
        100 * rate_oth
      ),
      y = "Cells"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8, color = "grey40")
    ) +
    coord_flip()
}

set.seed(42)
gp_alluvial_plots <- plot_gp_alluvial(gp = "GP3")
gp_alluvial_plots

# Diagnostic: mammary gland cell count and level-2 breakdown for GP37
mg_meta <- seurat_meta_filtered_no_thymocytes_healthy[
  seurat_meta_filtered_no_thymocytes_healthy$organ_simplified ==
    best_organ_per_gp["GP37"],
]
cat("Total mammary gland cells:", nrow(mg_meta), "\n")
print(sort(table(mg_meta$annotation_level2), decreasing = TRUE))

ggsave(
  filename = paste0(figure_path, "gp_decomposition.pdf"),
  plot = p_decomp,
  width = 12,
  height = 4 * length(gps_of_interest),
  dpi = 300
)

# ============================================================
# Heatmap: top positively regulated genes per GP
#   For the 7 organ-specific GPs, select the top N genes with
#   the highest positive loading. Gene selection is positive-only
#   (same logic as plot_factor_heatmap in Figure_CITEseq.R but
#   without the negative side). Loadings are column-scaled by
#   max(abs) so all GPs are on the same 0–1 scale.
#   -> organ_gp_gene_heatmap.pdf
# ============================================================

# Load and scale gene factor matrix (genes × GPs)
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F_pm_filtered) <- gsub("^F", "GP", colnames(F_pm_filtered))
D_scale <- diag(
  1 / apply(F_pm_filtered, 2, function(x) max(abs(x), na.rm = TRUE))
)
F_pm_scaled <- F_pm_filtered %*% D_scale
colnames(F_pm_scaled) <- colnames(F_pm_filtered)

# Select top N positive-loading genes per GP
n_top_genes <- 20
F_sub <- F_pm_scaled[, gps_of_interest, drop = FALSE]

min_loading <- 0.25
selected_genes <- lapply(gps_of_interest, function(gp) {
  vals <- F_sub[, gp]
  names(sort(vals[vals > min_loading], decreasing = TRUE))[seq_len(min(
    n_top_genes,
    sum(vals > min_loading)
  ))]
})
selected_genes <- unique(unlist(selected_genes))

# Diagonal gene ordering by dominant GP (highest loading within GP_orders)
# Used by both the heatmap and the dotplot
GP_orders <- c("GP37", "GP26", "GP6", "GP177", "GP3", "GP29", "GP11")
dominant_gp <- apply(F_sub[selected_genes, , drop = FALSE], 1, function(x) {
  GP_orders[which.max(x[GP_orders])]
})
dominant_loading <- mapply(
  function(g, gp) F_sub[g, gp],
  selected_genes,
  dominant_gp
)
gene_order_df <- data.frame(
  Gene = selected_genes,
  dominant_gp = factor(dominant_gp, levels = GP_orders),
  loading = dominant_loading,
  stringsAsFactors = FALSE
)
gene_order_df <- gene_order_df[
  order(gene_order_df$dominant_gp, -gene_order_df$loading),
]
heatmap_gene_order <- gene_order_df$Gene
organ_gp_genes <- heatmap_gene_order

# the cell ID of all the cells that we considered in this analysis
analysis_cells <- rownames(seurat_meta_filtered_no_thymocytes_healthy)


p_gp37_level1 <- plot_organ_higher_level1_rate(
  gp = "GP37",
  organ = "mammary gland",
  filename = "GP37_threshold_level1_rate_mammary_gland.pdf"
)

# Stacked barplot: same in-organ vs not-in-organ structure per level-1 lineage,
# but each bar is filled by level-2 proportion within GP+ cells.
# Height = GP+ rate; fill = level-2 contribution to that rate.
# Faceted by level-1; no legend (colors are for visual texture only).
{
  gp_sb <- "GP37"
  organ_sb <- "mammary gland"
  thr_sb <- organ_AUC_list$threshold[organ_sb, gp_sb]
  min_in_org <- 100

  meta_sb <- seurat_meta_filtered_no_thymocytes_healthy
  load_sb <- L_healthy[rownames(meta_sb), gp_sb]

  cell_df <- data.frame(
    level1 = meta_sb$annotation_level1,
    level2 = meta_sb$annotation_level2,
    group = ifelse(meta_sb$organ_simplified == organ_sb, organ_sb, "Other"),
    gp_pos = load_sb > thr_sb,
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(!is.na(level1), !is.na(level2))

  # Keep only level-1 lineages with >= min_in_org cells in target organ
  keep_l1 <- cell_df |>
    dplyr::filter(group == organ_sb) |>
    dplyr::count(level1) |>
    dplyr::filter(n >= min_in_org) |>
    dplyr::pull(level1)

  # Denominator: total cells per (level1, group)
  denom <- cell_df |>
    dplyr::filter(level1 %in% keep_l1) |>
    dplyr::count(level1, group, name = "n_total")

  # Numerator: GP+ cells per (level1, level2, group)
  numer <- cell_df |>
    dplyr::filter(level1 %in% keep_l1, gp_pos) |>
    dplyr::count(level1, level2, group, name = "n_pos")

  sb_df <- dplyr::left_join(numer, denom, by = c("level1", "group")) |>
    dplyr::mutate(
      rate_contrib = n_pos / n_total,
      group = factor(group, levels = c(organ_sb, "Other"))
    )

  level2_pal <- ZemmourLib::immgent_colors$level2

  # Use numeric x-axis so coordinate arithmetic works for arrow endpoints.
  # organ = x 1, Other = x 2; bar half-width = 0.4
  bar_w <- 0.4
  x_org <- 1L
  x_oth <- 2L
  x_text_org <- x_org - bar_w - 0.08 # text anchor left of organ bar
  x_text_oth <- x_oth + bar_w + 0.08 # text anchor right of other bar

  # Pre-stack the data manually so y coordinates are fully under our control —
  # avoids any mismatch with ggplot2's internal position_stack ordering.
  l2_levels <- sort(unique(sb_df$level2))
  sb_df$level2 <- factor(sb_df$level2, levels = l2_levels)
  sb_df$x_pos <- ifelse(sb_df$group == organ_sb, x_org, x_oth)

  sb_stacked <- sb_df |>
    dplyr::arrange(level1, group, level2) |>
    dplyr::group_by(level1, group) |>
    dplyr::mutate(
      ymax = cumsum(rate_contrib),
      ymin = ymax - rate_contrib,
      y_mid = (ymin + ymax) / 2
    ) |>
    dplyr::ungroup()

  # Labeled segments — organ bar labels go left, other bar labels go right
  label_df <- sb_stacked |>
    dplyr::filter(rate_contrib > 0.1) |>
    dplyr::mutate(
      x_arrow_end = ifelse(group == organ_sb, x_org - bar_w, x_oth + bar_w),
      x_label = ifelse(group == organ_sb, x_text_org, x_text_oth),
      hjust = ifelse(group == organ_sb, 1, 0)
    )

  p_gp37_stacked <- ggplot(sb_stacked) +
    geom_rect(
      aes(
        xmin = x_pos - bar_w,
        xmax = x_pos + bar_w,
        ymin = ymin,
        ymax = ymax,
        fill = level2
      )
    ) +
    geom_segment(
      data = label_df,
      aes(x = x_arrow_end, xend = x_label, y = y_mid, yend = y_mid),
      color = "grey40",
      linewidth = 0.3,
      arrow = arrow(length = unit(0.06, "cm"), type = "closed", ends = "first"),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = label_df,
      aes(x = x_label, y = y_mid, label = level2, hjust = hjust),
      size = 2.3,
      color = "black",
      inherit.aes = FALSE
    ) +
    facet_wrap(~level1, scales = "free_y", nrow = 1) +
    scale_x_continuous(
      breaks = c(x_org, x_oth),
      labels = c(organ_sb, "Other"),
      expand = c(0, 0.9)
    ) +
    scale_fill_manual(
      values = level2_pal,
      na.value = "grey60",
      guide = "none"
    ) +
    labs(
      title = paste0(
        gp_sb,
        " GP+ rate by level-2 composition — ",
        organ_sb,
        " vs other"
      ),
      x = NULL,
      y = "GP+ rate"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(face = "bold", size = 9)
    )

  p_gp37_stacked

  ggsave(
    filename = paste0(figure_path, "GP37_stacked_level2_mammary_gland.pdf"),
    plot = p_gp37_stacked,
    width = 14,
    height = 6,
    dpi = 300
  )
}

# Same stacked barplot for GP3
{
  gp_sb <- "GP3"
  organ_sb <- best_organ_per_gp["GP3"]
  thr_sb <- organ_AUC_list$threshold[organ_sb, gp_sb]
  min_in_org <- 100

  meta_sb <- seurat_meta_filtered_no_thymocytes_healthy
  load_sb <- L_healthy[rownames(meta_sb), gp_sb]

  cell_df <- data.frame(
    level1 = meta_sb$annotation_level1,
    level2 = meta_sb$annotation_level2,
    group = ifelse(meta_sb$organ_simplified == organ_sb, organ_sb, "Other"),
    gp_pos = load_sb > thr_sb,
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(!is.na(level1), !is.na(level2))

  keep_l1 <- cell_df |>
    dplyr::filter(group == organ_sb) |>
    dplyr::count(level1) |>
    dplyr::filter(n >= min_in_org) |>
    dplyr::pull(level1)

  denom <- cell_df |>
    dplyr::filter(level1 %in% keep_l1) |>
    dplyr::count(level1, group, name = "n_total")

  numer <- cell_df |>
    dplyr::filter(level1 %in% keep_l1, gp_pos) |>
    dplyr::count(level1, level2, group, name = "n_pos")

  sb_df <- dplyr::left_join(numer, denom, by = c("level1", "group")) |>
    dplyr::mutate(
      rate_contrib = n_pos / n_total,
      group = factor(group, levels = c(organ_sb, "Other"))
    )

  l2_levels <- sort(unique(sb_df$level2))
  sb_df$level2 <- factor(sb_df$level2, levels = l2_levels)
  sb_df$x_pos <- ifelse(sb_df$group == organ_sb, x_org, x_oth)

  sb_df_pos <- sb_df |>
    dplyr::arrange(level1, group, level2) |>
    dplyr::group_by(level1, group) |>
    dplyr::mutate(
      y_top = cumsum(rate_contrib),
      y_mid = y_top - rate_contrib / 2
    ) |>
    dplyr::ungroup()

  label_df <- sb_df_pos |>
    dplyr::filter(rate_contrib > 0.1) |>
    dplyr::mutate(
      x_pos = ifelse(group == organ_sb, x_org, x_oth),
      x_arrow_end = ifelse(group == organ_sb, x_org - bar_w, x_oth + bar_w),
      x_label = ifelse(group == organ_sb, x_text_org, x_text_oth),
      hjust = ifelse(group == organ_sb, 1, 0)
    )

  p_gp3_stacked <- ggplot(
    sb_df,
    aes(x = x_pos, y = rate_contrib, fill = level2)
  ) +
    geom_col(width = bar_w * 2) +
    geom_segment(
      data = label_df,
      aes(x = x_arrow_end, xend = x_label, y = y_mid, yend = y_mid),
      color = "grey40",
      linewidth = 0.3,
      arrow = arrow(length = unit(0.06, "cm"), type = "closed", ends = "first"),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = label_df,
      aes(x = x_label, y = y_mid, label = level2, hjust = hjust),
      size = 2.3,
      color = "black",
      inherit.aes = FALSE
    ) +
    facet_wrap(~level1, scales = "free_y", nrow = 1) +
    scale_x_continuous(
      breaks = c(x_org, x_oth),
      labels = c(organ_sb, "Other"),
      expand = c(0, 0.9)
    ) +
    scale_fill_manual(
      values = ZemmourLib::immgent_colors$level2,
      na.value = "grey60",
      guide = "none"
    ) +
    labs(
      title = paste0(gp_sb, " GP+ rate by level-2 — ", organ_sb, " vs other"),
      x = NULL,
      y = "GP+ rate"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(face = "bold", size = 9)
    )

  p_gp3_stacked

  ggsave(
    filename = paste0(figure_path, "GP3_stacked_level2_", organ_sb, ".pdf"),
    plot = p_gp3_stacked,
    width = 14,
    height = 6,
    dpi = 300
  )
}


library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(Matrix)
library(viridis)
library(cowplot)

expr <- readRDS(paste0(data_path, "shifted_log_counts_subset.rds"))
# rows = cells, cols = genes

tissue_order <- c(
  "mammary gland",
  "submandibular gland",
  "skin",
  "small intestine epi",
  "colon epi",
  "small intestine LP",
  "colon LP",
  "peritoneal cavity",
  "placenta",
  "liver",
  "lung",
  "kidney",
  "spleen",
  "LN"
)

# This plays the role of Seurat DotPlot(features = ...)
features <- rev(colnames(expr))

meta_use <- seurat_meta_filtered_no_thymocytes_healthy[
  rownames(expr),
  ,
  drop = FALSE
]

keep_cells <- meta_use$organ_simplified %in% tissue_order

expr_use <- expr[keep_cells, features, drop = FALSE]
meta_use <- meta_use[keep_cells, , drop = FALSE]

meta_use$organ_simplified <- factor(
  meta_use$organ_simplified,
  levels = tissue_order
)

# sparse-safe: returns both avg.exp and pct.exp in one pass (matches Seurat)
dot_stats <- function(mat) {
  avg_exp <- if (inherits(mat, "sparseMatrix")) {
    mat2 <- mat
    mat2@x <- expm1(mat2@x)
    Matrix::colMeans(mat2)
  } else {
    colMeans(expm1(mat))
  }
  list(
    avg.exp = as.numeric(avg_exp),
    pct.exp = as.numeric(Matrix::colMeans(mat > 0))
  )
}

tissues_present <- tissue_order[
  tissue_order %in% as.character(meta_use$organ_simplified)
]

dot_df <- map_dfr(tissues_present, function(tissue) {
  stats <- dot_stats(expr_use[
    meta_use$organ_simplified == tissue,
    features,
    drop = FALSE
  ])
  tibble(
    features.plot = features,
    id = tissue,
    avg.exp = stats$avg.exp,
    pct.exp = stats$pct.exp
  )
}) %>%
  mutate(
    avg.exp.scaled = log1p(avg.exp),
    pct.exp = pct.exp * 100,
    features.plot = factor(features.plot, levels = features),
    id = factor(id, levels = tissues_present)
  )

dot.scale <- 6

# Row-scaled version: Z-score computed across all displayed tissues
all_tissues <- unique(as.character(meta_use$organ_simplified))

dot_df_all <- map_dfr(all_tissues, function(tissue) {
  stats <- dot_stats(expr_use[
    meta_use$organ_simplified == tissue,
    features,
    drop = FALSE
  ])
  tibble(features.plot = features, id = tissue, avg.exp = stats$avg.exp)
})

global_stats <- dot_df_all |>
  dplyr::group_by(features.plot) |>
  dplyr::summarise(g_mean = mean(avg.exp), g_sd = sd(avg.exp), .groups = "drop")

dot_df_scaled <- dot_df |>
  dplyr::left_join(global_stats, by = "features.plot") |>
  dplyr::mutate(avg.exp.z = pmax(pmin((avg.exp - g_mean) / g_sd, 2.5), -2.5))

# Apply heatmap gene order to both dotplots; rev() because coord_flip() inverts
# factor level order (first level ends up at the bottom after flipping)
dot_df <- dot_df |>
  dplyr::mutate(
    features.plot = factor(features.plot, levels = rev(heatmap_gene_order))
  )
dot_df_scaled <- dot_df_scaled |>
  dplyr::mutate(
    features.plot = factor(features.plot, levels = rev(heatmap_gene_order))
  )

p <- ggplot(dot_df, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, dot.scale)) +
  scale_color_viridis_c(option = "C") +
  coord_flip() +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  guides(
    size = guide_legend(title = "Percent Expressed"),
    color = guide_colorbar(title = "Average Expression")
  ) +
  labs(x = "Features", y = "Identity")

pdf(
  paste0(figure_path, "gene_expression_dotplot_seurat_matched.pdf"),
  width = 8,
  height = 12,
  useDingbats = FALSE
)
print(p)
dev.off()

p_scaled <- ggplot(dot_df_scaled, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.z)) +
  scale_size(range = c(0, dot.scale)) +
  scale_color_distiller(
    palette = "RdBu",
    limits = c(-2.5, 2.5),
    direction = -1,
    name = "Avg Exp\n(Z-score)"
  ) +
  coord_flip() +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  guides(size = guide_legend(title = "Percent Expressed"))

pdf(
  paste0(figure_path, "gene_expression_dotplot_zscored.pdf"),
  width = 8,
  height = 12,
  useDingbats = FALSE
)
print(p_scaled)
dev.off()


# ============================================================
# Heatmap: top positively regulated genes per GP
#   Gene order: diagonal by dominant GP (highest loading)
#   -> organ_gp_gene_heatmap.pdf
# ============================================================

# Save genes and cells for downstream use
saveRDS(
  list(organ_gp_genes = organ_gp_genes, analysis_cells = analysis_cells),
  file = paste0(data_path, "organ_gp_genes_and_cells.rds")
)

plot_df_hm <- as.data.frame(F_sub[heatmap_gene_order, , drop = FALSE])
plot_df_hm$Gene <- rownames(plot_df_hm)
plot_df_hm <- tidyr::pivot_longer(
  plot_df_hm,
  cols = -Gene,
  names_to = "GP",
  values_to = "Loading"
)
plot_df_hm$GP <- factor(plot_df_hm$GP, levels = GP_orders)
plot_df_hm$Gene <- factor(plot_df_hm$Gene, levels = rev(heatmap_gene_order))

limit_hm <- max(abs(plot_df_hm$Loading), na.rm = TRUE)

p_gene_heatmap <- ggplot(plot_df_hm, aes(x = GP, y = Gene, fill = Loading)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-limit_hm, limit_hm),
    name = "Loading"
  ) +
  labs(title = "Top positive genes per organ GP", x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )

p_gene_heatmap

ggsave(
  filename = paste0(figure_path, "organ_gp_gene_heatmap.pdf"),
  plot = p_gene_heatmap,
  width = 6,
  height = 10,
  dpi = 300
)

# Side-by-side: gene vs tissue (dotplot, wider) | gene vs GP (heatmap, narrower)
# Gene order and labels appear only on the left; right panel suppresses y-axis.
library(patchwork)

p_combined <- (p_scaled +
  (p_gene_heatmap +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    ))) +
  plot_layout(widths = c(2, 1), guides = "collect") &
  theme(legend.position = "bottom")

pdf(
  paste0(figure_path, "gene_dotplot_heatmap_combined.pdf"),
  width = 10,
  height = 16,
  useDingbats = FALSE
)
print(p_combined)
dev.off()
