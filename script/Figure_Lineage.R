library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(fastTopics)
library(dplyr)
library(pheatmap)
library(scales)
library(scattermore)
library(Seurat)

#####################################################
#####################################################
####################################################
##### Defining directory and loading functions
####################################################
#####################################################
# data_path <- "../data/"
# code_path <- "../code/"
data_path <- "data/"
code_path <- "code/"
figure_path <- "figures/Figure2_Lineage/"
source(paste0(code_path, "ROC.R"))
source(paste0(code_path, "filtering_membership.R"))

#####################################################
#####################################################
#####################################################
### Loading Data
#####################################################
#####################################################
load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
seurat_meta <- readRDS(paste0(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
# flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
# protein_mat_normalized_lognorm <- readRDS(paste0(data_path, "protein_mat_normalized_lognorm.rds"))
# protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[rownames(L_pm_filtered) , "CD44"]
umap_result <- readRDS(paste0(data_path, "umap_result.rds"))
colnames(umap_result) <- c("UMAP_1", "UMAP_2")
umap_result <- umap_result[rownames(L_pm_filtered), ]
df_umap <- as.data.frame(umap_result)
level_1_AUC_list <- readRDS(
  file = paste0(data_path, "level_1_AUC_list_figure.rds")
)

# Rename K## to GP## for display consistency
colnames(L_pm_filtered) <- gsub("^K", "GP", colnames(L_pm_filtered))
colnames(level_1_AUC_list$auc) <- gsub(
  "^K",
  "GP",
  colnames(level_1_AUC_list$auc)
)


#####################################################
#####################################################
#####################################################
### Defining additional functions
#####################################################
#####################################################
# Identify Tukey outliers per group (matches geom_boxplot's outlier definition).
# Returns rows where value falls outside [Q1 - 1.5*IQR, Q3 + 1.5*IQR] within each
# group defined by `group_cols`. Used to plot outliers as a separate rasterized
# layer (smaller PDF), with `geom_boxplot(outlier.shape = NA)` for the boxes.
.tukey_outliers <- function(df, value_col, group_cols) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(
      .q1 = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
      .q3 = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
      .iqr = .q3 - .q1
    ) %>%
    filter(
      .data[[value_col]] < .q1 - 1.5 * .iqr |
        .data[[value_col]] > .q3 + 1.5 * .iqr
    ) %>%
    ungroup() %>%
    select(-.q1, -.q3, -.iqr)
}

# Lineage color palette (falls back to a default if ZemmourLib isn't loaded)
.lineage_colors <- function() {
  tryCatch(
    ZemmourLib::immgent_colors$level1,
    error = function(e) {
      c(
        CD4 = "blue",
        CD8 = "darkorange2",
        Treg = "deeppink",
        gdT = "chartreuse3",
        CD8aa = "darkorchid",
        Tz = "darkgoldenrod1",
        DN = "deepskyblue",
        DP = "red"
      )
    }
  )
}

# Unified GP-by-lineage swarm plot.
#
# Use cases:
#   - AUC or loading swarm with top-K labels (defaults)
#   - Loading-aware shape encoding (shape_by_loading = TRUE):
#       circles = up-regulated, triangles = down-regulated (Loading < Overall_Loading)
#   - Filter to up-regulated GPs only (filter_positive = TRUE),
#     with optional `forced_highlights` (named list lineage = c("GP##",...))
#     that bypass the filter and are always labeled.
#
# Down-regulated is defined as: Loading < Overall_Loading (strict).
plot_gp_swarm <- function(
  value_mat,
  loading_mat = NULL,
  overall_loading_vec = NULL,
  filter_positive = FALSE,
  shape_by_loading = FALSE,
  forced_highlights = NULL,
  top_k_labels = 3,
  threshold_line = NULL,
  threshold_color = "grey70",
  lineage_order = NULL,
  title = NULL,
  subtitle = NULL,
  y_label = "Value",
  bg_alpha = 0.6,
  italic_subtitle = FALSE,
  base_size = 12,
  seed = 42
) {
  df <- as.data.frame(as.table(as.matrix(value_mat)))
  colnames(df) <- c("Lineage", "GP", "Value")

  has_loading <- !is.null(loading_mat) && !is.null(overall_loading_vec)
  if (has_loading) {
    df_load <- as.data.frame(as.table(as.matrix(loading_mat)))
    colnames(df_load) <- c("Lineage", "GP", "Loading")
    df <- df %>%
      left_join(df_load, by = c("Lineage", "GP")) %>%
      left_join(
        data.frame(
          GP = names(overall_loading_vec),
          Overall_Loading = unname(overall_loading_vec)
        ),
        by = "GP"
      )
  }

  forced_keys <- character(0)
  if (!is.null(forced_highlights)) {
    forced_keys <- unlist(lapply(names(forced_highlights), function(lin) {
      paste(lin, forced_highlights[[lin]], sep = "::")
    }))
  }
  df$.key <- paste(df$Lineage, df$GP, sep = "::")

  if (filter_positive) {
    if (!has_loading) {
      stop("filter_positive requires loading_mat and overall_loading_vec")
    }
    df <- df %>% filter(Loading >= Overall_Loading | .key %in% forced_keys)
  }

  if (!is.null(lineage_order)) {
    existing <- intersect(lineage_order, unique(as.character(df$Lineage)))
    df$Lineage <- factor(df$Lineage, levels = existing)
  }

  df <- df %>%
    group_by(Lineage) %>%
    arrange(desc(Value)) %>%
    mutate(
      Rank = row_number(),
      is_forced = .key %in% forced_keys,
      is_top = (Rank <= top_k_labels) | is_forced,
      label_GP = ifelse(is_top, as.character(GP), "")
    ) %>%
    ungroup()

  if (has_loading && shape_by_loading) {
    df <- df %>%
      mutate(
        is_low_loading = Loading < Overall_Loading,
        plot_shape = case_when(
          is_top & !is_low_loading ~ 21L, # hollow circle
          is_top & is_low_loading ~ 24L, # hollow triangle
          !is_top & !is_low_loading ~ 16L, # solid circle
          !is_top & is_low_loading ~ 17L # solid triangle
        )
      )
  } else if (has_loading && filter_positive) {
    df <- df %>%
      mutate(
        is_low_loading = Loading < Overall_Loading,
        plot_shape = case_when(
          is_top & is_low_loading ~ 24L, # hollow triangle (forced down-reg)
          is_top & !is_low_loading ~ 21L, # hollow circle (top up-reg)
          TRUE ~ 16L # solid background dot
        )
      )
  } else {
    df$plot_shape <- ifelse(df$is_top, 21L, 16L)
  }

  pos <- position_jitter(width = 0.2, height = 0, seed = seed)

  p <- ggplot(df, aes(x = Lineage, y = Value, color = Lineage))

  if (!is.null(threshold_line)) {
    p <- p +
      geom_hline(
        yintercept = threshold_line,
        linetype = "dashed",
        color = threshold_color,
        linewidth = 0.8
      )
  }

  p +
    geom_jitter(
      aes(
        shape = I(plot_shape),
        size = I(ifelse(is_top, 2.5, 1.5)),
        alpha = I(ifelse(is_top, 1, bg_alpha)),
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
    scale_color_manual(values = .lineage_colors()) +
    labs(title = title, subtitle = subtitle, x = "Cell Lineage", y = y_label) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(
        hjust = 0.5,
        face = if (italic_subtitle) "italic" else "plain"
      ),
      panel.grid.major.x = element_blank(),
      legend.position = "none"
    )
}
plot_loadings_on_umap <- function(
  umap,
  loading,
  factor_num,
  size = 0.4,
  show_bg = TRUE,
  bg_alpha = 0.05,
  bg_color = "white",
  low_color = "lightgrey",
  high_color = "blue",
  clip_q = 0.99,
  lower_bound = 0.1,
  top_q = NULL
) {
  stopifnot(is.matrix(umap) || is.data.frame(umap))
  df <- as.data.frame(umap)
  stopifnot(all(c("UMAP_1", "UMAP_2") %in% colnames(df)))

  # Calibrate color scale on active cells only (ignore inactive cells below lower_bound)
  active_loadings <- loading[loading > lower_bound]

  if (length(active_loadings) > 0) {
    q_hi <- stats::quantile(active_loadings, probs = clip_q, na.rm = TRUE)
  } else {
    q_hi <- stats::quantile(loading, probs = clip_q, na.rm = TRUE)
  }

  loading_clip <- pmin(loading, q_hi)

  if (!is.null(top_q)) {
    thr <- stats::quantile(loading_clip, probs = top_q, na.rm = TRUE)
    loading_plot <- ifelse(loading_clip >= thr, loading_clip, NA_real_)
  } else {
    loading_plot <- loading_clip
  }

  df$loading_plot <- loading_plot

  # Z-order: plot higher loadings on top so they aren't masked by greys
  df <- df[order(is.na(df$loading_plot), df$loading_plot, decreasing = FALSE), ]

  p <- ggplot() +
    {
      if (show_bg) {
        scattermore::geom_scattermore(
          data = df,
          mapping = aes(x = UMAP_1, y = UMAP_2),
          pointsize = size,
          color = bg_color,
          alpha = bg_alpha
        )
      }
    } +
    scattermore::geom_scattermore(
      data = df[!is.na(df$loading_plot), ],
      mapping = aes(x = UMAP_1, y = UMAP_2, color = loading_plot),
      pointsize = size * 2,
      alpha = 1,
      na.rm = TRUE
    ) +
    scale_color_gradient(
      low = low_color,
      high = high_color,
      limits = c(0, q_hi),
      oob = scales::squish,
      breaks = scales::pretty_breaks(4),
      name = paste0("GP", factor_num)
    ) +
    coord_equal() +
    labs(
      title = paste0("Factor ", factor_num),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  return(p)
}


#####################################################
#####################################################
#####################################################
### Loading + AUC summary matrices and swarm plots
#####################################################
#####################################################
selected_lineage_in_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN")
# annotation_level1 is a factor; coerce to character so sapply propagates names
level_1_categories <- as.character(unique(
  seurat_meta_filtered$annotation_level1
))

# Per-lineage mean loading (rows: lineage, cols: GP)
mean_loading_matrix <- t(sapply(level_1_categories, function(cat) {
  cells_in_cat <- which(seurat_meta_filtered$annotation_level1 == cat)
  colMeans(L_pm_filtered[cells_in_cat, , drop = FALSE], na.rm = TRUE)
}))
mean_loading_matrix_selected <- mean_loading_matrix[selected_lineage_in_order, ]

# Overall mean loading per GP, computed across non-thymocyte cells only
non_thymocyte_cells <- which(
  seurat_meta_filtered$annotation_level1 != "thymocyte"
)
mean_loading_vec <- colMeans(
  L_pm_filtered[non_thymocyte_cells, , drop = FALSE],
  na.rm = TRUE
)

AUC_mat_selected <- level_1_AUC_list$auc[selected_lineage_in_order, ]

# (1) Mean loading swarm — y-axis is the lineage mean loading itself
p_loading_swarm <- plot_gp_swarm(
  mean_loading_matrix_selected,
  top_k_labels = 4,
  threshold_line = 0.1,
  threshold_color = "grey50",
  bg_alpha = 0.4,
  title = "GP Mean Loadings by Lineage",
  subtitle = "Top 4 highest loaded Gene Programs labeled per lineage",
  y_label = "Mean Loading"
)
ggsave(
  filename = paste0(figure_path, "Loading_swarm_level1.pdf"),
  plot = p_loading_swarm,
  width = 8,
  height = 6
)

# (2) AUC swarm
p_swarm <- plot_gp_swarm(
  AUC_mat_selected,
  top_k_labels = 4,
  threshold_line = 0.5,
  title = "GP Predictive Performance (AUC) by Lineage",
  subtitle = "Top 4 predictive Gene Programs labeled per lineage",
  y_label = "Area Under the Curve (AUC)"
)
ggsave(
  filename = paste0(figure_path, "AUC_swarm_level1.pdf"),
  plot = p_swarm,
  width = 8,
  height = 6
)

# (3) AUC + loading-aware shape encoding (triangles = down-regulated)
p_swarm_with_loading <- plot_gp_swarm(
  AUC_mat_selected,
  loading_mat = mean_loading_matrix,
  overall_loading_vec = mean_loading_vec,
  shape_by_loading = TRUE,
  top_k_labels = 4,
  threshold_line = 0.5,
  title = "GP Predictive Performance (AUC) by Lineage",
  subtitle = "Triangles indicate under-expression; Circles indicate over-expression",
  y_label = "Area Under the Curve (AUC)",
  italic_subtitle = TRUE
)
ggsave(
  filename = paste0(figure_path, "AUC_swarm_with_loading_level1.pdf"),
  plot = p_swarm_with_loading,
  width = 10,
  height = 6
)

# (4) AUC filtered to up-regulated GPs, with forced highlights
p_swarm_with_loading_positive <- plot_gp_swarm(
  AUC_mat_selected,
  loading_mat = mean_loading_matrix,
  overall_loading_vec = mean_loading_vec,
  filter_positive = TRUE,
  forced_highlights = list(
    gdT = c("GP22", "GP29", "GP3"),
    CD8aa = c("GP22", "GP29", "GP3"),
    DN = c("GP22", "GP29", "GP3")
  ),
  top_k_labels = 1,
  threshold_line = 0.5,
  title = "Up-regulated GP Predictive Performance (AUC)",
  subtitle = "Up-regulated GPs (Loading ≥ Overall Mean); triangles mark forced highlights that are down-regulated",
  y_label = "Area Under the Curve (AUC)",
  italic_subtitle = TRUE
)
ggsave(
  filename = paste0(figure_path, "AUC_swarm_with_loading_positive_level1.pdf"),
  plot = p_swarm_with_loading_positive,
  width = 10,
  height = 6
)


#####################################################
#####################################################
#####################################################
### Structure Plot by Lineage
#####################################################
#####################################################
set.seed(1234)
color_coding <- ZemmourLib::immgent_colors$level1
level1_category_to_factor <- c(
  "Treg" = "GP68",
  "CD8" = "GP58",
  "Tz" = "GP30",
  "DN" = "GP22",
  "CD8aa" = "GP29",
  "gdT" = "GP3"
)
fit2 <- L_pm_filtered[, level1_category_to_factor, drop = FALSE]
cell_type <- seurat_meta[rownames(fit2), "annotation_level1"]
cells <- which(cell_type == "CD4" | cell_type == "CD8")
cells <- sample(cells, 1e5)
cells <- sort(c(
  cells,
  which(
    cell_type != "CD4" &
      cell_type != "CD8" &
      cell_type != "thymocyte" &
      cell_type != "DP"
  )
))
structure_plot(
  fit2[cells, ],
  gap = 40,
  n = 10000,
  colors = color_coding[names(level1_category_to_factor)],
  grouping = cell_type[cells]
) +
  labs(y = "membership", color = "", fill = "") +
  guides(
    fill = guide_legend(nrow = 1),
    color = guide_legend(nrow = 1)
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )
ggsave(
  filename = paste0(figure_path, "structure_plot_level1.pdf"),
  width = 10,
  height = 6,
  dpi = 300
)


#####################################################
#####################################################
#####################################################
### MDE by Lineage
#####################################################
#####################################################
set.seed(1)
df_umap2 <- df_umap %>%
  tibble::rownames_to_column("cellID")
plot_df <- df_umap2 %>%
  inner_join(
    seurat_meta_filtered %>% select(cellID, annotation_level1),
    by = "cellID"
  ) %>%
  filter(annotation_level1 != "thymocyte")
max_total <- 1000000
min_per_group <- 300
cap_per_group <- 20000
group_sizes <- plot_df %>% count(annotation_level1, name = "n")
G <- nrow(group_sizes)
base_per_group <- ceiling(max_total / max(G, 1))
sample_plan <- group_sizes %>%
  mutate(
    n_take = pmin(n, pmax(min_per_group, pmin(cap_per_group, base_per_group)))
  )
plot_df_sub <- plot_df %>%
  group_by(annotation_level1) %>%
  group_modify(
    ~ dplyr::slice_sample(
      .x,
      n = sample_plan$n_take[
        sample_plan$annotation_level1 == .y$annotation_level1
      ]
    )
  ) %>%
  ungroup()
p_umap_level1_group <- ggplot(plot_df_sub, aes(x = UMAP_1, y = UMAP_2)) +
  scattermore::geom_scattermore(
    aes(color = annotation_level1),
    pointsize = 1.2
  ) +
  scale_color_manual(values = ZemmourLib::immgent_colors$level1) +
  coord_equal() +
  theme_classic() +
  labs(
    title = "MDE: Annotation Level 1",
    x = "MDE 1",
    y = "MDE 2",
    color = "Cell Type"
  ) +
  theme(
    legend.text = element_text(size = 10),
    legend.key.size = unit(1.5, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(
  filename = paste0(figure_path, "UMAP_level1_group.pdf"),
  width = 5,
  height = 5
)


#####################################################
#####################################################
#####################################################
### gdT-specific UMAP and GP loadings
#####################################################
#####################################################
library(Seurat)
gdT_seurat <- readRDS(paste0(
  data_path,
  "gdT_igt1_96_withtotalvi20260206_clean.Rds"
))
umap_result_gdT <- Embeddings(gdT_seurat, "mde_incremental")
colnames(umap_result_gdT) <- c("UMAP_1", "UMAP_2")

gdT_cells <- seurat_meta_filtered$cellID[
  seurat_meta_filtered$annotation_level1 == "gdT"
]
umap_gdT <- umap_result_gdT[gdT_cells, ]
df_umap_gdT <- as.data.frame(umap_gdT)
df_umap_gdT$cellID <- rownames(df_umap_gdT)
df_umap_gdT <- df_umap_gdT %>%
  inner_join(
    seurat_meta_filtered %>% select(cellID, annotation_level2),
    by = "cellID"
  )

p_umap_gdT <- ggplot(df_umap_gdT, aes(x = UMAP_1, y = UMAP_2)) +
  scattermore::geom_scattermore(
    aes(color = annotation_level2),
    pointsize = 1.2
  ) +
  coord_equal() +
  theme_classic() +
  labs(
    title = "UMAP of gdT cells",
    x = "UMAP 1",
    y = "UMAP 2",
    color = ""
  )
ggsave(
  filename = paste0(figure_path, "UMAP_gdT_level2.pdf"),
  plot = p_umap_gdT,
  width = 5,
  height = 4
)

df_umap_gdT$GP3_loading <- L_pm_filtered[df_umap_gdT$cellID, "GP3"]
df_umap_gdT$GP29_loading <- L_pm_filtered[df_umap_gdT$cellID, "GP29"]
df_umap_gdT$GP6_loading <- L_pm_filtered[df_umap_gdT$cellID, "GP6"]

p_umap_gdT_GP3 <- plot_loadings_on_umap(
  umap = df_umap_gdT[, c("UMAP_1", "UMAP_2")],
  loading = df_umap_gdT$GP3_loading,
  factor_num = 3,
  size = 0.6,
  bg_alpha = 0.1,
  bg_color = "grey90"
)
p_umap_gdT_GP29 <- plot_loadings_on_umap(
  umap = df_umap_gdT[, c("UMAP_1", "UMAP_2")],
  loading = df_umap_gdT$GP29_loading,
  factor_num = 29,
  size = 0.6,
  bg_alpha = 0.1,
  bg_color = "grey90"
)
p_umap_gdT_GP6 <- plot_loadings_on_umap(
  umap = df_umap_gdT[, c("UMAP_1", "UMAP_2")],
  loading = df_umap_gdT$GP6_loading,
  factor_num = 6,
  size = 0.6,
  bg_alpha = 0.1,
  bg_color = "grey90"
)

ggsave(
  filename = paste0(figure_path, "UMAP_gdT_GP3_loading.pdf"),
  plot = p_umap_gdT_GP3,
  width = 5,
  height = 4
)
ggsave(
  filename = paste0(figure_path, "UMAP_gdT_GP29_loading.pdf"),
  plot = p_umap_gdT_GP29,
  width = 5,
  height = 4
)
ggsave(
  filename = paste0(figure_path, "UMAP_gdT_GP6_loading.pdf"),
  plot = p_umap_gdT_GP6,
  width = 5,
  height = 4
)


### A boxplot figure comparing GP30 Loading across different subset of Tz cells:
### including iNKT (seurat_meta_filtered$iNKT)
### MAIT (seurat_meta_filtered$MAIT)
### other Tz cells (that are not iNKT nor MAIT)
### other T cells (that are not Tz)
GP30_df <- seurat_meta_filtered %>%
  select(cellID, annotation_level1, iNKT, MAIT) %>%
  filter(annotation_level1 != "thymocyte") %>%
  mutate(GP30_loading = L_pm_filtered[cellID, "GP30"]) %>%
  mutate(
    Group = case_when(
      annotation_level1 == "Tz" & iNKT ~ "iNKT",
      annotation_level1 == "Tz" & MAIT ~ "MAIT",
      annotation_level1 == "Tz" ~ "Other Tz",
      TRUE ~ "Other T cells"
    ),
    Group = factor(
      Group,
      levels = c("iNKT", "MAIT", "Other Tz", "Other T cells")
    )
  )

p_GP30_box <- ggplot(GP30_df, aes(x = Group, y = GP30_loading, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  ggrastr::rasterise(
    geom_point(
      data = .tukey_outliers(GP30_df, "GP30_loading", "Group"),
      size = 0.5,
      alpha = 0.3,
      show.legend = FALSE
    ),
    dpi = 300
  ) +
  scale_fill_manual(
    values = c(
      "iNKT" = "darkgoldenrod2",
      "MAIT" = "darkgoldenrod3",
      "Other Tz" = "darkgoldenrod1",
      "Other T cells" = "grey70"
    )
  ) +
  labs(
    title = "GP30 loading across Tz subsets and other T cells",
    x = NULL,
    y = "GP30 loading"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = paste0(figure_path, "GP30_boxplot_Tz_subsets.pdf"),
  plot = p_GP30_box,
  width = 6,
  height = 5
)

#### A Figure compare CD8 (resting) and CD8 (activated), and other T cells
#### loading of GP58.
non_CD8_lineages <- setdiff(selected_lineage_in_order, "CD8")
GP58_df <- seurat_meta_filtered %>%
  select(cellID, annotation_level1, annotation_level2_group) %>%
  mutate(GP58_loading = L_pm_filtered[cellID, "GP58"]) %>%
  mutate(
    Group = case_when(
      annotation_level1 == "CD8" & annotation_level2_group == "resting" ~
        "CD8 (Resting)",
      annotation_level1 == "CD8" & annotation_level2_group == "activated" ~
        "CD8 (Activated)",
      annotation_level1 %in% non_CD8_lineages ~ "Other T cells",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Group = factor(
      Group,
      levels = c("CD8 (Resting)", "CD8 (Activated)", "Other T cells")
    )
  )

p_GP58_box <- ggplot(GP58_df, aes(x = Group, y = GP58_loading, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  ggrastr::rasterise(
    geom_point(
      data = .tukey_outliers(GP58_df, "GP58_loading", "Group"),
      size = 0.5,
      alpha = 0.3,
      show.legend = FALSE
    ),
    dpi = 300
  ) +
  scale_fill_manual(
    values = c(
      "CD8 (Resting)" = "darkorange2",
      "CD8 (Activated)" = "orange",
      "Other T cells" = "grey70"
    )
  ) +
  labs(
    title = "GP58 loading across CD8 subsets and other T cells",
    x = NULL,
    y = "GP58 loading"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = paste0(figure_path, "GP58_boxplot_CD8_vs_other.pdf"),
  plot = p_GP58_box,
  width = 5,
  height = 5
)

### Plot compare protein CD8a and CD8b expression for
## between CD8 (resting) and CD8 (activated)
protein_mat_normalized_lognorm <- readRDS(paste0(
  data_path,
  "protein_mat_normalized_lognorm.rds"
))

cells_citeseq <- seurat_meta_filtered$cellID[seurat_meta_filtered$cite_seq]
CD8_citeseq_cells <- seurat_meta_filtered$cellID[
  seurat_meta_filtered$annotation_level1 == "CD8" &
    seurat_meta_filtered$cite_seq &
    seurat_meta_filtered$annotation_level2_group %in% c("resting", "activated")
]
CD8_citeseq_cells <- intersect(
  CD8_citeseq_cells,
  rownames(protein_mat_normalized_lognorm)
)

CD8_protein_df <- data.frame(
  cellID = CD8_citeseq_cells,
  Group = ifelse(
    seurat_meta_filtered[CD8_citeseq_cells, "annotation_level2_group"] ==
      "resting",
    "CD8 (Resting)",
    "CD8 (Activated)"
  ),
  CD8A = protein_mat_normalized_lognorm[CD8_citeseq_cells, "CD8A"],
  CD8B = protein_mat_normalized_lognorm[CD8_citeseq_cells, "CD8B"]
) %>%
  tidyr::pivot_longer(
    cols = c("CD8A", "CD8B"),
    names_to = "Protein",
    values_to = "Expression"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("CD8 (Resting)", "CD8 (Activated)"))
  )

p_CD8ab_box <- ggplot(
  CD8_protein_df,
  aes(x = Protein, y = Expression, fill = Group)
) +
  geom_boxplot(outlier.shape = NA) +
  ggrastr::rasterise(
    geom_point(
      data = .tukey_outliers(
        CD8_protein_df,
        "Expression",
        c("Protein", "Group")
      ),
      aes(group = Group),
      size = 0.5,
      alpha = 0.3,
      position = position_dodge(width = 0.75),
      show.legend = FALSE
    ),
    dpi = 300
  ) +
  scale_fill_manual(
    values = c(
      "CD8 (Resting)" = "darkorange2",
      "CD8 (Activated)" = "orange"
    )
  ) +
  labs(
    title = "CD8A / CD8B protein expression in CD8 cells",
    x = NULL,
    y = "Log-normalized protein expression",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

ggsave(
  filename = paste0(
    figure_path,
    "CD8ab_protein_boxplot_resting_vs_activated.pdf"
  ),
  plot = p_CD8ab_box,
  width = 5,
  height = 5
)

#####################################################
#####################################################
#####################################################
### gdT / CD8aa / DN: lineage UMAP and GP22/29/3 loadings
#####################################################
#####################################################
gcd_lineages <- c("gdT", "CD8aa", "DN")
# Exclude proliferating / miniverse subsets (annotation_level2_group)
gcd_cells <- seurat_meta_filtered$cellID[
  seurat_meta_filtered$annotation_level1 %in%
    gcd_lineages &
    !(seurat_meta_filtered$annotation_level2_group %in%
      c("proliferating", "miniverse"))
]
df_umap_gcd <- df_umap[gcd_cells, ]
df_umap_gcd$annotation_level1 <- factor(
  seurat_meta_filtered[gcd_cells, "annotation_level1"],
  levels = gcd_lineages
)

# (1) UMAP colored by lineage
p_umap_gcd_lineage <- ggplot(df_umap_gcd, aes(x = UMAP_1, y = UMAP_2)) +
  scattermore::geom_scattermore(
    aes(color = annotation_level1),
    pointsize = 1.2
  ) +
  scale_color_manual(values = .lineage_colors()[gcd_lineages]) +
  coord_equal() +
  theme_classic() +
  labs(
    title = "UMAP: gdT / CD8aa / DN",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cell Type"
  ) +
  theme(
    legend.text = element_text(size = 10),
    legend.key.size = unit(1.5, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(
  filename = paste0(figure_path, "UMAP_gdTCD8aaDN_lineage.pdf"),
  plot = p_umap_gcd_lineage,
  width = 5,
  height = 5
)

# (2-4) GP loadings on the gdT/CD8aa/DN-restricted UMAP
for (gp_num in c(22, 29, 3, 6)) {
  gp_name <- paste0("GP", gp_num)
  p_loading <- plot_loadings_on_umap(
    umap = df_umap_gcd[, c("UMAP_1", "UMAP_2")],
    loading = L_pm_filtered[gcd_cells, gp_name],
    factor_num = gp_num,
    size = 0.6,
    bg_alpha = 0.1,
    bg_color = "grey90"
  )
  ggsave(
    filename = paste0(
      figure_path,
      "UMAP_gdTCD8aaDN_",
      gp_name,
      "_loading.pdf"
    ),
    plot = p_loading,
    width = 5,
    height = 4
  )
}
