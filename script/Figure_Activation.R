# Bioconductor packages load first so their generics (rename, select, ...)
# don't mask dplyr/tidyverse versions used throughout the script.
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(fastTopics)
library(pheatmap)
library(scales)
library(scattermore)
library(tidygraph)
library(ggraph)
library(igraph)
library(stringr)
# dplyr / tidyverse last so dplyr::rename, dplyr::select win the masking war.
library(dplyr)
library(tidyverse)

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
figure_path <- "figures/Figure3_Activation/"
source(paste0(code_path, "ROC.R"))
source(paste0(code_path, "filtering_membership.R"))

#####################################################
#####################################################
#####################################################
### Loading Data
#####################################################
#####################################################
load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
# flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
colnames(L_pm_filtered) <- paste0("GP", seq_len(ncol(L_pm_filtered)))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F_pm_filtered) <- paste0("GP", seq_len(ncol(F_pm_filtered)))
# Importance vector D: the max abs of each column of F_pm_filtered
D <- apply(F_pm_filtered, 2, function(col) max(abs(col), na.rm = TRUE))
# Make D into a diagonal matrix for later use
D_matrix <- diag(D)
seurat_meta <- readRDS(paste0(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
df_sig = read.csv(
  "data/GSEA_signatures_select_toplot.csv",
  header = T,
  sep = ","
)


#####################################################
#####################################################
#####################################################
### Additional Functions needed
#####################################################
#####################################################
top_quadrant_labels <- function(df, x, y, label, n = 10) {
  x <- rlang::ensym(x)
  y <- rlang::ensym(y)
  label <- rlang::ensym(label)

  df %>%
    mutate(
      quad = case_when(
        !!x >= 0 & !!y >= 0 ~ "++",
        !!x < 0 & !!y < 0 ~ "--",
        !!x >= 0 & !!y < 0 ~ "+-",
        TRUE ~ "-+"
      ),
      score = case_when(
        quad == "++" ~ (!!x) + (!!y),
        quad == "--" ~ -(!!x) - (!!y),
        quad == "+-" ~ (!!x) - (!!y),
        quad == "-+" ~ -(!!x) + (!!y)
      ),
      lab_col = recode(
        quad,
        "++" = "darkgreen",
        "--" = "darkred",
        "+-" = "darkblue",
        "-+" = "purple"
      )
    ) %>%
    group_by(quad) %>%
    slice_max(score, n = n, with_ties = FALSE) %>%
    ungroup()
}
FlashierDGE_corrected <- function(F1, L1, group1, group2, title_plot = "") {
  loadings_group1 = colMeans(L1[group1, ])
  loadings_group2 = colMeans(L1[group2, ])
  loadings_groups = colMeans(L1[c(group1, group2), ])
  mean_genes_group1 = F1 %*% loadings_group1
  mean_genes_group2 = F1 %*% loadings_group2
  mean_genes = F1 %*% colMeans(L1[c(group1, group2), ])
  mean_change_loadings = loadings_group1 - loadings_group2
  fc_genes = F1 %*% mean_change_loadings %>% as.data.frame()

  vplot = data.frame(
    SYMBOL = names(mean_change_loadings),
    mean_change_loadings = mean_change_loadings,
    AveExpr = loadings_groups
  )
  max_mean_change = ceiling(max(abs(vplot$mean_change_loadings)))
  top_genes = vplot %>%
    dplyr::arrange(dplyr::desc(abs(mean_change_loadings))) %>%
    utils::head(50)
  p1 = ggplot2::ggplot(data = vplot) +
    ggplot2::geom_point(
      ggplot2::aes(x = mean_change_loadings, y = AveExpr),
      colour = "black",
      alpha = I(1),
      size = I(1)
    ) +
    ggplot2::xlim(-max_mean_change, max_mean_change) +
    ggrepel::geom_text_repel(
      data = top_genes,
      ggplot2::aes(x = mean_change_loadings, y = AveExpr, label = SYMBOL),
      size = 3,
      color = "red",
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = "grey50",
      max.overlaps = 20
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Difference in Mean Loading",
      y = "Average Loading",
      title = title_plot
    )
  vplot_genes = data.frame(
    SYMBOL = rownames(fc_genes),
    log2FC = fc_genes[,
      1
    ] /
      log(2),
    AveExpr = mean_genes[, 1]
  )
  max_fc_genes = ceiling(max(abs(vplot_genes$log2FC)))
  top_genes_genes <- vplot_genes %>%
    dplyr::arrange(dplyr::desc(abs(log2FC))) %>%
    utils::head(50)
  p2 <- ggplot2::ggplot(data = vplot_genes) +
    scattermore::geom_scattermore(
      ggplot2::aes(x = log2FC, y = AveExpr),
      colour = "black",
      alpha = I(1),
      size = I(1),
      pixels = c(512, 512)
    ) +
    ggplot2::xlim(-max_fc_genes, max_fc_genes) +
    ggrepel::geom_text_repel(
      data = top_genes_genes,
      ggplot2::aes(x = log2FC, y = AveExpr, label = SYMBOL),
      size = 3,
      color = "red",
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = "grey50",
      max.overlaps = 20
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Fold Change (log2)",
      y = "Average Expression",
      title = title_plot
    )
  list(p1 = p1, p2 = p2, diff_factors = vplot, diff_genes = vplot_genes)
}
plot_target_gps <- function(
  df,
  x_var,
  y_var,
  label_var,
  target_gps,
  highlight_color = "darkorange",
  x_limits = c(-0.1, 0.2),
  y_limits = c(-0.2, 0.2),
  background_alpha = 1,
  xlab = "Difference in Mean Loading (CD4 Activated vs Resting)",
  ylab = "Difference in Mean Loading (CD8 Activated vs Resting)",
  title = "Comparison of Specific GP Loadings",
  point_size = 2,
  label_size = 3.5
) {
  # Filter the dataframe for only the specific GPs you care about
  highlight_df <- df %>%
    filter({{ label_var }} %in% target_gps)

  # Check if any GPs were found to prevent ggplot errors
  if (nrow(highlight_df) == 0) {
    warning("None of the target_gps were found in the dataframe.")
  }

  # Build the scatter plot
  p <- ggplot(df, aes(x = {{ x_var }}, y = {{ y_var }})) +

    # 1. Plot background points (unselected GPs) in black
    geom_point(color = "black", alpha = background_alpha) +

    # 2. Add reference lines matching your original plot
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +

    # 3. Plot the targeted GPs on top, mapping color inside aes()
    geom_point(
      data = highlight_df,
      aes(color = {{ label_var }}),
      size = point_size
    ) +

    # 4. Label only the targeted GPs, mapping color inside aes()
    ggrepel::geom_text_repel(
      data = highlight_df,
      aes(label = {{ label_var }}, color = {{ label_var }}),
      max.overlaps = Inf,
      size = label_size,
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = "grey50",
      show.legend = FALSE
    ) +

    # 5. Force ggplot to use your named color vector!
    scale_color_manual(values = highlight_color, guide = "none") +

    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal()

  return(p)
}
plot_target_gps_with_quadrants <- function(
  df,
  x_var,
  y_var,
  label_var,
  target_gps,
  highlight_color_vec,
  top_n = 3,
  x_limits = c(-0.1, 0.2),
  y_limits = c(-0.2, 0.2),
  background_alpha = 0.5,
  xlab = "Difference in Mean Loading (CD4 Activated vs Resting)",
  ylab = "Difference in Mean Loading (CD8 Activated vs Resting)",
  title = "Manual Highlights vs. Top Automated Discoveries"
) {
  # 1. Identify Top Quadrant GPs EXCLUDING the manual list
  # We use the raw column names for the quadrant calculation
  df_for_quads <- df %>% filter(!({{ label_var }} %in% target_gps))

  quad_hits <- top_quadrant_labels(
    df_for_quads,
    !!enquo(x_var),
    !!enquo(y_var),
    !!enquo(label_var),
    n = top_n
  )

  # Ensure we only take quadrants that actually had GPs to offer
  quad_gps <- quad_hits %>% pull({{ label_var }})

  # 2. Create the combined highlight dataframe
  # Grouping ensures we can apply different colors to manual vs auto hits
  highlight_df <- df %>%
    filter({{ label_var }} %in% c(target_gps, quad_gps)) %>%
    mutate(
      group = case_when(
        {{ label_var }} %in% target_gps ~ as.character({{ label_var }}),
        TRUE ~ "Top Automated Discovery"
      )
    )

  # 3. Update the color vector
  # Preserves your manual colors and adds 'Black' for the new discoveries
  final_colors <- c(highlight_color_vec, "Top Automated Discovery" = "black")

  # 4. Build Plot
  p <- ggplot(df, aes(x = {{ x_var }}, y = {{ y_var }})) +
    # Background
    geom_point(color = "grey85", alpha = background_alpha, size = 1) +

    # Reference lines
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey70"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +

    # Highlight points
    geom_point(data = highlight_df, aes(color = group), size = 2.5) +

    # Labels
    ggrepel::geom_text_repel(
      data = highlight_df,
      aes(label = {{ label_var }}, color = group),
      max.overlaps = Inf,
      size = 3.2,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50",
      fontface = "bold",
      show.legend = FALSE
    ) +

    scale_color_manual(values = final_colors) +

    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    labs(
      x = xlab,
      y = ylab,
      title = title,
      subtitle = paste0(
        "Black labels indicate top ",
        top_n,
        " discoveries per quadrant (excluding manual hits)"
      )
    ) +
    theme_minimal() +
    theme(legend.position = "none") # Keeping it clean since labels are self-explanatory

  return(p)
}


#####################################################
#####################################################
#####################################################
### Scatter plot comparing differential loadings between CD4 and CD8 activation
#####################################################
#####################################################
# DGE of CD4 activated vs resting
CD4_cells <- seurat_meta_filtered$cellID[
  seurat_meta_filtered$annotation_level1 == "CD4"
]
CD4_resting_cells <- CD4_cells[
  seurat_meta_filtered$annotation_level2_group[match(
    CD4_cells,
    seurat_meta_filtered$cellID
  )] ==
    "resting"
]
CD4_activated_cells <- CD4_cells[
  seurat_meta_filtered$annotation_level2_group[match(
    CD4_cells,
    seurat_meta_filtered$cellID
  )] ==
    "activated"
]
CD4_DGE <- FlashierDGE_corrected(
  F1 = F_pm_filtered,
  L1 = L_pm_filtered,
  group1 = CD4_activated_cells,
  group2 = CD4_resting_cells
)
# DGE of CD8 activated vs resting
CD8_cells <- seurat_meta_filtered$cellID[
  seurat_meta_filtered$annotation_level1 == "CD8"
]
CD8_resting_cells <- CD8_cells[
  seurat_meta_filtered$annotation_level2_group[match(
    CD8_cells,
    seurat_meta_filtered$cellID
  )] ==
    "resting"
]
CD8_activated_cells <- CD8_cells[
  seurat_meta_filtered$annotation_level2_group[match(
    CD8_cells,
    seurat_meta_filtered$cellID
  )] ==
    "activated"
]
CD8_DGE <- FlashierDGE_corrected(
  F1 = F_pm_filtered,
  L1 = L_pm_filtered,
  group1 = CD8_activated_cells,
  group2 = CD8_resting_cells
)
diff_factors_CD4 <- CD4_DGE$diff_factors %>%
  rename(mean_change_loadings_CD4 = mean_change_loadings, AveExpr_CD4 = AveExpr)
diff_factors_CD8 <- CD8_DGE$diff_factors %>%
  rename(mean_change_loadings_CD8 = mean_change_loadings, AveExpr_CD8 = AveExpr)
diff_factors_merged <- diff_factors_CD4 %>%
  inner_join(diff_factors_CD8, by = "SYMBOL")
# Define the specific GPs you care about
GPs_of_interest <- paste0(
  "GP",
  c(
    25,
    26,
    10,
    12,
    58,
    171,
    9,
    79,
    35,
    177,
    161,
    152,
    162,
    36,
    56,
    181,
    32,
    80,
    49,
    41,
    57,
    176,
    11,
    13,
    159
  )
)
highlight_colors <- c(
  "GP56" = "blue",
  "GP162" = "blue",
  "GP36" = "blue",
  "GP152" = "blue",
  "GP161" = "blue",
  "GP177" = "blue",
  "GP79" = "blue",
  "GP12" = "blue",
  "GP13" = "blue",
  "GP159" = "blue",
  "GP10" = "darkorange2",
  "GP58" = "darkorange2",
  "GP181" = "darkorange2",
  "GP176" = "darkorange2",
  "GP25" = "darkred",
  "GP26" = "darkred",
  "GP35" = "darkred",
  "GP32" = "darkred",
  "GP80" = "darkred",
  "GP57" = "darkred",
  "GP9" = "darkgreen",
  "GP171" = "darkgreen",
  "GP49" = "darkgreen",
  "GP41" = "darkgreen",
  "GP11" = "darkgreen"
)

#####################################################
#####################################################
#####################################################
### Standardized mean difference (activated vs resting) per GP.
### d_k = (mean_act - mean_rest) / sd_pooled
### sd_pooled = SD of GP k's loading over ALL activated+resting CD4 and CD8
### cells (one per-GP denominator shared by both lineages). Unlike a Welch z,
### the denominator is a loading SD, not a standard error, so d does NOT depend
### on the (unequal) CD4 vs CD8 cell counts and is directly comparable across
### lineages. Units are "loading SDs" (Cohen's-d style). Because the denominator
### is shared, d_CD8 / d_CD4 still equals the mean-change Ratio_CD8_CD4.
#####################################################
#####################################################
std_mean_diff <- function(mat, idx_a, idx_b, sd_ref) {
  d <- (colMeans(mat[idx_a, , drop = FALSE]) -
    colMeans(mat[idx_b, , drop = FALSE])) /
    sd_ref
  d[!is.finite(d)] <- 0
  d
}

# Shared per-GP denominator: SD over all activated+resting CD4 and CD8 cells.
sd_pooled_act_rest <- apply(
  L_pm_filtered[
    c(
      CD4_activated_cells,
      CD4_resting_cells,
      CD8_activated_cells,
      CD8_resting_cells
    ),
    ,
    drop = FALSE
  ],
  2,
  sd
)
d_CD4 <- std_mean_diff(
  L_pm_filtered,
  CD4_activated_cells,
  CD4_resting_cells,
  sd_pooled_act_rest
)
d_CD8 <- std_mean_diff(
  L_pm_filtered,
  CD8_activated_cells,
  CD8_resting_cells,
  sd_pooled_act_rest
)
d_factors_merged <- data.frame(
  SYMBOL = colnames(L_pm_filtered),
  d_CD4 = d_CD4,
  d_CD8 = d_CD8
)


#####################################################
#####################################################
#####################################################
### Scatter plot comparing differential loadings between CD4 and CD8 activation
#####################################################
#####################################################
# Highlight every GP with |d| > d_thr in CD4 or CD8 (effect-size criterion).
# Color by abs(Ratio_CD8_CD4) = |mean_change_CD8 / mean_change_CD4|:
#   > ratio_cutoff   -> CD8 dominates  (darkorange2)
#   < 1/ratio_cutoff -> CD4 dominates  (blue)
#   otherwise both contribute, color by sign of the mean changes:
#     both up   -> darkred
#     both down -> darkgreen
#     mixed     -> grey50  (rare — opposite-sign effects)
# Unified magnitude gate: drives the p_quadrant auto-highlight filter, the
# four count summaries below, and the manual_curated_df colouring rule.
d_thr <- 0.15
ratio_cutoff <- 3

highlight_df_sel <- diff_factors_merged %>%
  inner_join(
    d_factors_merged %>% select(SYMBOL, d_CD4, d_CD8),
    by = "SYMBOL"
  ) %>%
  # Union: GPs passing the auto d-threshold OR any manually-curated GP in
  # GPs_of_interest (so the manual list is always shown even if d < thr).
  filter(
    abs(d_CD4) > d_thr |
      abs(d_CD8) > d_thr |
      SYMBOL %in% GPs_of_interest
  ) %>%
  mutate(
    abs_ratio = abs(mean_change_loadings_CD8 / mean_change_loadings_CD4),
    color_group = case_when(
      abs_ratio > ratio_cutoff ~ "darkorange2",
      abs_ratio < 1 / ratio_cutoff ~ "blue",
      mean_change_loadings_CD4 > 0 & mean_change_loadings_CD8 > 0 ~ "darkred",
      mean_change_loadings_CD4 < 0 & mean_change_loadings_CD8 < 0 ~ "darkgreen",
      TRUE ~ "grey50"
    )
  )
target_gps_diff_mean <- highlight_df_sel$SYMBOL
highlight_colors_diff_mean <- setNames(
  highlight_df_sel$color_group,
  highlight_df_sel$SYMBOL
)

p_quadrant <- plot_target_gps(
  df = diff_factors_merged,
  x_var = mean_change_loadings_CD4,
  y_var = mean_change_loadings_CD8,
  label_var = SYMBOL,
  target_gps = target_gps_diff_mean,
  background_alpha = 1,
  highlight_color = highlight_colors_diff_mean,
  title = paste0(
    "GPs with |d| > ",
    d_thr,
    " — colored by abs(Ratio_CD8_CD4) (cutoff ",
    ratio_cutoff,
    ")"
  )
)
print(p_quadrant)
ggsave(
  filename = paste0(figure_path, "DGE_diff_mean.pdf"),
  plot = p_quadrant,
  width = 8,
  height = 7
)

# Per-GP summary table — mean loading differences, average loadings,
# standardized mean differences (d), and the CD8/CD4 mean-difference ratio.
GP_activation_summary <- diff_factors_merged %>%
  inner_join(
    d_factors_merged %>% select(SYMBOL, d_CD4, d_CD8),
    by = "SYMBOL"
  ) %>%
  mutate(
    Ratio_CD8_CD4 = mean_change_loadings_CD8 / mean_change_loadings_CD4
  ) %>%
  select(
    GP = SYMBOL,
    mean_change_loadings_CD4,
    mean_change_loadings_CD8,
    AveExpr_CD4,
    AveExpr_CD8,
    d_CD4,
    d_CD8,
    Ratio_CD8_CD4
  )

GP_activation_summary$Category <- case_when(
  abs(GP_activation_summary$Ratio_CD8_CD4) > ratio_cutoff &
    abs(GP_activation_summary$d_CD8) > d_thr ~ "CD8-dominant",
  abs(GP_activation_summary$Ratio_CD8_CD4) < 1 / ratio_cutoff &
    abs(GP_activation_summary$d_CD4) > d_thr ~ "CD4-dominant",
  abs(GP_activation_summary$Ratio_CD8_CD4) > 1 / ratio_cutoff &
    abs(GP_activation_summary$Ratio_CD8_CD4) < ratio_cutoff &
    GP_activation_summary$d_CD4 > d_thr &
    GP_activation_summary$d_CD8 > d_thr ~ "Both up",
  abs(GP_activation_summary$Ratio_CD8_CD4) > 1 / ratio_cutoff &
    abs(GP_activation_summary$Ratio_CD8_CD4) < ratio_cutoff &
    GP_activation_summary$d_CD4 < -d_thr &
    GP_activation_summary$d_CD8 < -d_thr ~ "Both down",
  TRUE ~ "Other"
)


write.csv(
  GP_activation_summary,
  file = paste0(figure_path, "GP_activation_summary.csv"),
  row.names = FALSE
)


# Number of GPs strongly CD8-dominant (|Ratio_CD8_CD4| > ratio_cutoff) with
# strong activation signal in CD8 (|d_CD8| > d_thr).
sum(
  abs(GP_activation_summary$Ratio_CD8_CD4) > ratio_cutoff &
    abs(GP_activation_summary$d_CD8) > d_thr
)

# Number of GPs strongly CD4-dominant (|Ratio_CD8_CD4| < 1/ratio_cutoff) with
# strong activation signal in CD4 (|d_CD4| > d_thr).
sum(
  abs(GP_activation_summary$Ratio_CD8_CD4) < 1 / ratio_cutoff &
    abs(GP_activation_summary$d_CD4) > d_thr
)

# Number of GPs strongly up in both CD4 and CD8 (d_CD4 > d_thr &
# d_CD8 > d_thr), |Ratio_CD8_CD4| within (1/ratio_cutoff, ratio_cutoff).
sum(
  abs(GP_activation_summary$Ratio_CD8_CD4) > 1 / ratio_cutoff &
    abs(GP_activation_summary$Ratio_CD8_CD4) < ratio_cutoff &
    GP_activation_summary$d_CD4 > d_thr &
    GP_activation_summary$d_CD8 > d_thr
)

# Number of GPs strongly down in both CD4 and CD8 (d_CD4 < -d_thr &
# d_CD8 < -d_thr), |Ratio_CD8_CD4| within (1/ratio_cutoff, ratio_cutoff).
sum(
  abs(GP_activation_summary$Ratio_CD8_CD4) > 1 / ratio_cutoff &
    abs(GP_activation_summary$Ratio_CD8_CD4) < ratio_cutoff &
    GP_activation_summary$d_CD4 < -d_thr &
    GP_activation_summary$d_CD8 < -d_thr
)

# Colour every GP using the same four-category rule as the counts above
# (ratio + sign + magnitude gate via d_thr). Curated GPs (GPs_of_interest)
# override with their fixed manual highlight_colors; non-curated GPs are
# classified automatically:
#   |Ratio_CD8_CD4| > ratio_cutoff,   |d_CD8| > d_thr      -> darkorange2 (CD8-dominant)
#   |Ratio_CD8_CD4| < 1/ratio_cutoff, |d_CD4| > d_thr      -> blue        (CD4-dominant)
#   |Ratio_CD8_CD4| in (1/ratio_cutoff, ratio_cutoff), both d > d_thr  -> darkred   (both up)
#   |Ratio_CD8_CD4| in (1/ratio_cutoff, ratio_cutoff), both d < -d_thr -> darkgreen (both down)
#   anything else (fails magnitude gate / opposite signs)        -> black
# Only the curated GPs are labelled, to keep the plot readable.
manual_curated_df <- GP_activation_summary %>%
  mutate(
    auto_color = case_when(
      abs(Ratio_CD8_CD4) > ratio_cutoff &
        abs(d_CD8) > d_thr ~ "darkorange2",
      abs(Ratio_CD8_CD4) < 1 / ratio_cutoff &
        abs(d_CD4) > d_thr ~ "blue",
      abs(Ratio_CD8_CD4) > 1 / ratio_cutoff &
        abs(Ratio_CD8_CD4) < ratio_cutoff &
        d_CD4 > d_thr &
        d_CD8 > d_thr ~ "darkred",
      abs(Ratio_CD8_CD4) > 1 / ratio_cutoff &
        abs(Ratio_CD8_CD4) < ratio_cutoff &
        d_CD4 < -d_thr &
        d_CD8 < -d_thr ~ "darkgreen",
      TRUE ~ "black"
    ),
    point_color = ifelse(
      GP %in% GPs_of_interest,
      highlight_colors[GP],
      auto_color
    )
  )

p_quadrant_manual <- ggplot(manual_curated_df, aes(x = d_CD4, y = d_CD8)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
  geom_point(aes(color = point_color), size = 2) +
  ggrepel::geom_text_repel(
    data = filter(manual_curated_df, GP %in% GPs_of_interest),
    aes(label = GP, color = point_color),
    max.overlaps = Inf,
    size = 3.5,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = "grey50"
  ) +
  scale_color_identity() +
  # Signed (pseudo-)log axes: d is signed, so a plain log drops negatives/zeros.
  # pseudo_log is symmetric about 0, spreads out the crowded near-zero points
  # and compresses the tails. sigma sets the half-width of the ~linear zone
  # around 0 (smaller sigma -> more expansion near 0). y=x and the 0 axes stay
  # straight since both axes share the transform.
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.15),
    breaks = c(-1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1)
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.15),
    breaks = c(-1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1)
  ) +
  coord_equal(xlim = c(-1.6, 1.6), ylim = c(-1.6, 1.6)) +
  labs(
    x = "Standardized Mean Difference (d) for CD4 Activated vs Resting",
    y = "Standardized Mean Difference (d) for CD8 Activated vs Resting",
    title = "GPs colored by semantic category (curated set labeled)"
  ) +
  theme_minimal()
print(p_quadrant_manual)
ggsave(
  filename = paste0(figure_path, "DGE_diff_mean_curated.pdf"),
  plot = p_quadrant_manual,
  width = 8,
  height = 7
)

# Proportion of activated cells with highly active GP26
prop_active_GP26_activated <- mean(
  L_pm_filtered[c(CD4_activated_cells, CD8_activated_cells), "GP26"] > 0.1
)
prop_active_GP26_activated

# Proportion of activated cells with highly active GP80
prop_active_GP80_activated <- mean(
  L_pm_filtered[c(CD4_activated_cells, CD8_activated_cells), "GP80"] > 0.1
)
prop_active_GP80_activated

# Proportion of activated cells with both highly active GP26 and GP80
prop_active_both_GP26_GP80_activated <- mean(
  L_pm_filtered[c(CD4_activated_cells, CD8_activated_cells), "GP26"] > 0.1 &
    L_pm_filtered[c(CD4_activated_cells, CD8_activated_cells), "GP80"] > 0.1
)
prop_active_both_GP26_GP80_activated

# Proportion of activated cells with either highly active GP26 or GP80
prop_active_either_GP26_GP80_activated <- mean(
  L_pm_filtered[c(CD4_activated_cells, CD8_activated_cells), "GP26"] > 0.1 |
    L_pm_filtered[c(CD4_activated_cells, CD8_activated_cells), "GP80"] > 0.1
)
prop_active_either_GP26_GP80_activated


#####################################################
#####################################################
#####################################################
### Heatmap of GP log2 fold changes for CD4 vs CD8 activated vs resting comparisons,
### with nested annotation for cell type and activation status
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
# 1. Define target GPs and subset the matrix
L_pm_filtered_subset <- L_pm_filtered[
  c(CD4_cells, CD8_cells),
  GPs_of_interest,
  drop = FALSE
]
L <- L_pm_filtered_subset # cells x GP matrix
cell_ids <- rownames(L)
# 2. Match metadata
cell_type <- seurat_meta_filtered$annotation_level2[
  match(cell_ids, seurat_meta_filtered$cellID)
]
# 3. Downsample cells based on budget
set.seed(123)
total_budget <- 8000
small_keep_all <- 150 # clusters <= this keep all cells
min_per_cluster <- 150 # clusters larger than small_keep_all get at least this many
max_per_cluster <- 600 # optional cap per cluster
cells_by_ct <- split(cell_ids, cell_type)
sizes <- sapply(cells_by_ct, length)
alloc <- ifelse(sizes <= small_keep_all, sizes, pmin(min_per_cluster, sizes))
remaining <- total_budget - sum(alloc)
if (remaining > 0) {
  room <- pmax(pmin(sizes, max_per_cluster) - alloc, 0)
  if (sum(room) > 0) {
    extra <- floor(remaining * room / sum(room))
    alloc <- alloc + extra

    # Fill greedily if there is leftover budget due to flooring
    leftover <- total_budget - sum(alloc)
    if (leftover > 0) {
      idx <- order(room, decreasing = TRUE)
      for (j in idx) {
        if (leftover <= 0) {
          break
        }
        addable <- pmin(room[j] - extra[j], leftover)
        if (addable > 0) {
          alloc[j] <- alloc[j] + addable
          leftover <- leftover - addable
        }
      }
    }
  }
}
sampled_cells <- unlist(
  mapply(
    function(v, m) {
      if (length(v) <= m) v else sample(v, m)
    },
    cells_by_ct,
    pmin(alloc, sizes),
    SIMPLIFY = FALSE
  ),
  use.names = FALSE
)
# 4. Clean up metadata for the sampled cells
cell_group_s <- seurat_meta_filtered$annotation_level2_group[match(
  sampled_cells,
  seurat_meta_filtered$cellID
)]
cell_level2_s <- seurat_meta_filtered$annotation_level2[match(
  sampled_cells,
  seurat_meta_filtered$cellID
)]
cell_group_s <- trimws(tolower(as.character(cell_group_s)))
cell_level2_s <- trimws(as.character(cell_level2_s))
cell_level1_s <- sapply(strsplit(cell_level2_s, "[_.]"), `[`, 1)
is_w_cell <- grepl("\\.w", cell_level2_s)
valid_idx <- which(
  !is.na(cell_level2_s) &
    !tolower(cell_level2_s) %in% c("", "na", "nan") &
    cell_group_s %in% c("resting", "activated") &
    !is_w_cell
)
final_cells <- sampled_cells[valid_idx]
final_group <- cell_group_s[valid_idx]
final_level1 <- cell_level1_s[valid_idx]
final_level2 <- cell_level2_s[valid_idx]
# Order columns by group, then level1, then level2
col_order <- order(final_group, final_level1, final_level2)
final_cells <- final_cells[col_order]
final_group <- final_group[col_order]
final_level1 <- final_level1[col_order]
final_level2 <- final_level2[col_order]
# 5. Heatmap constants and GP ordering
pc <- 1e-10
cap <- 2
# Define the exact order you want the GPs to appear in the heatmap
ordered_GPs <- c(
  "GP56",
  "GP162",
  "GP36",
  "GP152",
  "GP161",
  "GP177",
  "GP79",
  "GP12",
  "GP13",
  "GP159", # "Blue" group
  "GP10",
  "GP58",
  "GP181",
  "GP176", # "Orange" group
  "GP25",
  "GP26",
  "GP35",
  "GP32",
  "GP80",
  "GP57", # "Darkred" group (both up)
  "GP9",
  "GP171",
  "GP49",
  "GP41",
  "GP11" # "Green" group
)
# Calculate where to put physical gaps (white lines) between the groups
group_counts <- sapply(
  list(
    c(
      "GP56",
      "GP162",
      "GP36",
      "GP152",
      "GP161",
      "GP177",
      "GP79",
      "GP12",
      "GP13",
      "GP159"
    ),
    c("GP10", "GP58", "GP181", "GP176"),
    c("GP25", "GP26", "GP35", "GP32", "GP80", "GP57"),
    c("GP9", "GP171", "GP49", "GP41", "GP11")
  ),
  function(gp_group) sum(ordered_GPs %in% gp_group)
)
group_counts <- group_counts[group_counts > 0]
gaps_row <- cumsum(group_counts)[-length(group_counts)]


#####################################################
### Activated-only heatmap
#####################################################
act_idx <- which(final_group == "activated")
act_cells <- final_cells[act_idx]
act_level1 <- final_level1[act_idx]
act_level2 <- final_level2[act_idx]

# Re-order by level1 then level2 (no Group split needed)
act_order <- order(act_level1, act_level2)
act_cells <- act_cells[act_order]
act_level1 <- act_level1[act_order]
act_level2 <- act_level2[act_order]

L_sub_act <- L[act_cells, , drop = FALSE]
M_raw_act <- t(L_sub_act)

# Per-lineage resting baseline: each activated cell's log2FC is computed
# against the mean loading in resting cells of its own lineage
# (CD4 cell -> mean across CD4 resting; CD8 cell -> mean across CD8 resting).
mu_resting_CD4 <- colMeans(
  L_pm_filtered[CD4_resting_cells, GPs_of_interest, drop = FALSE]
)
mu_resting_CD8 <- colMeans(
  L_pm_filtered[CD8_resting_cells, GPs_of_interest, drop = FALSE]
)
baseline_mat <- matrix(NA_real_, nrow = nrow(M_raw_act), ncol = ncol(M_raw_act))
rownames(baseline_mat) <- rownames(M_raw_act)
baseline_mat[, act_level1 == "CD4"] <- mu_resting_CD4[rownames(M_raw_act)]
baseline_mat[, act_level1 == "CD8"] <- mu_resting_CD8[rownames(M_raw_act)]
M_fc_act <- log2((M_raw_act + pc) / (baseline_mat + pc))
M_fc_act_cap <- pmax(pmin(M_fc_act, cap), -cap)
M_fc_act_cap <- M_fc_act_cap[ordered_GPs, , drop = FALSE]

ann_col_act <- data.frame(
  Level2 = factor(act_level2),
  Level1 = factor(act_level1)
)
rownames(ann_col_act) <- colnames(M_fc_act_cap)

present_level1_act <- levels(ann_col_act$Level1)
present_level2_act <- levels(ann_col_act$Level2)
level1_cols_act <- ZemmourLib::immgent_colors$level1
level1_cols_act <- level1_cols_act[
  names(level1_cols_act) %in% present_level1_act
]
level2_cols_act <- ZemmourLib::immgent_colors$level2
level2_cols_act <- level2_cols_act[
  names(level2_cols_act) %in% present_level2_act
]
missing_l2_act <- setdiff(present_level2_act, names(level2_cols_act))
if (length(missing_l2_act) > 0) {
  level2_cols_act <- c(
    level2_cols_act,
    setNames(rep("grey80", length(missing_l2_act)), missing_l2_act)
  )
}
ann_colors_act <- list(Level1 = level1_cols_act, Level2 = level2_cols_act)

# Column gaps between Level1 groups
rle_l1_act <- rle(as.character(ann_col_act$Level1))
gaps_col_act <- cumsum(rle_l1_act$lengths)
gaps_col_act <- gaps_col_act[-length(gaps_col_act)]

pheatmap(
  M_fc_act_cap,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gaps_row,
  gaps_col = gaps_col_act,
  color = colorRampPalette(c("#7A0177", "black", "#FFD700"))(101),
  breaks = seq(-cap, cap, length.out = 102),
  show_colnames = FALSE,
  fontsize_row = 7,
  fontface_row = "bold",
  border_color = NA,
  annotation_col = ann_col_act,
  annotation_colors = ann_colors_act,
  annotation_names_col = TRUE,
  useRaster = TRUE,
  main = "GP Log2FC vs per-lineage RESTING baseline (CD4 cells vs CD4 resting; CD8 vs CD8 resting)",
  filename = paste0(figure_path, "DGE_heatmap_level12_meanref.pdf"),
  width = 12,
  height = 4
)


#####################################################
### Structure plot of a curated subset of highlighted GPs across the 4
### lineage-x-state groups (CD4_resting, CD4_activated, CD8_resting,
### CD8_activated).
###
### Cells: CD4 + CD8 with annotation_level2_group in {resting, activated}.
### Factors: a small subset of GPs (see `structure_GPs`) — keeps the plot
### readable.
#####################################################
# Curated subset for the structure plot (small enough to read easily).
structure_GPs <- paste0(
  "GP",
  c(25, 26, 80, 11, 171, 9, 12, 56, 159, 58, 10)
)

set.seed(1234)
struct_idx <- seurat_meta_filtered$annotation_level1 %in%
  c("CD4", "CD8") &
  seurat_meta_filtered$annotation_level2_group %in% c("resting", "activated")
struct_meta <- seurat_meta_filtered[struct_idx, ]
struct_group <- factor(
  paste0(
    struct_meta$annotation_level1,
    "_",
    struct_meta$annotation_level2_group
  ),
  levels = c("CD4_resting", "CD4_activated", "CD8_resting", "CD8_activated")
)
struct_mat <- L_pm_filtered[struct_meta$cellID, structure_GPs, drop = FALSE]

# Drop the very smallest groups if any (need >= 50 cells for a sane bar)
keep_groups <- names(table(struct_group))[table(struct_group) >= 50]
keep_cells_idx <- which(as.character(struct_group) %in% keep_groups)
struct_mat <- struct_mat[keep_cells_idx, , drop = FALSE]
struct_group <- factor(
  as.character(struct_group)[keep_cells_idx],
  levels = intersect(
    c("CD4_resting", "CD4_activated", "CD8_resting", "CD8_activated"),
    keep_groups
  )
)

# Each GP gets a unique shade within its color family (so different GPs in
# the same group are still distinguishable in the stacked bars). Endpoints
# are dark -> mid of the family; order within family follows the order in
# `highlight_colors` (which is grouped by color).
gp_family_endpoints <- list(
  blue = c("#08306b", "#6baed6"),
  darkorange2 = c("#7f2704", "#fdae6b"),
  darkred = c("#67000d", "#fb6a4a"),
  darkgreen = c("#00441b", "#74c476")
)
# Build unique per-GP colors within the subset (so each family's GPs are
# visually distinguishable; ramp length is based on how many subset GPs
# belong to that family).
subset_colors <- highlight_colors[structure_GPs]
gp_unique_colors <- character(length(structure_GPs))
names(gp_unique_colors) <- structure_GPs
for (fam in names(gp_family_endpoints)) {
  gps_fam <- names(subset_colors)[subset_colors == fam]
  if (!length(gps_fam)) {
    next
  }
  ramp <- colorRampPalette(gp_family_endpoints[[fam]])(length(gps_fam))
  gp_unique_colors[gps_fam] <- ramp
}

# Reorder struct_mat columns by canonical family order
# (blue -> orange -> red -> green), preserving the user's order within
# each family.
fam_order <- c("blue", "darkorange2", "darkred", "darkgreen")
struct_mat <- struct_mat[,
  structure_GPs[order(match(subset_colors, fam_order))],
  drop = FALSE
]

p_structure <- fastTopics::structure_plot(
  struct_mat,
  grouping = struct_group,
  colors = gp_unique_colors[colnames(struct_mat)],
  gap = 40,
  n = 8000
) +
  labs(y = "GP loading", fill = "GP", color = "GP") +
  guides(
    fill = guide_legend(nrow = 3),
    color = guide_legend(nrow = 3)
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 9)
  )
ggsave(
  filename = paste0(figure_path, "StructurePlot_highlightedGPs.pdf"),
  plot = p_structure,
  width = 12,
  height = 6,
  dpi = 300
)


#####################################################
#####################################################
#####################################################
### Network plot of top GP-gene associations colored by manual GP classification (Blue, Orange, Green groups)
#####################################################
#####################################################
set.seed(42)
# 1. Normalize and subset F matrix
F_pm_filtered_norm <- scale_cols(
  F_pm_filtered,
  1 / apply(abs(F_pm_filtered), 2, max)
)
colnames(F_pm_filtered_norm) <- paste0("GP", seq_len(ncol(F_pm_filtered_norm)))
F_pm_filtered_norm_subset <- F_pm_filtered_norm[, GPs_of_interest, drop = FALSE]
# 2. Extract top positive and negative genes
top_5_pos <- apply(F_pm_filtered_norm_subset, 2, function(x) {
  idx <- order(abs(x), decreasing = TRUE)[1:5]
  idx <- idx[x[idx] > 0]
  names(x)[idx]
})
top_5_neg <- apply(F_pm_filtered_norm_subset, 2, function(x) {
  idx <- order(abs(x), decreasing = TRUE)[1:5]
  idx <- idx[x[idx] < 0]
  names(x)[idx]
})
names(top_5_pos) <- GPs_of_interest
names(top_5_neg) <- GPs_of_interest
# 3. Create Edge List
pos_edges <- stack(top_5_pos) %>%
  rename(Gene = values, GP = ind) %>%
  mutate(Type = "Positive", Color = "red")
neg_edges <- stack(top_5_neg) %>%
  rename(Gene = values, GP = ind) %>%
  mutate(Type = "Negative", Color = "blue")
all_edges <- bind_rows(pos_edges, neg_edges) %>%
  filter(Gene != "" & !is.na(Gene))
all_edges_sorted <- all_edges %>%
  arrange(Type, GP, Gene)
# A. Create a lookup table for the nodes (uses the global highlight_colors)
gp_group_df <- data.frame(
  name = names(highlight_colors),
  ManualGroup = case_when(
    highlight_colors == "blue" ~ "CD4 only",
    highlight_colors == "darkorange2" ~ "CD8 only",
    highlight_colors == "darkgreen" ~ "both down",
    highlight_colors == "darkred" ~ "both up",
  )
)
# C. Define the exact color hex/names for ggplot to use
manual_colors_palette <- c(
  "CD4 only" = "blue",
  "CD8 only" = "darkorange2",
  "both down" = "darkgreen",
  "both up" = "darkred",
  "Gene" = "#666666" # Keep the genes grey
)
# 4. Build the Graph Object
graph <- as_tbl_graph(all_edges_sorted) %>%
  activate(nodes) %>%
  mutate(
    NodeGroup = ifelse(name %in% all_edges$GP, "GP", "Gene"),
    Importance = centrality_degree()
  ) %>%
  left_join(gp_group_df, by = "name") %>%
  mutate(
    # Assign genes the "Gene" group, otherwise use the manual group we just joined
    ColorGroup = ifelse(NodeGroup == "Gene", "Gene", ManualGroup),
    gp_label = ifelse(NodeGroup == "GP", name, "")
  )
# 5. Generate the Network Plot
# Use a stress-majorization layout (graphlayouts::layout_with_stress, via
# ggraph) which spreads GP nodes more evenly than the default "nicely"
# heuristic and is deterministic given a seed.
set.seed(2)
p_network <- ggraph(graph, layout = "stress") +
  # Edges
  geom_edge_link(aes(color = Color), alpha = 0.4, width = 0.6) +
  # Gene Nodes
  geom_node_point(
    aes(filter = (NodeGroup == "Gene"), color = ColorGroup),
    shape = 16,
    size = 2,
    alpha = 0.8
  ) +
  # GP Nodes (smaller squares to reduce visual overlap when GPs are close)
  geom_node_point(
    aes(filter = (NodeGroup == "GP"), color = ColorGroup),
    shape = 15,
    size = 10,
    alpha = 0.7
  ) +
  # GP Labels
  geom_node_text(
    aes(filter = (NodeGroup == "GP"), label = gp_label),
    color = "white",
    fontface = "bold",
    size = 3
  ) +
  # Gene Labels
  geom_node_text(
    aes(filter = (NodeGroup == "Gene"), label = name),
    repel = TRUE,
    size = 2.5,
    color = "black",
    max.overlaps = 20
  ) +

  scale_edge_color_identity() +
  # NEW: Use our manual color palette
  scale_color_manual(
    name = "GP Types",
    values = manual_colors_palette,
    breaks = c("CD4 only", "CD8 only", "both down", "both up") # Hides 'Gene' from legend
  ) +
  theme_void() +
  labs(
    title = "GP-Gene Signature Network",
    subtitle = "Nodes colored by manual GP classification",
    caption = "Red edges: Positive | Blue edges: Negative"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15)))
# 6. Print and Save
print(p_network)
ggsave(
  filename = paste0(figure_path, "Networkplot.pdf"),
  plot = p_network,
  width = 11,
  height = 11
)

### Heatmap of a curated set of genes of interest.
genes_of_interest <- c(
  "Areg",
  "Bhlhe40",
  "Ccl5",
  "Ccr7",
  "Ctla4",
  "Cxcl3",
  "Dusp1",
  "Fos",
  "Foxp3",
  "Gata3",
  "Gzma",
  "Gzmk",
  "Ifng",
  "Ifngr1",
  "Il4",
  "Il10",
  "Il13",
  "Il17a",
  "Il18rap",
  "Il1rl1",
  "Il22",
  "Il7r",
  "Izumo1r",
  "Jun",
  "Klf2",
  "Lag3",
  "Lef1",
  "Nkg7",
  "Nr4a1",
  "Nr4a2",
  "Odc1",
  "Pdcd1",
  "Rgs1",
  "Rora",
  "Stat4",
  "Tbx21",
  "Tcf7",
  "Tigit",
  "Tmem176a",
  "Tmem176b",
  "Tox",
  "Tox2"
)
missing_genes <- setdiff(genes_of_interest, rownames(F_pm_filtered_norm_subset))
if (length(missing_genes) > 0) {
  message("Genes not found in matrix: ", paste(missing_genes, collapse = ", "))
}
important_genes <- intersect(
  genes_of_interest,
  rownames(F_pm_filtered_norm_subset)
)
length(important_genes)
mat_to_plot <- F_pm_filtered_norm_subset[important_genes, , drop = FALSE]

# Group the GP columns by their earlier CD4/CD8 categorization (derived from
# highlight_colors). Columns are ordered by group and separated with gaps.
group_colors <- c(
  "CD4 only" = "blue",
  "CD8 only" = "darkorange2",
  "both up" = "darkred",
  "both down" = "darkgreen"
)
color_to_group <- setNames(names(group_colors), unname(group_colors))
gp_group <- factor(
  color_to_group[unname(highlight_colors[colnames(mat_to_plot)])],
  levels = names(group_colors)
)
names(gp_group) <- colnames(mat_to_plot)

col_order <- order(gp_group)
mat_to_plot <- mat_to_plot[, col_order, drop = FALSE]
gp_group <- gp_group[col_order]

annotation_col <- data.frame(Group = gp_group)
rownames(annotation_col) <- colnames(mat_to_plot)
gaps_col <- head(cumsum(table(gp_group)), -1)

# Non-linear color scale: weights with |value| <= 0.1 stay near white, while
# |value| > 0.1 jumps to a clearly saturated color. Built from custom breaks so
# the central white band is narrow and the colored sides ramp to full intensity.
thr <- 0.1
m <- max(abs(mat_to_plot), na.rm = TRUE)
m <- max(m, thr + 1e-6)
n_side <- 50 # colors on each saturated side
n_core <- 20 # colors in the near-white central band
heat_breaks <- c(
  seq(-m, -thr, length.out = n_side + 1),
  seq(-thr, thr, length.out = n_core + 1)[-1],
  seq(thr, m, length.out = n_side + 1)[-1]
)
heat_colors <- c(
  colorRampPalette(c("blue", "#8a8aff"))(n_side), # -m .. -thr: deep -> medium blue
  colorRampPalette(c("#eef0ff", "white", "#ffeeee"))(n_core), # central near-white band
  colorRampPalette(c("#ff8a8a", "red"))(n_side) # thr .. m: medium -> deep red
)

pheatmap(
  mat_to_plot,
  color = heat_colors,
  breaks = heat_breaks,
  fontsize_row = 8,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = list(Group = group_colors),
  gaps_col = gaps_col,
  # filename = paste0(figure_path, "Heatmap_genes_above_0.1.pdf"),
  width = 8,
  height = max(6, nrow(mat_to_plot) * 0.12)
)





#### Comprehensive heatmap: all genes with |weight| > 0.1, arranged in a block-diagonal
#### layout so that genes whose peak weight falls in each GP group cluster together.
#### Only genes_of_interest are labeled on the rows; a side bar marks their positions.
# Select genes where at least one GP has |score| > 0.1
gene_mask <- apply(abs(F_pm_filtered_norm_subset), 1, max) > 0.1
mat_huge <- F_pm_filtered_norm_subset[gene_mask, , drop = FALSE]

# Column ordering by GP group (same as the curated heatmap)
gp_group_huge <- factor(
  color_to_group[unname(highlight_colors[colnames(mat_huge)])],
  levels = names(group_colors)
)
names(gp_group_huge) <- colnames(mat_huge)
col_order_huge <- order(gp_group_huge)
mat_huge <- mat_huge[, col_order_huge, drop = FALSE]
gp_group_huge <- gp_group_huge[col_order_huge]

# Row ordering: assign each gene to the GP group where its max |weight| falls,
# then sort by (primary_group, sign_of_max desc, column_of_max) so that within
# each block red genes (positive peak) come before blue genes (negative peak).
gene_max_col_idx <- apply(abs(mat_huge), 1, which.max)
gene_primary_group <- gp_group_huge[gene_max_col_idx]
gene_max_sign <- sign(mat_huge[cbind(seq_len(nrow(mat_huge)), gene_max_col_idx)])
row_order_huge <- order(
  gene_max_col_idx,  # one block per GP column (follows the GP-group column order)
  -gene_max_sign     # +1 (red) before -1 (blue) within each GP's genes
)
mat_huge <- mat_huge[row_order_huge, ]
gene_primary_group <- gene_primary_group[row_order_huge]

annotation_col_huge <- data.frame(Group = gp_group_huge)
rownames(annotation_col_huge) <- colnames(mat_huge)
gaps_col_huge <- head(cumsum(table(gp_group_huge)), -1)
gaps_row_huge <- head(cumsum(table(gene_primary_group)), -1)

# Row annotation bar: marks which rows are curated genes of interest
annotation_row_huge <- data.frame(
  Curated = factor(
    ifelse(rownames(mat_huge) %in% genes_of_interest, "yes", "no"),
    levels = c("yes", "no")
  )
)
rownames(annotation_row_huge) <- rownames(mat_huge)

# Row labels: curated genes get their name; all others are blank.
# This shows roughly where each gene of interest sits without cluttering the full list.
row_labels_block <- ifelse(
  rownames(mat_huge) %in% genes_of_interest,
  rownames(mat_huge),
  ""
)

# Nonlinear color scale (same scheme as the curated heatmap)
thr_h <- 0.1
m_h <- max(abs(mat_huge), na.rm = TRUE)
m_h <- max(m_h, thr_h + 1e-6)
n_side_h <- 50
n_core_h <- 20
heat_breaks_huge <- c(
  seq(-m_h, -thr_h, length.out = n_side_h + 1),
  seq(-thr_h, thr_h, length.out = n_core_h + 1)[-1],
  seq(thr_h, m_h, length.out = n_side_h + 1)[-1]
)
heat_colors_huge <- c(
  colorRampPalette(c("blue", "#8a8aff"))(n_side_h),
  colorRampPalette(c("#eef0ff", "white", "#ffeeee"))(n_core_h),
  colorRampPalette(c("#ff8a8a", "red"))(n_side_h)
)

pheatmap(
  mat_huge,
  color = heat_colors_huge,
  breaks = heat_breaks_huge,
  cellheight = 0.5,
  fontsize_row = 6,
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  labels_row = row_labels_block,
  annotation_col = annotation_col_huge,
  annotation_row = annotation_row_huge,
  annotation_colors = list(
    Group = group_colors,
    Curated = c(yes = "black", no = "grey90")
  ),
  filename = paste0(figure_path, "Heatmap_diagonal.pdf"),
  width = 7
)

#####################################################
#####################################################
#####################################################
### Dot plot of GSEA results for CD4 vs CD8 activated vs resting comparisons,
### colored by -log10(p-adj) and sized by NES (Matching Heatmap Order/Colors)
#####################################################
#####################################################
# Reuses the global `ordered_GPs` and `highlight_colors` defined upstream.
# 2. Format the GP names in the data
df_sig <- df_sig %>%
  mutate(factor = str_replace(factor, "^F", "GP"))
# 3. Find which of our ordered GPs are actually present in this GSEA dataframe
# intersect() keeps the exact order of `ordered_GPs`, but drops any missing ones
present_GPs <- intersect(ordered_GPs, unique(df_sig$factor))
# Reverse the order so the first GP (GP56) appears at the TOP of the y-axis
y_levels <- rev(present_GPs)
# 4. Apply the strict ordering to the dataframe
df_plot <- df_sig %>%
  filter(factor %in% present_GPs) %>% # Drops any extra GPs not in our custom list
  mutate(
    pathway = factor(pathway, levels = unique(pathway)),
    factor = factor(factor, levels = y_levels),
    log10padj = -log10(padj)
  )
# 5. Extract the newly ordered levels and map them to our manual colors
y_factors <- levels(df_plot$factor)
y_colors <- highlight_colors[y_factors]
# 6. Generate the Plot
p_dotplot <- ggplot(df_plot, aes(x = pathway, y = factor)) +
  geom_point(aes(size = NES, color = log10padj)) +
  scale_color_viridis_c(name = "-log10(p-adj)") +
  scale_size(range = c(3, 10)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    # Apply our custom, ordered colors to the Y-axis text
    axis.text.y = element_text(color = y_colors, face = "bold")
  ) +
  labs(size = "NES", x = "Pathway", y = "Gene Program")
# 7. Print and Save
print(p_dotplot)
ggsave(
  filename = paste0(figure_path, "GSEA_dotplot.pdf"),
  plot = p_dotplot,
  width = 8,
  height = 10
)


#####################################################
### Shared GP grouping setup (derived from highlight_colors)
###   CD4 only  -> blue
###   CD8 only  -> darkorange2
###   both up   -> darkred
###   both down -> darkgreen
#####################################################
group_colors <- c(
  "CD4 only" = "blue",
  "CD8 only" = "darkorange2",
  "both up" = "darkred",
  "both down" = "darkgreen"
)
color_to_group <- setNames(names(group_colors), unname(group_colors))
gp_to_group <- setNames(
  color_to_group[unname(highlight_colors)],
  names(highlight_colors)
)
gp_groups <- split(
  names(gp_to_group),
  factor(gp_to_group, levels = names(group_colors))
)
lineages <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN")
L_subset <- L_pm_filtered[, GPs_of_interest, drop = FALSE]


#####################################################
### Heatmap showing average GP loadings for selected GPs
###
### Rows: the Figure 3 GPs in semantic group order.
### Cols: annotation_level2 sub-types (parent in `lineages`, >= 50 cells),
###       ordered by parent lineage.
### Color groups (row + column annotation bars) match Figure_Activation.R.
#####################################################
keep_cells <- seurat_meta_filtered$annotation_level1 %in% lineages
meta_sub <- seurat_meta_filtered[keep_cells, ]
L_keep <- L_subset[keep_cells, , drop = FALSE]

l2_counts <- table(meta_sub$annotation_level2)
l2_keep <- names(l2_counts)[l2_counts >= 50]

# Drop the "P" cluster and any "w..." clusters (wM, wW, etc.) across all
# lineages — matches the filter used in Figure_CD4.R. Strips an
# arbitrary lineage prefix ("CD4.", "CD8_", "Treg.", ...) before testing.
l2_stripped <- sub("^[^._]+[._]", "", l2_keep)
exclude_l2 <- l2_stripped == "P" |
  grepl("^w", l2_stripped, ignore.case = TRUE) |
  grepl("[._]w", l2_keep, ignore.case = TRUE)
l2_keep <- l2_keep[!exclude_l2]

# Row order: GPs in group order (so gaps_row aligns with gp_groups)
gp_row_order <- unlist(gp_groups, use.names = FALSE)

mean_mat <- vapply(
  l2_keep,
  function(l2) {
    colMeans(L_keep[meta_sub$annotation_level2 == l2, , drop = FALSE])
  },
  numeric(ncol(L_keep))
)

# Map each level2 column to its parent level1, then order columns by lineage
l2_to_l1 <- vapply(
  l2_keep,
  function(l2) {
    as.character(meta_sub$annotation_level1[meta_sub$annotation_level2 == l2][
      1
    ])
  },
  character(1)
)
col_order <- order(match(l2_to_l1, lineages), l2_keep)
mean_mat <- mean_mat[gp_row_order, col_order]
l2_to_l1 <- l2_to_l1[col_order]

immgen_cols <- ZemmourLib::immgent_colors

# Column annotation bar: parent level1 (uses ZemmourLib$level1)
col_anno <- data.frame(
  Lineage = factor(l2_to_l1, levels = lineages),
  row.names = colnames(mean_mat)
)
anno_colors_mean <- list(Lineage = immgen_cols$level1[lineages])

# Per-label text colors
#   rows: GP group (matches Figure_Activation.R)
#   cols: level2 colors from ZemmourLib (fall back to black if missing)
row_label_cols <- group_colors[gp_to_group[gp_row_order]]
col_label_cols <- immgen_cols$level2[colnames(mean_mat)]
col_label_cols[is.na(col_label_cols)] <- "black"

ph <- pheatmap(
  mean_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(200),
  annotation_col = col_anno,
  annotation_colors = anno_colors_mean,
  gaps_row = head(cumsum(lengths(gp_groups)), -1),
  gaps_col = head(cumsum(rle(l2_to_l1)$lengths), -1),
  main = "Average loading of Figure 3 GPs per Level-2 sub-lineage",
  silent = TRUE
)

# Recolor row + column label text (pheatmap doesn't expose per-label colors)
row_idx <- which(ph$gtable$layout$name == "row_names")
col_idx <- which(ph$gtable$layout$name == "col_names")
ph$gtable$grobs[[row_idx]]$gp$col <- row_label_cols
ph$gtable$grobs[[col_idx]]$gp$col <- col_label_cols

pdf(
  paste0(figure_path, "MeanLoading_Figure3_GPs_by_Level2.pdf"),
  width = 11,
  height = 5.5
)
grid::grid.draw(ph$gtable)
invisible(dev.off())


#####################################################
### log2FC heatmap of activated CD4 + CD8 cells, grouped by
### `condition_detailed_simplified`, ordered by `condition_broad`.
###
### Rows: Figure 3 GPs in semantic group order.
### Cols: condition_detailed_simplified categories with >= min_cells
###       activated CD4 AND >= min_cells activated CD8.
###       Ordered by parent condition_broad (majority vote per category).
### Reference: per-GP MEAN across all CD4/CD8 cells (resting + activated,
### every condition). Capped at +-cap_lfc.
#####################################################
act_keep <- seurat_meta_filtered$annotation_level1 %in%
  c("CD4", "CD8") &
  seurat_meta_filtered$annotation_level2_group == "activated"
meta_act <- seurat_meta_filtered[act_keep, ]
L_act <- L_subset[act_keep, , drop = FALSE]

min_cells_cond <- 50
cd_lin <- table(
  meta_act$condition_detailed_simplified,
  meta_act$annotation_level1
)
cond_keep <- rownames(cd_lin)[
  cd_lin[, "CD4"] >= min_cells_cond & cd_lin[, "CD8"] >= min_cells_cond
]

# Majority-vote condition_broad per condition_detailed_simplified category.
cd_br <- table(
  meta_act$condition_detailed_simplified,
  meta_act$condition_broad
)
cd_to_broad <- setNames(
  colnames(cd_br)[apply(cd_br, 1, which.max)],
  rownames(cd_br)
)[cond_keep]

# Column order: `healthy` broad first (with `baseline` as its first
# condition); other broads / conditions fall through to alphabetical.
broad_rank <- ifelse(cd_to_broad == "healthy", 0L, 1L)
within_broad_rank <- ifelse(
  cd_to_broad == "healthy" & cond_keep == "baseline",
  0L,
  1L
)
col_order_cond <- order(broad_rank, cd_to_broad, within_broad_rank, cond_keep)
cond_keep <- cond_keep[col_order_cond]
cd_to_broad <- cd_to_broad[cond_keep]

mean_mat_cond <- vapply(
  cond_keep,
  function(cond) {
    colMeans(
      L_act[meta_act$condition_detailed_simplified == cond, , drop = FALSE]
    )
  },
  numeric(ncol(L_act))
)
mean_mat_cond <- mean_mat_cond[gp_row_order, , drop = FALSE]

broad_levels <- unique(cd_to_broad)
col_anno_cond <- data.frame(
  condition_broad = factor(cd_to_broad, levels = broad_levels),
  row.names = colnames(mean_mat_cond)
)

row_label_cols <- group_colors[gp_to_group[gp_row_order]]

# log2FC vs per-GP MEAN across all CD4/CD8 cells
pc_lfc <- 1e-10
cap_lfc <- 2
cd4cd8_idx <- seurat_meta_filtered$annotation_level1 %in% c("CD4", "CD8")
L_cd4cd8 <- L_subset[cd4cd8_idx, , drop = FALSE]
mu_lfc_mean <- colMeans(L_cd4cd8, na.rm = TRUE)[rownames(mean_mat_cond)]
lfc_mat_mean <- log2((mean_mat_cond + pc_lfc) / (mu_lfc_mean + pc_lfc))
lfc_mat_mean <- pmax(pmin(lfc_mat_mean, cap_lfc), -cap_lfc)

ph_cond_lfc_mean <- pheatmap(
  lfc_mat_mean,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#7A0177", "black", "#FFD700"))(101),
  breaks = seq(-cap_lfc, cap_lfc, length.out = 102),
  annotation_col = col_anno_cond,
  gaps_row = head(cumsum(lengths(gp_groups)), -1),
  gaps_col = head(cumsum(rle(as.character(cd_to_broad))$lengths), -1),
  main = "log2FC vs per-GP MEAN across all CD4/CD8 (activated CD4+CD8 by condition_detailed_simplified)",
  silent = TRUE
)
row_idx_lfc_m <- which(ph_cond_lfc_mean$gtable$layout$name == "row_names")
ph_cond_lfc_mean$gtable$grobs[[row_idx_lfc_m]]$gp$col <- row_label_cols

pdf(
  paste0(figure_path, "Log2FC_ActivatedCD4CD8_by_condition_meanref.pdf"),
  width = max(8, 0.18 * ncol(lfc_mat_mean) + 4),
  height = 6
)
grid::grid.draw(ph_cond_lfc_mean$gtable)
invisible(dev.off())


#####################################################
### TF-GP network plot (ported from Figure_TF_and_Activation.R).
###
### Builds the TF list from GO:0003700 (DNA-binding TF activity) plus the
### Tox family, filters TFs whose max normalized score in any highlighted
### GP exceeds `tf_gp_threshold`, and draws a bipartite TF-GP network with
### barycenter-optimized layout. Top genes per GP are listed next to each
### GP node in italics.
#####################################################
mm <- org.Mm.eg.db
go2eg <- as.list(org.Mm.egGO2ALLEGS)
tf_symbols <- AnnotationDbi::select(
  mm,
  keys = unique(unlist(go2eg)),
  columns = "SYMBOL",
  keytype = "ENTREZID"
)
tf <- c(
  sort(tf_symbols$SYMBOL[
    tf_symbols$ENTREZID %in% unique(go2eg[["GO:0003700"]])
  ]),
  "Tox",
  "Tox2",
  "Tox3",
  "Tox4"
) %>%
  sort() %>%
  unique()

# Reuse F_pm_filtered_norm from the earlier network section; restrict to the
# highlighted GPs.
F_sub_tf <- F_pm_filtered_norm[, GPs_of_interest, drop = FALSE]

# Filter TFs: keep those with max normalized score > threshold in any
# highlighted GP.
tf_gp_threshold <- 0.25
tf_in_F <- intersect(tf, rownames(F_sub_tf))
tf_max_score <- apply(F_sub_tf[tf_in_F, , drop = FALSE], 1, max, na.rm = TRUE)
selected_tfs <- sort(names(tf_max_score)[tf_max_score > tf_gp_threshold])
message(sprintf(
  "TF filter: %d / %d TFs (from %d in F) clear score > %.2f in any of the %d GPs",
  length(selected_tfs),
  length(tf),
  length(tf_in_F),
  tf_gp_threshold,
  length(GPs_of_interest)
))

# Top-to-bottom group order in the network: both-up -> CD8-only -> both-down
# -> CD4-only. Matches the prior TF network output.
gp_color_group_order <- c("darkred", "darkorange2", "darkgreen", "blue")

# Barycenter-based bipartite ordering. Keeps each color group contiguous;
# inside each group, TF and GP positions are iteratively reordered so edges
# straighten out.
optimize_bipartite_order <- function(
  adj,
  gp_groups = NULL,
  gp_group_order = NULL,
  n_iter = 12
) {
  if (is.null(gp_groups)) {
    gp_groups <- setNames(rep("_all", ncol(adj)), colnames(adj))
  } else {
    missing_g <- setdiff(colnames(adj), names(gp_groups))
    if (length(missing_g)) gp_groups[missing_g] <- "_unassigned"
  }
  groups_present <- unique(gp_groups[colnames(adj)])
  if (!is.null(gp_group_order)) {
    groups_present <- c(
      intersect(gp_group_order, groups_present),
      setdiff(groups_present, gp_group_order)
    )
  }
  reorder_gps_within_groups <- function(adj_mat) {
    tf_pos <- setNames(seq_len(nrow(adj_mat)), rownames(adj_mat))
    unlist(lapply(groups_present, function(g) {
      gps_g <- colnames(adj_mat)[gp_groups[colnames(adj_mat)] == g]
      if (length(gps_g) <= 1) {
        return(gps_g)
      }
      sub <- adj_mat[, gps_g, drop = FALSE]
      bary <- apply(sub, 2, function(col) {
        if (sum(col) == 0) {
          return(mean(tf_pos))
        }
        sum(tf_pos * col) / sum(col)
      })
      names(sort(bary))
    }))
  }
  reorder_tfs <- function(adj_mat) {
    gp_pos <- setNames(seq_len(ncol(adj_mat)), colnames(adj_mat))
    bary <- apply(adj_mat, 1, function(row) {
      if (sum(row) == 0) {
        return(mean(gp_pos))
      }
      sum(gp_pos * row) / sum(row)
    })
    names(sort(bary))
  }
  adj <- adj[, reorder_gps_within_groups(adj), drop = FALSE]
  for (i in seq_len(n_iter)) {
    new_tf <- reorder_tfs(adj)
    adj <- adj[new_tf, , drop = FALSE]
    new_gp <- reorder_gps_within_groups(adj)
    if (
      identical(new_gp, colnames(adj)) &&
        identical(new_tf, rownames(adj))
    ) {
      adj <- adj[, new_gp, drop = FALSE]
      break
    }
    adj <- adj[, new_gp, drop = FALSE]
  }
  list(tf_order = rownames(adj), gp_order = colnames(adj))
}

plot_tf_gp_network_v2 <- function(
  F,
  selected_tfs,
  tf_gp_threshold = 0.25,
  top_genes_per_gp = 5,
  tf_colors = NULL,
  gp_colors = NULL,
  gp_group_order = NULL,
  default_gp_fill = "#4daf4a",
  optimize_layout = TRUE,
  barycenter_iter = 12,
  node_size_tf = 5,
  node_size_gp = 3.5,
  edge_alpha = 0.55,
  edge_curvature = 0.2,
  label_size_tf = 3.5,
  label_size_gp = 2.8,
  label_size_gene = 2.3,
  gp_spacing = 1.5
) {
  if (is.null(colnames(F))) {
    colnames(F) <- paste0("GP", seq_len(ncol(F)))
  }
  selected_tfs <- intersect(selected_tfs, rownames(F))
  if (!length(selected_tfs)) {
    stop("No selected TFs found in F.")
  }

  tf_gp_edges <- do.call(
    rbind,
    lapply(selected_tfs, function(tf_name) {
      vals <- setNames(as.numeric(F[tf_name, ]), colnames(F))
      idx <- which(is.finite(vals) & vals > tf_gp_threshold)
      if (!length(idx)) {
        return(NULL)
      }
      data.frame(
        from = tf_name,
        to = names(vals)[idx],
        weight = vals[idx],
        stringsAsFactors = FALSE
      )
    })
  )
  if (is.null(tf_gp_edges)) {
    stop("No TF-GP edges above threshold.")
  }

  connected_gps <- unique(tf_gp_edges$to)
  connected_gps <- connected_gps[order(as.numeric(sub(
    "GP",
    "",
    connected_gps
  )))]
  selected_tfs <- intersect(selected_tfs, unique(tf_gp_edges$from))

  adj <- matrix(
    0,
    nrow = length(selected_tfs),
    ncol = length(connected_gps),
    dimnames = list(selected_tfs, connected_gps)
  )
  adj[cbind(tf_gp_edges$from, tf_gp_edges$to)] <- tf_gp_edges$weight

  if (optimize_layout) {
    gp_groups <- if (!is.null(gp_colors)) gp_colors else NULL
    ord <- optimize_bipartite_order(
      adj,
      gp_groups = gp_groups,
      gp_group_order = gp_group_order,
      n_iter = barycenter_iter
    )
    selected_tfs <- ord$tf_order
    connected_gps <- ord$gp_order
  }
  n_gps <- length(connected_gps)
  n_tfs <- length(selected_tfs)

  gp_top_genes <- setNames(
    lapply(connected_gps, function(gp) {
      vals <- setNames(as.numeric(F[, gp]), rownames(F))
      vals <- sort(
        vals[is.finite(vals) & vals > 0 & !names(vals) %in% selected_tfs],
        decreasing = TRUE
      )
      names(vals)[seq_len(min(top_genes_per_gp, length(vals)))]
    }),
    connected_gps
  )

  gp_y <- rev(seq_len(n_gps)) * gp_spacing
  gp_fill_vec <- if (!is.null(gp_colors)) {
    out <- gp_colors[connected_gps]
    out[is.na(out)] <- default_gp_fill
    out
  } else {
    setNames(rep(default_gp_fill, n_gps), connected_gps)
  }
  gp_nodes <- data.frame(
    node = connected_gps,
    x = 1,
    y = gp_y,
    fill_color = unname(gp_fill_vec),
    stringsAsFactors = FALSE
  )
  tf_nodes <- data.frame(
    node = selected_tfs,
    x = 0,
    y = seq(max(gp_y), min(gp_y), length.out = n_tfs),
    stringsAsFactors = FALSE
  )
  if (is.null(tf_colors)) {
    tf_colors <- setNames(scales::hue_pal()(n_tfs), selected_tfs)
  }
  tf_nodes$fill_color <- unname(tf_colors[selected_tfs])

  edge_df <- tf_gp_edges %>%
    left_join(tf_nodes %>% rename(x0 = x, y0 = y), by = c("from" = "node")) %>%
    left_join(gp_nodes %>% rename(x1 = x, y1 = y), by = c("to" = "node")) %>%
    filter(is.finite(x0), is.finite(y0), is.finite(x1), is.finite(y1))

  gene_label_df <- do.call(
    rbind,
    lapply(connected_gps, function(gp) {
      data.frame(
        gp = gp,
        label = paste(gp_top_genes[[gp]], collapse = ", "),
        x = 1.05,
        y = gp_nodes$y[gp_nodes$node == gp],
        stringsAsFactors = FALSE
      )
    })
  )

  ggplot() +
    geom_curve(
      data = edge_df,
      aes(
        x = x0,
        y = y0,
        xend = x1,
        yend = y1,
        color = from,
        linewidth = weight
      ),
      curvature = edge_curvature,
      alpha = edge_alpha,
      arrow = arrow(length = unit(0.07, "inches"), type = "closed")
    ) +
    geom_point(
      data = gp_nodes,
      aes(x = x, y = y, fill = fill_color),
      shape = 22,
      color = "black",
      size = node_size_gp
    ) +
    geom_point(
      data = tf_nodes,
      aes(x = x, y = y, fill = fill_color),
      shape = 21,
      color = "black",
      size = node_size_tf,
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    geom_text(
      data = tf_nodes,
      aes(x = x, y = y, label = node),
      hjust = 1.3,
      size = label_size_tf,
      fontface = "bold"
    ) +
    geom_text(
      data = gp_nodes,
      aes(x = x, y = y, label = node),
      hjust = 1.3,
      size = label_size_gp,
      color = "grey25"
    ) +
    geom_text(
      data = gene_label_df,
      aes(x = x, y = y, label = label),
      hjust = 0,
      size = label_size_gene,
      color = "grey40",
      fontface = "italic"
    ) +
    scale_color_manual(values = tf_colors, name = "TF") +
    scale_linewidth_continuous(range = c(0.3, 1.8), guide = "none") +
    coord_cartesian(xlim = c(-0.25, 4.5), clip = "off") +
    theme_void(base_size = 12) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(20, 20, 20, 80)
    )
}

tf_network_plot <- plot_tf_gp_network_v2(
  F = F_sub_tf,
  selected_tfs = selected_tfs,
  tf_gp_threshold = tf_gp_threshold,
  top_genes_per_gp = 5,
  gp_colors = highlight_colors,
  gp_group_order = gp_color_group_order,
  optimize_layout = TRUE,
  barycenter_iter = 12,
  gp_spacing = 1.5,
  node_size_tf = 6,
  node_size_gp = 5,
  label_size_tf = 4.5,
  label_size_gp = 4,
  label_size_gene = 3.4
)
plot_height_tf <- min(
  60,
  max(12, length(selected_tfs) * 0.35, length(GPs_of_interest) * 1.5 * 0.55 + 2)
)
ggsave(
  filename = paste0(figure_path, "TF_GP_network_activation.pdf"),
  plot = tf_network_plot,
  width = 18,
  height = plot_height_tf,
  limitsize = FALSE
)


#####################################################
### Supplemental centered TF-GP network for Gata3, Rorc, and Tbx21
###
### Uses the per-GP normalized F matrix. Edges connect each TF to every GP
### with abs(normalized F value) >= 0.1; edge color indicates sign and edge
### labels show the signed normalized value.
#####################################################
tf_focus <- c("Gata3", "Rorc", "Tbx21")
tf_focus_threshold <- 0.1
missing_tf_focus <- setdiff(tf_focus, rownames(F_pm_filtered_norm))
if (length(missing_tf_focus) > 0) {
  stop(
    "Missing focus TFs from F_pm_filtered_norm: ",
    paste(missing_tf_focus, collapse = ", ")
  )
}

tf_focus_edges <- do.call(
  rbind,
  lapply(tf_focus, function(tf_name) {
    vals <- setNames(
      as.numeric(F_pm_filtered_norm[tf_name, ]),
      colnames(F_pm_filtered_norm)
    )
    idx <- which(is.finite(vals) & abs(vals) >= tf_focus_threshold)
    data.frame(
      TF = tf_name,
      GP = names(vals)[idx],
      value = vals[idx],
      stringsAsFactors = FALSE
    )
  })
) %>%
  mutate(
    GP_number = as.numeric(sub("^GP", "", GP)),
    edge_sign = if_else(value < 0, "Negative", "Positive"),
    abs_value = abs(value),
    edge_label = sprintf("%+.2f", value)
  ) %>%
  arrange(match(TF, tf_focus), GP_number)

write.csv(
  tf_focus_edges %>%
    select(TF, GP, value, edge_sign, abs_value),
  file = paste0(figure_path, "TF_GP_network_Gata3_Rorc_Tbx21_edges.csv"),
  row.names = FALSE
)

tf_focus_nodes <- data.frame(
  name = tf_focus,
  x = c(0, 2.7, 1.65),
  y = c(0, -2.75, 2.35),
  node_type = "TF",
  stringsAsFactors = FALSE
)
tf_focus_colors <- c(
  Gata3 = "#E41A1C",
  Rorc = "#1F78B4",
  Tbx21 = "#33A02C"
)

gp_tf_membership <- tf_focus_edges %>%
  group_by(GP) %>%
  summarise(
    tf_members = list(sort(unique(TF))),
    degree_tf = n_distinct(TF),
    .groups = "drop"
  )

shared_gp_nodes <- gp_tf_membership %>%
  filter(degree_tf > 1) %>%
  rowwise() %>%
  mutate(
    x = mean(tf_focus_nodes$x[match(tf_members, tf_focus_nodes$name)]),
    y = mean(tf_focus_nodes$y[match(tf_members, tf_focus_nodes$name)])
  ) %>%
  ungroup() %>%
  transmute(name = GP, x, y, node_type = "GP")

tf_arc_range <- list(
  Gata3 = c(170, -105),
  Rorc = c(205, -20),
  Tbx21 = c(165, -15)
)
tf_arc_radius <- c(Gata3 = 1.65, Rorc = 1.35, Tbx21 = 1.30)
make_tf_arc_nodes <- function(tf_name, gp_names) {
  if (length(gp_names) == 0) {
    return(NULL)
  }
  center <- tf_focus_nodes[tf_focus_nodes$name == tf_name, ]
  angles <- seq(
    tf_arc_range[[tf_name]][1],
    tf_arc_range[[tf_name]][2],
    length.out = length(gp_names)
  ) *
    pi /
    180
  radius <- tf_arc_radius[[tf_name]]
  data.frame(
    name = gp_names,
    x = center$x + radius * cos(angles),
    y = center$y + radius * sin(angles),
    node_type = "GP",
    stringsAsFactors = FALSE
  )
}

exclusive_gp_nodes <- do.call(
  rbind,
  lapply(tf_focus, function(tf_name) {
    gp_names <- gp_tf_membership %>%
      filter(
        degree_tf == 1,
        vapply(tf_members, identical, logical(1), tf_name)
      ) %>%
      pull(GP)
    gp_names <- gp_names[order(as.numeric(sub("^GP", "", gp_names)))]
    make_tf_arc_nodes(tf_name, gp_names)
  })
)

tf_focus_plot_nodes <- bind_rows(
  tf_focus_nodes,
  shared_gp_nodes,
  exclusive_gp_nodes
)

tf_focus_plot_edges <- tf_focus_edges %>%
  left_join(
    tf_focus_plot_nodes %>%
      select(name, x, y) %>%
      rename(x0 = x, y0 = y),
    by = c("TF" = "name")
  ) %>%
  left_join(
    tf_focus_plot_nodes %>%
      select(name, x, y) %>%
      rename(x1 = x, y1 = y),
    by = c("GP" = "name")
  ) %>%
  mutate(
    label_x = x0 + 0.55 * (x1 - x0),
    label_y = y0 + 0.55 * (y1 - y0),
    label_angle = atan2(y1 - y0, x1 - x0) * 180 / pi,
    label_angle = case_when(
      label_angle > 90 ~ label_angle - 180,
      label_angle < -90 ~ label_angle + 180,
      TRUE ~ label_angle
    )
  )

tf_focus_network_plot <- ggplot() +
  geom_segment(
    data = tf_focus_plot_edges,
    aes(
      x = x0,
      y = y0,
      xend = x1,
      yend = y1,
      color = edge_sign,
      linewidth = abs_value
    ),
    alpha = 0.65,
    lineend = "round"
  ) +
  geom_text(
    data = tf_focus_plot_edges,
    aes(
      x = label_x,
      y = label_y,
      label = edge_label,
      angle = label_angle
    ),
    size = 3.0,
    color = "grey20"
  ) +
  geom_point(
    data = tf_focus_plot_nodes %>% filter(node_type == "GP"),
    aes(x = x, y = y),
    shape = 21,
    size = 4.0,
    fill = "grey78",
    color = "grey45",
    stroke = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = tf_focus_plot_nodes %>% filter(node_type == "GP"),
    aes(x = x, y = y, label = name),
    size = 3.4,
    color = "grey20",
    max.overlaps = Inf,
    min.segment.length = Inf,
    box.padding = 0.2,
    point.padding = 0.25
  ) +
  geom_label(
    data = tf_focus_plot_nodes %>% filter(node_type == "TF"),
    aes(x = x, y = y, label = name, fill = name),
    color = "white",
    fontface = "bold",
    size = 4.6,
    label.padding = unit(0.22, "lines"),
    label.r = unit(0.16, "lines"),
    linewidth = 0,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(Negative = "#2B6CB0", Positive = "#D62728"),
    name = "Edge sign"
  ) +
  scale_fill_manual(values = tf_focus_colors) +
  scale_linewidth_continuous(range = c(0.45, 3.2), guide = "none") +
  coord_equal(clip = "off") +
  labs(
    title = paste0(
      "TF <-> GP network (",
      length(tf_focus),
      " TFs, ",
      nrow(tf_focus_edges),
      " edges)"
    )
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(15, 25, 15, 25)
  )

print(tf_focus_network_plot)
ggsave(
  filename = paste0(figure_path, "TF_GP_network_Gata3_Rorc_Tbx21.pdf"),
  plot = tf_focus_network_plot,
  width = 10,
  height = 8.5
)


#####################################################
### GP57 high-loading proportion in activated CD4/CD8 cells
###
### For each lineage, compare condition_broad == "cancer" against all
### other conditions combined.
#####################################################
gp57_condition_df <- data.frame(
  lineage = seurat_meta_filtered$annotation_level1,
  activation_group = seurat_meta_filtered$annotation_level2_group,
  condition_broad = as.character(seurat_meta_filtered$condition_broad),
  gp57_loading = L_pm_filtered[, "GP57"],
  stringsAsFactors = FALSE
) %>%
  filter(
    lineage %in% c("CD8", "CD4"),
    activation_group == "activated"
  ) %>%
  mutate(
    lineage = factor(lineage, levels = c("CD8", "CD4")),
    condition_group = if_else(
      condition_broad == "cancer",
      "Cancer",
      "Other conditions"
    ),
    condition_group = factor(
      condition_group,
      levels = c("Cancer", "Other conditions")
    ),
    gp57_high = gp57_loading > 0.1
  )

gp57_condition_summary <- gp57_condition_df %>%
  group_by(lineage, condition_group) %>%
  summarise(
    n_cells = n(),
    n_gp57_high = sum(gp57_high),
    proportion_gp57_high = mean(gp57_high),
    .groups = "drop"
  )

write.csv(
  gp57_condition_summary,
  file = paste0(figure_path, "GP57_condition_broad_summary.csv"),
  row.names = FALSE
)

gp57_ymax <- max(gp57_condition_summary$proportion_gp57_high, na.rm = TRUE)
gp57_ymax <- ifelse(gp57_ymax > 0, gp57_ymax * 1.25, 0.05)

gp57_condition_barplot <- ggplot(
  gp57_condition_summary,
  aes(x = lineage, y = proportion_gp57_high, fill = condition_group)
) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62,
    color = "grey20",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      label = paste0(
        scales::percent(proportion_gp57_high, accuracy = 0.1),
        "\n",
        n_gp57_high,
        "/",
        n_cells
      )
    ),
    position = position_dodge(width = 0.72),
    vjust = -0.25,
    size = 3.4,
    lineheight = 0.9
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, gp57_ymax),
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_fill_manual(
    values = c("Cancer" = "#C44E52", "Other conditions" = "#4C72B0")
  ) +
  labs(
    x = NULL,
    y = "Proportion of activated cells",
    fill = "Condition",
    title = "GP57 loading > 0.1 in activated CD8 and CD4 cells"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(gp57_condition_barplot)
ggsave(
  filename = paste0(figure_path, "GP57_condition_broad_barplot.pdf"),
  plot = gp57_condition_barplot,
  width = 5.4,
  height = 4.2
)


#####################################################
### GP79 high-loading proportion across conditions in activated CD4 cells
###
### Summary CSV keeps every condition_detailed_simplified category. The PDF
### filters out very small categories so one-cell conditions do not dominate
### the visual comparison.
#####################################################
min_cells_gp79_condition <- 50

gp79_cd4_condition_df <- data.frame(
  condition_broad = as.character(seurat_meta_filtered$condition_broad),
  condition_detailed_simplified = as.character(
    seurat_meta_filtered$condition_detailed_simplified
  ),
  lineage = seurat_meta_filtered$annotation_level1,
  activation_group = seurat_meta_filtered$annotation_level2_group,
  gp79_loading = L_pm_filtered[, "GP79"],
  stringsAsFactors = FALSE
) %>%
  filter(
    lineage == "CD4",
    activation_group == "activated",
    !is.na(condition_detailed_simplified),
    condition_detailed_simplified != ""
  ) %>%
  mutate(gp79_high = gp79_loading > 0.1)

gp79_cd4_condition_summary <- gp79_cd4_condition_df %>%
  group_by(condition_detailed_simplified) %>%
  summarise(
    condition_broad = names(sort(table(condition_broad), decreasing = TRUE))[1],
    n_cells = n(),
    n_gp79_high = sum(gp79_high),
    proportion_gp79_high = mean(gp79_high),
    .groups = "drop"
  ) %>%
  mutate(
    included_in_plot = n_cells >= min_cells_gp79_condition
  ) %>%
  arrange(desc(proportion_gp79_high), condition_detailed_simplified)

write.csv(
  gp79_cd4_condition_summary,
  file = paste0(
    figure_path,
    "GP79_activated_CD4_condition_detailed_summary.csv"
  ),
  row.names = FALSE
)

gp79_cd4_condition_plot_df <- gp79_cd4_condition_summary %>%
  filter(included_in_plot) %>%
  arrange(proportion_gp79_high, condition_detailed_simplified) %>%
  mutate(
    condition_label = factor(
      condition_detailed_simplified,
      levels = condition_detailed_simplified
    )
  )

gp79_broad_levels <- sort(unique(gp79_cd4_condition_plot_df$condition_broad))
gp79_broad_colors <- setNames(
  scales::hue_pal()(length(gp79_broad_levels)),
  gp79_broad_levels
)
gp79_xmax <- max(gp79_cd4_condition_plot_df$proportion_gp79_high, na.rm = TRUE)
gp79_xmax <- ifelse(gp79_xmax > 0, gp79_xmax, 0.05)
gp79_label_pad <- gp79_xmax * 0.015
gp79_plot_width <- 9
gp79_plot_height <- max(7, 0.18 * nrow(gp79_cd4_condition_plot_df) + 2)

gp79_cd4_condition_barplot <- ggplot(
  gp79_cd4_condition_plot_df,
  aes(
    x = proportion_gp79_high,
    y = condition_label,
    fill = condition_broad
  )
) +
  geom_col(width = 0.75, color = "grey25", linewidth = 0.2) +
  geom_text(
    aes(
      x = proportion_gp79_high + gp79_label_pad,
      label = scales::percent(proportion_gp79_high, accuracy = 1)
    ),
    hjust = 0,
    size = 2.7
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 10),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = gp79_broad_colors) +
  coord_cartesian(xlim = c(0, gp79_xmax * 1.15), clip = "off") +
  labs(
    x = "Proportion of activated CD4 cells with GP79 loading > 0.1",
    y = NULL,
    fill = "Condition broad",
    title = "GP79-high fraction across activated CD4 conditions",
    subtitle = paste0(
      "condition_detailed_simplified categories with >= ",
      min_cells_gp79_condition,
      " activated CD4 cells"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(10, 45, 10, 10)
  )

print(gp79_cd4_condition_barplot)
ggsave(
  filename = paste0(
    figure_path,
    "GP79_activated_CD4_condition_detailed_barplot.pdf"
  ),
  plot = gp79_cd4_condition_barplot,
  width = gp79_plot_width,
  height = gp79_plot_height,
  limitsize = FALSE
)
