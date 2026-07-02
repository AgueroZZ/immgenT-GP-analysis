library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(fastTopics)
library(dplyr)
library(pheatmap)
library(scales)
library(scattermore)
library(tibble)
library(pROC)
library(parallel)
library(doParallel)
library(tidyr)

gp_label <- function(x) sub("^K(\\d+)$", "GP\\1", x)

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
figure_path <- "figures/Batch_Effects/"
source(paste0(code_path, "ROC.R"))
source(paste0(code_path, "filtering_membership.R"))


#####################################################
#####################################################
#####################################################
### Loading Data
#####################################################
#####################################################
load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
seurat_meta <- readRDS(paste0(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
seurat_meta_filtered_spleen <- seurat_meta_filtered %>%
  filter(spleen_standard == TRUE)


#####################################################
#####################################################
#####################################################
### Plot GP loading by IGT (Boxplot, Mean+SE, or Activated Cell Count)
#####################################################
#####################################################
plot_gp_loading <- function(
  L_pm_filtered,
  seurat_meta,
  gp_name,
  igt_vector,
  igt_col = "IGT", # Added this to make it flexible
  plot_type = "boxplot",
  loading_threshold = 0.1
) {
  # 1. Match Cells
  common_cells <- intersect(rownames(L_pm_filtered), rownames(seurat_meta))

  if (length(common_cells) == 0) {
    stop("No matching Cell IDs found between Loadings and Metadata!")
  }

  # 2. Build Dataframe
  # Using [[ ]] helps catch dynamic column names
  plot_df <- data.frame(
    cellID = common_cells,
    loading = L_pm_filtered[common_cells, gp_name],
    IGT = seurat_meta[common_cells, igt_col]
  )

  # 3. Filter by IGT
  plot_df <- plot_df %>% filter(IGT %in% igt_vector)

  # DEBUG PRINT: Check if rows exist after filtering
  message(paste("Rows to plot:", nrow(plot_df)))
  if (nrow(plot_df) == 0) {
    message(
      "Warning: Your igt_vector doesn't match values in the metadata column."
    )
    return(NULL)
  }

  # 4. Plotting
  p <- ggplot(
    plot_df,
    aes(x = factor(IGT, levels = igt_vector), y = loading, fill = IGT)
  ) +
    theme_cowplot() +
    labs(
      title = paste("Activity of", gp_label(gp_name)),
      y = "Loading Value",
      x = "IGT"
    ) +
    theme(legend.position = "none")

  if (plot_type == "mean_se") {
    p <- ggplot(
      plot_df,
      aes(x = factor(IGT, levels = igt_vector), y = loading)
    ) +
      stat_summary(fun = "mean", geom = "point", size = 3) +
      stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
      theme_cowplot() +
      labs(title = paste("Mean Loading:", gp_label(gp_name)))
  } else if (plot_type == "boxplot") {
    p <- p + geom_boxplot(outlier.size = 0.5, alpha = 0.7)
  } else if (plot_type == "activated_count") {
    summary_df <- plot_df %>%
      group_by(IGT) %>%
      summarise(n_activated = sum(loading >= loading_threshold))

    p <- ggplot(
      summary_df,
      aes(x = factor(IGT, levels = igt_vector), y = n_activated, fill = IGT)
    ) +
      geom_col() +
      theme_cowplot() +
      geom_text(aes(label = n_activated), vjust = -0.5) +
      labs(
        title = paste("Cells Activated in", gp_label(gp_name)),
        subtitle = paste("Threshold >=", loading_threshold)
      )
  }

  return(p)
}


# select IGTs that at least 500 cells in the spleen standard subset
selected_igts <- names(table(seurat_meta_filtered_spleen$IGT))[
  table(seurat_meta_filtered_spleen$IGT) >= 500
]


#####################################################
#####################################################
#####################################################
### Batch Effect on specific organ through Heatmap
#####################################################
#####################################################
plot_gp_heatmap <- function(
  L_pm_filtered,
  seurat_meta,
  gp_vector,
  igt_vector,
  igt_col = "IGT",
  igt_stat = c("median", "mean"),
  top_n_var = NULL, # If set (e.g., 50), keeps most variable GPs
  floor_threshold = 0, # Values below this become 0 (White)
  scale_row = FALSE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color_scheme = c("white", "red"),
  main_title = NULL,
  ...
) {
  igt_stat <- match.arg(igt_stat)

  # 1. Align cells
  common_cells <- intersect(rownames(L_pm_filtered), rownames(seurat_meta))

  # 2. Extract and Aggregate
  df_subset <- data.frame(L_pm_filtered[common_cells, gp_vector, drop = FALSE])
  df_subset[[igt_col]] <- seurat_meta[common_cells, igt_col]

  agg_fn <- if (igt_stat == "mean") mean else median
  agg_df <- df_subset %>%
    filter(!!sym(igt_col) %in% igt_vector) %>%
    group_by(!!sym(igt_col)) %>%
    summarise(across(everything(), agg_fn, na.rm = TRUE)) %>%
    as.data.frame()

  rownames(agg_df) <- agg_df[[igt_col]]
  agg_mat <- as.matrix(agg_df[,
    -which(colnames(agg_df) == igt_col),
    drop = FALSE
  ])

  if (is.null(main_title)) {
    main_title <- paste0(
      if (igt_stat == "mean") "Mean" else "Median",
      " GP Loading"
    )
  }

  # 3. Filter by Variance (Pre-Thresholding)
  # We want to find GPs that vary significantly across IGTs
  if (!is.null(top_n_var)) {
    gp_vars <- apply(agg_mat, 2, var)
    # Sort and take top N
    top_gps <- names(sort(gp_vars, decreasing = TRUE))[
      1:min(top_n_var, length(gp_vars))
    ]
    agg_mat <- agg_mat[, top_gps, drop = FALSE]
  }

  # 4. Apply Floor Threshold
  agg_mat[agg_mat < floor_threshold] <- 0

  # 5. Transpose (GPs as rows)
  plot_mat <- t(agg_mat)
  rownames(plot_mat) <- gp_label(rownames(plot_mat))

  # 6. Optional Z-score scaling
  if (scale_row) {
    nonzero_rows <- rowSums(plot_mat) > 0
    if (any(nonzero_rows)) {
      plot_mat[nonzero_rows, ] <- t(apply(
        plot_mat[nonzero_rows, , drop = FALSE],
        1,
        scale
      ))
    }
    colnames(plot_mat) <- rownames(agg_df)
    main_title <- paste(main_title, "(Row Z-score)")
  }

  # 7. Plot
  pheatmap(
    plot_mat,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    main = main_title,
    color = colorRampPalette(color_scheme)(100),
    border_color = "white",
    ...
  )
}
selected_igts <- names(table(seurat_meta_filtered_spleen$IGT))[
  table(seurat_meta_filtered_spleen$IGT) >= 500
]
all_gps <- paste0("K", 1:200)


#####################################################
#####################################################
#####################################################
### GP Overall Loading vs. Between-IGT SD (Spleen Standard)
#####################################################
#####################################################
plot_gp_igt_dispersion <- function(
  L_pm_filtered,
  seurat_meta,
  igt_col = "IGT",
  x_stat = c("mean", "median", "max_igt"),
  igt_stat = c("mean", "median"),
  y_type = c("sd", "cv", "var"),
  log = "",
  min_cells_per_igt = NULL,
  top_n_label = 10
) {
  x_stat <- match.arg(x_stat)
  igt_stat <- match.arg(igt_stat)
  y_type <- match.arg(y_type)

  common_cells <- intersect(rownames(L_pm_filtered), rownames(seurat_meta))
  L_sub <- L_pm_filtered[common_cells, ]
  igt_vec <- seurat_meta[common_cells, igt_col]

  if (!is.null(min_cells_per_igt)) {
    valid_igts <- names(table(igt_vec))[table(igt_vec) >= min_cells_per_igt]
    keep <- igt_vec %in% valid_igts
    L_sub <- L_sub[keep, ]
    igt_vec <- igt_vec[keep]
  }

  # Per-IGT mean or median for each GP  ->  IGTs x GPs matrix
  igt_levels <- unique(igt_vec)
  igt_mat <- do.call(
    rbind,
    lapply(igt_levels, function(igt) {
      cells_i <- igt_vec == igt
      if (igt_stat == "mean") {
        colMeans(L_sub[cells_i, , drop = FALSE])
      } else {
        apply(L_sub[cells_i, , drop = FALSE], 2, median)
      }
    })
  )

  gp_igt_sd <- apply(igt_mat, 2, sd)
  gp_igt_var <- apply(igt_mat, 2, var)
  gp_overall <- switch(
    x_stat,
    mean = colMeans(igt_mat),
    median = apply(igt_mat, 2, median),
    max_igt = apply(igt_mat, 2, max)
  )

  gp_y <- switch(
    y_type,
    sd = gp_igt_sd,
    cv = gp_igt_sd / gp_overall,
    var = gp_igt_var
  )

  gp_stats <- data.frame(GP = colnames(L_sub), x = gp_overall, y = gp_y)
  gp_stats <- gp_stats %>%
    arrange(desc(y)) %>%
    mutate(
      label = ifelse(
        row_number() <= top_n_label,
        gp_label(as.character(GP)),
        ""
      )
    )

  igt_stat_label <- if (igt_stat == "mean") "Mean" else "Median"
  x_lab <- switch(
    x_stat,
    mean = paste0("Mean of IGT ", igt_stat_label, " Loading"),
    median = paste0("Median of IGT ", igt_stat_label, " Loading"),
    max_igt = paste0("Max IGT ", igt_stat_label, " Loading")
  )
  y_lab <- switch(
    y_type,
    sd = paste0("SD of IGT ", igt_stat_label, " Loading"),
    cv = paste0("CV of IGT ", igt_stat_label, " Loading (SD / ", x_lab, ")"),
    var = paste0("Variance of IGT ", igt_stat_label, " Loading")
  )

  p <- ggplot(gp_stats, aes(x = x, y = y, label = label)) +
    geom_point(size = 1.5, alpha = 0.7, color = "steelblue") +
    geom_text_repel(
      size = 3,
      box.padding = 0.4,
      max.overlaps = Inf,
      segment.color = "grey50"
    ) +
    theme_cowplot() +
    labs(
      title = paste0(
        "GP ",
        x_lab,
        " vs. Between-IGT ",
        toupper(y_type)
      ),
      x = x_lab,
      y = y_lab
    )

  if (grepl("x", log)) {
    p <- p + scale_x_log10()
  }
  if (grepl("y", log)) {
    p <- p + scale_y_log10()
  }
  p
}
p_igt_disp <- plot_gp_igt_dispersion(
  L_pm_filtered,
  igt_stat = "mean",
  x_stat = "mean",
  # log = "x",
  seurat_meta_filtered_spleen,
  y_type = "var",
  top_n_label = 10
)
p_igt_disp
ggsave(
  paste0(figure_path, "GP_IGT_Dispersion_Spleen.pdf"),
  width = 6,
  height = 5,
  dpi = 300
)


#####################################################
#####################################################
#####################################################
### GP Heatmap: Top GPs by Variance of IGT Mean Loading
#####################################################
#####################################################
spleen_cells_act <- intersect(
  rownames(L_pm_filtered),
  rownames(seurat_meta_filtered_spleen)
)
igt_vec_act <- seurat_meta_filtered_spleen[spleen_cells_act, "IGT"]
igt_mean_mat_act <- do.call(
  rbind,
  lapply(selected_igts, function(igt) {
    cells_i <- spleen_cells_act[igt_vec_act == igt]
    colMeans(L_pm_filtered[cells_i, all_gps, drop = FALSE])
  })
)
top10_var_gps <- names(sort(
  apply(igt_mean_mat_act, 2, var),
  decreasing = TRUE
))[1:10]
message(paste(
  "Top 10 GPs by var of IGT mean loading:",
  paste(top10_var_gps, collapse = ", ")
))

pdf(
  paste0(figure_path, "GP_Heatmap_Spleen_Top10VarIGTMean.pdf"),
  width = 5,
  height = 5
)
plot_gp_heatmap(
  L_pm_filtered,
  seurat_meta_filtered_spleen,
  gp_vector = top10_var_gps,
  igt_vector = selected_igts,
  igt_stat = "mean",
  top_n_var = NULL,
  floor_threshold = 0,
  color_scheme = c("white", "red"),
  main_title = "Top 10 GPs by Variance of IGT Mean Loading",
  fontsize_row = 8,
  angle_col = 45
)
dev.off()
