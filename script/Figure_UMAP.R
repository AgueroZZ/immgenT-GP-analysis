library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(fastTopics)
library(dplyr)
library(pheatmap)
library(scales)
library(scattermore)

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
figure_path <- "figures/"
source(paste0(code_path, "ROC.R"))
source(paste0(code_path, "filtering_membership.R"))

#####################################################
#####################################################
#####################################################
### Loading Data
#####################################################
#####################################################
load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
# L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered_500.rds"))
# F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered_500.rds"))
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered),]
umap_result <- readRDS(paste0(data_path, "umap_result.rds"))
colnames(umap_result) <- c("UMAP_1", "UMAP_2")
umap_result <- umap_result[rownames(L_pm_filtered),]
df_umap <- as.data.frame(umap_result)




#####################################################
#####################################################
#####################################################
### Defining functions for plotting
#####################################################
#####################################################
plot_loadings_on_umap <- function(umap, loading, factor_num,
                                  size = 0.4,
                                  show_bg = TRUE,
                                  bg_alpha = 0.05,
                                  bg_color = "white",
                                  low_color = "lightgrey",
                                  high_color = "blue",
                                  clip_q = 0.99,
                                  lower_bound = 0.1,     # New: threshold for quantile calculation
                                  top_q = NULL) {

  stopifnot(is.matrix(umap) || is.data.frame(umap))
  df <- as.data.frame(umap)
  stopifnot(all(c("UMAP_1","UMAP_2") %in% colnames(df)))

  # 1. Calculate clipping threshold based only on cells above lower_bound
  # This addresses sparsity by ignoring 'inactive' cells for the color scale calibration
  active_loadings <- loading[loading > lower_bound]

  if(length(active_loadings) > 0) {
    q_hi <- stats::quantile(active_loadings, probs = clip_q, na.rm = TRUE)
  } else {
    # Fallback if no cells exceed the lower_bound
    q_hi <- stats::quantile(loading, probs = clip_q, na.rm = TRUE)
  }

  # 2. Clip values and handle top_q filtering if requested
  loading_clip <- pmin(loading, q_hi)

  if (!is.null(top_q)) {
    thr <- stats::quantile(loading_clip, probs = top_q, na.rm = TRUE)
    loading_plot <- ifelse(loading_clip >= thr, loading_clip, NA_real_)
  } else {
    loading_plot <- loading_clip
  }

  df$loading_plot <- loading_plot

  # 3. Sorting logic: Move higher values to the end of the dataframe
  # This ensures blue points are plotted last (on top), preventing 'grey-out'
  df <- df[order(is.na(df$loading_plot), df$loading_plot, decreasing = FALSE), ]

  # 4. Build Plot
  p <- ggplot() +
    # Background Layer: Structural context of all cells
    {
      if (show_bg)
        scattermore::geom_scattermore(
          data = df,
          mapping = aes(x = UMAP_1, y = UMAP_2),
          pointsize = size,
          color = bg_color,
          alpha = bg_alpha
        )
    } +
    # Foreground Layer: Active factor loadings with prioritized z-order
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
      name = paste0("K", factor_num)
    ) +
    coord_equal() +
    labs(
      title = paste0("Factor ", factor_num),
      x = "UMAP 1", y = "UMAP 2"
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
### Plotting UMAPs
#####################################################
#####################################################
selected_factors <- 1:ncol(L_pm_filtered)
sample_n <- 80000
set.seed(1)
sample_index <- sample(seq_len(nrow(L_pm_filtered)), size = min(sample_n, nrow(L_pm_filtered)))
umap_sub <- umap_result[sample_index, , drop = FALSE]
L_sub    <- L_pm_filtered[sample_index, selected_factors, drop = FALSE]

pdf(paste0(figure_path, "UMAP_factors.pdf"), width = 8, height = 8)
# pdf(paste0(figure_path, "UMAP_factors_500.pdf"), width = 8, height = 8)
for (f in selected_factors) {
  p <-  plot_loadings_on_umap(
    umap    = umap_sub,
    loading = L_sub[, f],
    factor_num = f,
    size = 0.8,
    show_bg = TRUE,
    bg_color = "white", # Lighter background
    low_color = "lightgrey",
    high_color = "blue",
    bg_alpha = 0.6,
    clip_q = 0.9,
    top_q = NULL
  )
  print(p)
}
dev.off()
