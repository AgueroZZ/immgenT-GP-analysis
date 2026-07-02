library(ggplot2)
library(fastTopics)
library(dplyr)
library(ggrepel)
data_path <- "../data/"
code_path <- "../code/"
data_path <- "data/"
code_path <- "code/"
source(paste0(code_path, "ROC.R"))
seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))

cells_flashier <- rownames(flashier_snmf_summary$L_pm)
cells_seurat <- seurat_meta$cellID
L_pm <- flashier_snmf_summary$L_pm[cells_flashier %in% cells_seurat, ]


### Filter cells with high total membership
filter_cells_by_total_membership <- function (L, max_val = 10, numiter = 10) {
  n <- nrow(L)
  rows <- 1:n
  for (iter in 1:numiter) {
    x <- rowSums(L)
    cat(sprintf("%d. Filtered out %d cells.\n",iter,sum(x > max_val)))
    i <- which(x <= max_val)
    L <- L[i,]
    rows <- rows[i]
    d <- apply(L,2,max)
    L <- scale_cols(L,1/d)
  }
  return(rows)
}
scale_cols <- function (A, b) t(t(A) * b)
D <- diag(1 / apply(L_pm, 2, function(x) max(x)))
L <- L_pm %*% D
cells <- filter_cells_by_total_membership(L,numiter = 12)
seurat_meta_filtered <- seurat_meta[cells,]
L_pm_filtered <- L_pm[cells,]
d <- apply(L_pm_filtered,2,max)
L_pm_filtered <- scale_cols(L_pm_filtered,1/d)
D <- (1/d) * D


# Plot all UMAPs
umap_result <- readRDS(paste0(data_path, "umap_result.rds"))
colnames(umap_result) <- c("UMAP_1", "UMAP_2")
umap_result <- umap_result[rownames(L_pm_filtered),]

df_umap <- as.data.frame(umap_result)
df_umap$CellType <- seurat_meta_filtered$annotation_level1
stopifnot(nrow(df_umap) == length(seurat_meta_filtered$annotation_level1))

plot_loadings_on_umap_v2 <- function(umap, loading, factor_num,
                                     size = 0.25,
                                     show_bg = TRUE,
                                     bg_alpha = 0.05,
                                     bg_color = "grey90",
                                     clip_q = 0.99,
                                     top_q = NULL,
                                     viridis_option = "viridis",
                                     viridis_direction = 1) {
  stopifnot(is.matrix(umap) || is.data.frame(umap))
  df <- as.data.frame(umap)
  stopifnot(all(c("UMAP_1","UMAP_2") %in% colnames(df)))
  stopifnot(nrow(df) == length(loading))

  q_hi <- stats::quantile(loading, probs = clip_q, na.rm = TRUE)
  loading_clip <- pmin(loading, q_hi)

  if (!is.null(top_q)) {
    stopifnot(top_q > 0 && top_q < 1)
    thr <- stats::quantile(loading_clip, probs = top_q, na.rm = TRUE)
    loading_plot <- ifelse(loading_clip >= thr, loading_clip, NA_real_)
  } else {
    loading_plot <- loading_clip
  }

  df$loading_plot <- loading_plot

  p <- ggplot2::ggplot() +
    {
      if (show_bg) ggplot2::geom_point(
        data = df, ggplot2::aes(UMAP_1, UMAP_2),
        color = bg_color, size = size, alpha = bg_alpha
      )
    } +
    ggplot2::geom_point(
      data = df, ggplot2::aes(UMAP_1, UMAP_2, color = loading_plot),
      size = size, alpha = 0.9, na.rm = TRUE
    ) +
    ggplot2::scale_color_viridis_c(
      option = viridis_option,
      direction = viridis_direction,
      limits = c(0, q_hi),
      oob = scales::squish,
      breaks = scales::pretty_breaks(4),
      name = paste0("K", factor_num, " loading")
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0("Factor ", factor_num),
      x = "", y = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )

  return(p)
}

selected_factors <- 1:ncol(L_pm_filtered)

sample_n <- 80000

set.seed(1)
sample_index <- sample(seq_len(nrow(L_pm_filtered)), size = min(sample_n, nrow(L_pm_filtered)))
umap_sub <- umap_result[sample_index, , drop = FALSE]
L_sub    <- L_pm_filtered[sample_index, selected_factors, drop = FALSE]

for (f in selected_factors) {
  p <-  plot_loadings_on_umap_v2(
    umap    = umap_sub,
    loading = L_sub[, f],
    factor_num = f,
    size = 0.4,
    show_bg = TRUE,
    bg_alpha = 0.06,
    clip_q = 1,
    top_q = NULL,
    viridis_option = "viridis",
    viridis_direction = 1
  )

  ggsave(
    filename = paste0("figures/UMAP_fixed_color/UMAP_F", f, ".png"),
    plot     = p,
    width    = 8,
    height   = 8,
    dpi      = 300
  )
}


for (f in selected_factors) {
  p <-  plot_loadings_on_umap_v2(
    umap    = umap_sub,
    loading = L_sub[, f],
    factor_num = f,
    size = 0.4,
    show_bg = TRUE,
    bg_alpha = 0.06,
    clip_q = 0.99,
    top_q = NULL,
    viridis_option = "viridis",
    viridis_direction = 1
  )

  ggsave(
    filename = paste0("figures/UMAP/UMAP_F", f, ".png"),
    plot     = p,
    width    = 8,
    height   = 8,
    dpi      = 300
  )
}

