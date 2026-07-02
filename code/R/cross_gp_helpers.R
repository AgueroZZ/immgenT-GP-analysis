# Cross-GP comparison heatmap: rows = top signature features (by max |weight|
# across the selected GPs), columns = selected GPs, diverging color scale.
#
# NOT from script/ or code/ -- like plot_gp_signature_volcano() in
# volcano_helpers.R, this is ported from the separate "immgen-signature"
# Shiny app (.../Immgen/webapps/immgen-signature/R/mod_cross_gp.R,
# `.build_heatmap()`), used by script/Figure2.R for panel 2M.

# The app's 4 selectable diverging palettes; "bwr" (Blue-White-Red) is the
# module's default and what the published panel uses.
cross_gp_color_scale <- function(colorscheme = c("bwr", "pbg", "gwr", "rdylbu"), cap = 1, feat_label = "Gene") {
  colorscheme <- match.arg(colorscheme)
  switch(colorscheme,
    bwr = ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-cap, cap), oob = scales::squish,
      name = paste0(feat_label, "\nweight")
    ),
    pbg = ggplot2::scale_fill_gradientn(
      colours = c("#7A0177", "#3A003A", "black", "#6B3000", "#FFD700"),
      values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
      limits = c(-cap, cap), oob = scales::squish,
      name = paste0(feat_label, "\nweight")
    ),
    gwr = ggplot2::scale_fill_gradient2(
      low = "#1A9641", mid = "white", high = "#D7191C",
      midpoint = 0, limits = c(-cap, cap), oob = scales::squish,
      name = paste0(feat_label, "\nweight")
    ),
    rdylbu = ggplot2::scale_fill_distiller(
      palette = "RdYlBu", direction = -1,
      limits = c(-cap, cap), oob = scales::squish,
      name = paste0(feat_label, "\nweight")
    )
  )
}

# mat: features x GPs (already max|w|=1 normalized, e.g. F_pm from
# normalize_maxabs() in volcano_helpers.R). Matches
# mod_cross_gp.R::.build_heatmap() exactly, EXCEPT for the optional
# `pin_top` argument below, which is our own addition (not part of the
# app) for manually forcing specific rows to the top of the heatmap.
plot_cross_gp_heatmap <- function(
  mat,
  gps,
  feat_label = "Gene",
  n_genes = 50,
  direction = c("both", "pos", "neg"),
  threshold = 0.05,
  colorscheme = c("bwr", "pbg", "gwr", "rdylbu"),
  cluster_r = TRUE,
  cluster_c = FALSE,
  pin_top = NULL # character vector of feature names to force to the top
                 # rows (in the given order), after the usual top-n_genes
                 # selection and clustering -- not part of the original app.
) {
  direction <- match.arg(direction)
  colorscheme <- match.arg(colorscheme)

  max_abs <- apply(mat, 1, function(x) max(abs(x), na.rm = TRUE))
  keep <- max_abs >= threshold
  if (direction == "pos") {
    keep <- keep & apply(mat, 1, max, na.rm = TRUE) > 0
  } else if (direction == "neg") {
    keep <- keep & apply(mat, 1, min, na.rm = TRUE) < 0
  }
  mat_f <- mat[keep, , drop = FALSE]

  max_abs_f <- apply(mat_f, 1, function(x) max(abs(x), na.rm = TRUE))
  top_idx <- order(max_abs_f, decreasing = TRUE)[seq_len(min(n_genes, nrow(mat_f)))]
  mat_top <- mat_f[top_idx, , drop = FALSE]

  if (cluster_r && nrow(mat_top) > 2) {
    mat_top <- mat_top[hclust(dist(mat_top))$order, , drop = FALSE]
  }
  if (cluster_c && ncol(mat_top) > 2) {
    mat_top <- mat_top[, hclust(dist(t(mat_top)))$order, drop = FALSE]
  }

  if (!is.null(pin_top)) {
    pin_top <- intersect(pin_top, rownames(mat_top)) # keep only rows actually present
    rest <- setdiff(rownames(mat_top), pin_top)
    # ggplot's discrete y-scale draws the FIRST factor level at the bottom
    # and the LAST at the top, so "pin to top, in this order" means placing
    # them last and reversed (so pin_top[1] ends up the very top row).
    mat_top <- mat_top[c(rest, rev(pin_top)), , drop = FALSE]
  }

  df_long <- as.data.frame(mat_top) |>
    tibble::rownames_to_column("Feature") |>
    tidyr::pivot_longer(cols = -Feature, names_to = "GP", values_to = "Weight") |>
    dplyr::mutate(
      Feature = factor(Feature, levels = rownames(mat_top)),
      GP = factor(GP, levels = colnames(mat_top))
    )

  cap <- 1 # weights are already max|w|=1 normalized

  ggplot2::ggplot(df_long, ggplot2::aes(x = GP, y = Feature, fill = Weight)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    cross_gp_color_scale(colorscheme, cap, feat_label) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(
      title = paste0("Top ", nrow(mat_top), " ", tolower(feat_label), " weights across ", length(gps), " GPs"),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "bold", angle = 0, hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = max(6, min(11, 350 / nrow(mat_top)))),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
}
