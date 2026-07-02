# Protein-gate vs. GP-loading comparison on an MDE embedding. Ported
# from gated_protein_loading_plot.R; shared by script/Figure6.R
# (panels c-f, the 4 main-figure GPs) and FigureS6.R (the full gallery).

gp_label <- function(x) sub("^K", "GP", x)

# Density-colored highlight of a subset of cells on a 2D embedding, with all
# other cells shown as a grey background layer.
MyDimPlotHighlightDensity_df <- function(
  emb,
  highlight,
  split = NULL,
  dim_names = c("Dim1", "Dim2"),
  raster = TRUE,
  highlight_size = 0.5,
  highlight_alpha = 0.5,
  highlight_pointsize = 0L,
  base_pixels = c(512, 512),
  highlight_pixels = c(512, 512),
  cols = rev(rainbow(10, end = 4 / 6)),
  nbin = 500
) {
  requireNamespace("scattermore", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)

  if (missing(emb) || missing(highlight)) {
    stop("Please provide emb (n x 2) and highlight.")
  }
  emb <- as.data.frame(emb)
  if (ncol(emb) < 2) stop("emb must have at least 2 columns (dim1, dim2).")
  emb <- emb[, 1:2, drop = FALSE]

  n <- nrow(emb)
  if (length(highlight) != n) stop("highlight must have length nrow(emb).")
  if (!is.null(split) && length(split) != n) stop("split must have length nrow(emb).")
  highlight <- as.logical(highlight)

  df <- data.frame(feature1 = as.numeric(emb[[1]]), feature2 = as.numeric(emb[[2]]), highlight = highlight, stringsAsFactors = FALSE)
  if (!is.null(split)) df$split <- split

  keep <- !is.na(df$feature1) & !is.na(df$feature2)
  if (!is.null(split)) keep <- keep & !is.na(df$split)
  df <- df[keep, , drop = FALSE]
  df2 <- df[!is.na(df$highlight) & df$highlight, , drop = FALSE]

  p1 <- ggplot2::ggplot(df) +
    scattermore::geom_scattermore(ggplot2::aes(feature1, feature2), color = "grey", pixels = base_pixels)

  if (nrow(df2) > 0) {
    # densCols needs enough unique points to form bin breaks; fall back to a
    # solid colour (top of the ramp) when there are too few highlighted cells.
    df2$density_col <- tryCatch(
      grDevices::densCols(df2$feature1, df2$feature2, colramp = grDevices::colorRampPalette(cols), nbin = nbin),
      error = function(e) rep(tail(cols, 1), nrow(df2))
    )
    if (isTRUE(raster)) {
      p2 <- scattermore::geom_scattermore(data = df2, ggplot2::aes(feature1, feature2, color = density_col), pixels = highlight_pixels, pointsize = highlight_pointsize)
    } else {
      p2 <- ggplot2::geom_point(data = df2, ggplot2::aes(feature1, feature2, color = density_col), size = highlight_size, alpha = highlight_alpha)
    }
  } else {
    p2 <- NULL
  }

  p <- p1 + p2 +
    ggplot2::xlab(dim_names[1]) + ggplot2::ylab(dim_names[2]) +
    ggplot2::scale_colour_identity() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 15), axis.text.y = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 10),
      axis.title.x = ggplot2::element_text(size = 20), axis.title.y = ggplot2::element_text(size = 20),
      legend.title = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()
    )
  if (!is.null(split)) p <- p + ggplot2::facet_wrap(~split)
  p
}

# For one GP, builds a protein gate from its curated marker signature
# (positive markers above threshold, negative markers at/below), and
# separately gates the same-sized top-loading cells; returns the two panels
# side by side (protein gate | GP-loading gate) on the same embedding.
plot_gated_gp_vs_protein <- function(
  gp_name,
  df_markers,
  protein_mat,
  loading_mat,
  mde_emb,
  threshold_df,
  selected_proteins = NULL,
  exclude_cells = NULL,
  loading_q = 0.9,
  min_cells = 2,
  missing_threshold_action = "median",
  min_pointsize = 0L,
  save_path = NULL
) {
  if (!is.null(exclude_cells)) {
    cells_to_keep <- setdiff(rownames(loading_mat), exclude_cells)
    loading_mat <- loading_mat[cells_to_keep, , drop = FALSE]
    common <- intersect(cells_to_keep, rownames(protein_mat))
    protein_mat <- protein_mat[common, , drop = FALSE]
    common_mde <- intersect(common, rownames(mde_emb))
    mde_emb <- mde_emb[common_mde, , drop = FALSE]
  }

  row_idx <- which(rownames(df_markers) == gp_label(gp_name))
  if (length(row_idx) == 0) stop("GP name not found in df_markers.")

  pos_markers <- strsplit(df_markers[row_idx, "Positive"], ", ")[[1]]
  neg_markers <- strsplit(df_markers[row_idx, "Negative"], ", ")[[1]]
  if (!is.null(selected_proteins)) {
    pos_markers <- intersect(pos_markers, selected_proteins)
    neg_markers <- intersect(neg_markers, selected_proteins)
  }
  pos_markers <- pos_markers[pos_markers != ""]
  neg_markers <- neg_markers[neg_markers != ""]

  is_protein_gated <- rep(TRUE, nrow(protein_mat))
  names(is_protein_gated) <- rownames(protein_mat)

  apply_gate <- function(marker_list, is_positive = TRUE) {
    for (p in marker_list) {
      if (p %in% colnames(protein_mat)) {
        if (p %in% threshold_df$Protein) {
          thresh <- threshold_df$Threshold[threshold_df$Protein == p]
          if (is_positive) {
            is_protein_gated <<- is_protein_gated & (protein_mat[, p] > thresh)
          } else {
            is_protein_gated <<- is_protein_gated & (protein_mat[, p] <= thresh)
          }
        } else {
          if (missing_threshold_action == "skip") {
            message(sprintf("[%s] Warning: %s threshold missing. Skipping this marker.", gp_name, p))
          } else {
            thresh <- median(protein_mat[, p], na.rm = TRUE)
            message(sprintf("[%s] Warning: %s threshold missing. Using median (%.2f).", gp_name, p, thresh))
            if (is_positive) {
              is_protein_gated <<- is_protein_gated & (protein_mat[, p] > thresh)
            } else {
              is_protein_gated <<- is_protein_gated & (protein_mat[, p] <= thresh)
            }
          }
        }
      }
    }
  }
  apply_gate(pos_markers, is_positive = TRUE)
  apply_gate(neg_markers, is_positive = FALSE)
  n_prot <- sum(is_protein_gated)

  loadings <- loading_mat[, gp_name]
  if (is.null(loading_q)) {
    if (n_prot <= 1) {
      loading_q_val <- 0.999
    } else {
      loading_q_val <- 1 - (n_prot / length(loadings))
      loading_q_val <- max(0, min(0.9999, loading_q_val))
    }
    q_label <- paste0("Matched n=", n_prot)
  } else {
    loading_q_val <- loading_q
    q_label <- paste0("Top ", round((1 - loading_q_val) * 100), "%")
  }
  loading_cutoff <- quantile(loadings, loading_q_val)
  is_loading_gated <- loadings >= loading_cutoff
  n_load <- sum(is_loading_gated)

  message(sprintf("[%s] Final Gate: Protein=%d, Loading=%d", gp_name, n_prot, n_load))

  common_cells <- intersect(rownames(loading_mat), rownames(mde_emb))
  emb_subset <- mde_emb[common_cells, ]
  pos_subtitle <- if (length(pos_markers) > 0) paste0(paste(pos_markers, collapse = "+ "), "+") else "None+"
  neg_subtitle <- if (length(neg_markers) > 0) paste0(paste(neg_markers, collapse = "- "), "-") else "None-"

  render_plot <- function(highlight_vec, title_text, cell_count) {
    if (cell_count < min_cells) {
      label_text <- if (cell_count == 0) "Zero cells in gate" else paste0("Too few cells (n=", cell_count, ")")
      return(
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5, label = label_text, size = 5, fontface = "italic") +
          ggplot2::ggtitle(title_text) +
          ggplot2::theme_void() +
          ggplot2::theme(plot.title = ggplot2::element_text(size = 9, face = "bold"))
      )
    }
    hl_pointsize <- if (cell_count < 100) 8L else if (cell_count < 1000) 3L else 0L
    hl_pointsize <- max(hl_pointsize, min_pointsize)

    MyDimPlotHighlightDensity_df(
      emb = emb_subset, highlight = highlight_vec[common_cells], dim_names = c("MDE1", "MDE2"),
      cols = grDevices::colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(10),
      highlight_pointsize = hl_pointsize
    ) +
      ggplot2::ggtitle(title_text) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 9, face = "bold"),
        axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank()
      )
  }

  gp_disp <- gp_label(gp_name)
  p1 <- render_plot(is_protein_gated, paste0(gp_disp, " (", q_label, ")\n", pos_subtitle, "\n", neg_subtitle), n_prot)
  p2 <- render_plot(is_loading_gated, "", n_load)
  combined_plot <- p1 + p2 + patchwork::plot_layout(ncol = 2)

  if (!is.null(save_path)) {
    ggplot2::ggsave(filename = save_path, plot = combined_plot, width = 12, height = 6)
  }
  combined_plot
}
