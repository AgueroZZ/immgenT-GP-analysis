library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(fastTopics)
library(dplyr)
library(pheatmap)
library(scales)
library(scattermore)
library(tidyr)
library(ggalluvial)
library(knitr)
library(kableExtra)
library(ggbeeswarm)


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
figure_path <- "figures/Figure6_CITEseq/"


#####################################################
#####################################################
####################################################
##### Defining additional functions
####################################################
#####################################################
plot_scatterplot <- function(factor_index, n_higlight = 10) {
  f_values <- citeseq_F[, factor_index]
  f_values <- f_values / max(abs(f_values))
  gene_names <- rownames(citeseq_F)
  df <- data.frame(
    Gene = gene_names,
    FactorValue = f_values,
    MeanLogExpr = mean_shifted_log_expr[gene_names]
  )
  p <- ggplot(df, aes(x = FactorValue, y = MeanLogExpr))
  p <- p + geom_point(alpha = 0.5)
  p <- p +
    labs(
      title = colnames(citeseq_F)[factor_index],
      x = "Increase in Shifted Log Expression",
      y = "Mean Shifted Log Expression"
    ) +
    theme_minimal(base_size = 13)
  top_genes <- df %>%
    slice_max(order_by = abs(FactorValue), n = n_higlight)
  p <- p +
    geom_point(
      data = top_genes,
      aes(x = FactorValue, y = MeanLogExpr),
      color = "skyblue"
    )
  p <- p +
    geom_text_repel(
      data = top_genes,
      aes(label = Gene),
      color = "skyblue",
      box.padding = 0.35,
      point.padding = 0.15,
      max.time = 1,
      max.iter = 5e3
    )
  return(p)
}
marker_list_to_df <- function(gp_marker_sign_list, collapse = ", ") {
  stopifnot(is.list(gp_marker_sign_list), length(gp_marker_sign_list) > 0)

  sets <- names(gp_marker_sign_list)
  if (is.null(sets)) {
    stop("gp_marker_sign_list must have names, e.g. 'GP68'.")
  }

  df <- data.frame(
    Set = sets,
    Positive = vapply(
      gp_marker_sign_list,
      function(x) {
        pos <- x$pos
        if (length(pos) == 0) {
          return("")
        }
        paste(pos, collapse = collapse)
      },
      character(1)
    ),
    Negative = vapply(
      gp_marker_sign_list,
      function(x) {
        neg <- x$neg
        if (length(neg) == 0) {
          return("")
        }
        paste(neg, collapse = collapse)
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )

  df
}
MyFeatureScatter_df <- function(
  x,
  y,
  highlight,
  split = NULL,
  feature1 = "feature1",
  feature2 = "feature2",
  raster = TRUE,
  cols = rev(rainbow(10, end = 4 / 6)),
  highlight_size = 1,
  highlight_alpha = 1,
  base_pixels = c(512, 512),
  highlight_pixels = c(216, 216),
  nbin = 500
) {
  requireNamespace("scattermore", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)

  # ---- checks ----
  if (missing(x) || missing(y) || missing(highlight)) {
    stop("Please provide x, y, and highlight.")
  }
  if (length(x) != length(y) || length(x) != length(highlight)) {
    stop("x, y, and highlight must have the same length.")
  }
  if (!is.null(split) && length(split) != length(x)) {
    stop("split must have the same length as x/y.")
  }

  # coerce highlight to logical safely
  highlight <- as.logical(highlight)

  df <- data.frame(
    feature1 = as.numeric(x),
    feature2 = as.numeric(y),
    highlight = highlight,
    stringsAsFactors = FALSE
  )
  if (!is.null(split)) {
    df$split <- split
  }

  # drop rows with NA in x/y (and split if present)
  keep <- !is.na(df$feature1) & !is.na(df$feature2)
  if (!is.null(split)) {
    keep <- keep & !is.na(df$split)
  }
  df <- df[keep, , drop = FALSE]

  df2 <- df[!is.na(df$highlight) & df$highlight, , drop = FALSE]

  # If nothing to highlight, still return the base plot
  p1 <- ggplot2::ggplot(df) +
    scattermore::geom_scattermore(
      ggplot2::aes(feature1, feature2),
      color = "grey",
      pixels = base_pixels
    )

  if (nrow(df2) > 0) {
    df2$density_col <- grDevices::densCols(
      df2$feature1,
      df2$feature2,
      colramp = grDevices::colorRampPalette(cols),
      nbin = nbin
    )

    if (isTRUE(raster)) {
      p2 <- scattermore::geom_scattermore(
        data = df2,
        ggplot2::aes(feature1, feature2, color = density_col),
        pointsize = highlight_size,
        pixels = highlight_pixels
      )
    } else {
      p2 <- ggplot2::geom_point(
        data = df2,
        ggplot2::aes(feature1, feature2, color = density_col),
        size = highlight_size,
        alpha = highlight_alpha
      )
    }
  } else {
    p2 <- NULL
  }

  p <- p1 +
    p2 +
    ggplot2::xlab(feature1) +
    ggplot2::ylab(feature2) +
    ggplot2::scale_colour_identity() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 10),
      axis.title.x = ggplot2::element_text(size = 20),
      axis.title.y = ggplot2::element_text(size = 20),
      legend.title = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(split)) {
    p <- p + ggplot2::facet_wrap(~split)
  }

  return(p)
}
MyDimPlotHighlightDensity_df <- function(
  emb, # n x 2 matrix/data.frame: columns are dim1, dim2
  highlight, # length n logical (or coercible): which cells to highlight
  split = NULL, # optional length n vector for facet
  dim_names = c("Dim1", "Dim2"),
  raster = TRUE,
  highlight_size = 0.5,
  highlight_alpha = 0.5,
  base_pixels = c(512, 512),
  highlight_pixels = c(512, 512),
  cols = rev(rainbow(10, end = 4 / 6)),
  nbin = 500
) {
  requireNamespace("scattermore", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)

  # ---- checks ----
  if (missing(emb) || missing(highlight)) {
    stop("Please provide emb (n x 2) and highlight.")
  }

  emb <- as.data.frame(emb)
  if (ncol(emb) < 2) {
    stop("emb must have at least 2 columns (dim1, dim2).")
  }
  emb <- emb[, 1:2, drop = FALSE]

  n <- nrow(emb)
  if (length(highlight) != n) {
    stop("highlight must have length nrow(emb).")
  }
  if (!is.null(split) && length(split) != n) {
    stop("split must have length nrow(emb).")
  }

  highlight <- as.logical(highlight)

  df <- data.frame(
    feature1 = as.numeric(emb[[1]]),
    feature2 = as.numeric(emb[[2]]),
    highlight = highlight,
    stringsAsFactors = FALSE
  )
  if (!is.null(split)) {
    df$split <- split
  }

  # drop NA coords (and split if present)
  keep <- !is.na(df$feature1) & !is.na(df$feature2)
  if (!is.null(split)) {
    keep <- keep & !is.na(df$split)
  }
  df <- df[keep, , drop = FALSE]

  df2 <- df[!is.na(df$highlight) & df$highlight, , drop = FALSE]

  # ---- base plot (all cells in grey) ----
  p1 <- ggplot2::ggplot(df) +
    scattermore::geom_scattermore(
      ggplot2::aes(feature1, feature2),
      color = "grey",
      pixels = base_pixels
    )

  # ---- highlight layer (density-colored) ----
  if (nrow(df2) > 0) {
    df2$density_col <- grDevices::densCols(
      df2$feature1,
      df2$feature2,
      colramp = grDevices::colorRampPalette(cols),
      nbin = nbin
    )

    if (isTRUE(raster)) {
      p2 <- scattermore::geom_scattermore(
        data = df2,
        ggplot2::aes(feature1, feature2, color = density_col),
        pixels = highlight_pixels
      )
    } else {
      p2 <- ggplot2::geom_point(
        data = df2,
        ggplot2::aes(feature1, feature2, color = density_col),
        size = highlight_size,
        alpha = highlight_alpha
      )
    }
  } else {
    p2 <- NULL
  }

  p <- p1 +
    p2 +
    ggplot2::xlab(dim_names[1]) +
    ggplot2::ylab(dim_names[2]) +
    ggplot2::scale_colour_identity() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 10),
      axis.title.x = ggplot2::element_text(size = 20),
      axis.title.y = ggplot2::element_text(size = 20),
      legend.title = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(split)) {
    p <- p + ggplot2::facet_wrap(~split)
  }

  return(p)
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
  background_color = "black",
  x_limits = c(-0.5, 0.5),
  y_limits = c(-0.5, 0.5),
  background_alpha = 0.5,
  xlab = "Difference in Mean Loading",
  ylab = "Difference in Mean Loading",
  title = "Comparison of Specific GP Loadings"
) {
  # Filter the dataframe for only the specific GPs you care about
  highlight_df <- df %>%
    filter({{ label_var }} %in% target_gps) %>%
    mutate(.label_display = as.character({{ label_var }}))

  # Check if any GPs were found to prevent ggplot errors
  if (nrow(highlight_df) == 0) {
    warning("None of the target_gps were found in the dataframe.")
  }

  # Build the scatter plot
  p <- ggplot(df, aes(x = {{ x_var }}, y = {{ y_var }})) +

    # 1. Plot background points (unselected GPs) in muted grey
    geom_point(color = background_color, alpha = background_alpha) +

    # 2. Add reference lines matching your original plot
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +

    # 3. Plot the targeted GPs on top, mapping color inside aes()
    geom_point(data = highlight_df, aes(color = {{ label_var }}), size = 2) +

    # 4. Label only the targeted GPs, mapping color inside aes()
    ggrepel::geom_text_repel(
      data = highlight_df,
      aes(label = .label_display, color = {{ label_var }}),
      max.overlaps = Inf,
      size = 3.5,
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
get_klrg1_split <- function(cell_type_label, meta, protein_data, threshold) {
  cells <- meta$cellID[meta$annotation_level1 == cell_type_label]
  cells <- intersect(cells, rownames(protein_data))

  pos <- cells[protein_data[cells, "KLRG1"] >= threshold]
  neg <- cells[protein_data[cells, "KLRG1"] < threshold]

  return(list(pos = pos, neg = neg))
}
run_checked_dge <- function(group_list, F_mat, L_mat, label) {
  if (length(group_list$pos) < 3 || length(group_list$neg) < 3) {
    stop(paste("Insufficient data:", label))
  }

  dge_res <- FlashierDGE_corrected(
    F1 = F_mat,
    L1 = L_mat,
    group1 = group_list$pos,
    group2 = group_list$neg,
    title_plot = paste("KLRG1+ vs - in", label)
  )

  df <- dge_res$diff_factors
  if (!"SYMBOL" %in% colnames(df)) {
    df$SYMBOL <- rownames(df)
  }
  return(df)
}
my_swarm_plot <- function(
  input_matrix,
  title = "Predictive Power across all GPs",
  y_label = "Statistic Value",
  label_threshold_quantile = 0.99,
  top_hit = 1,
  global_threshold = NULL,
  manual_colors = c("steelblue", "darkorange", "forestgreen", "purple"),
  n_colors = 2,
  use_beeswarm = TRUE, # TRUE = beeswarm offsets, FALSE = jitter
  beeswarm_cex = 1.5, # spacing between beeswarm points
  jitter_width = 0.2, # spread width for jitter fallback
  label_direction = "both",
  label_nudge_y = 0.005,
  label_force = 1.5,
  label_force_pull = 1.0
) {
  # ── 1. Prepare data ──────────────────────────────────────────────────────────
  unique_features <- rownames(input_matrix)
  n_features <- length(unique_features)

  plot_df <- as.data.frame(input_matrix) %>%
    mutate(Feature = rownames(.)) %>%
    pivot_longer(cols = -Feature, names_to = "GP", values_to = "Value")

  # ── 2. Color groups ───────────────────────────────────────────────────────────
  color_map <- manual_colors[seq_len(n_colors)]
  feature_info <- data.frame(
    Feature = unique_features,
    ColorGroup = as.factor(rep(seq_len(n_colors), length.out = n_features)),
    ActualColor = color_map[rep(seq_len(n_colors), length.out = n_features)]
  )

  plot_df <- plot_df %>%
    left_join(feature_info, by = "Feature") %>%
    mutate(Feature = factor(Feature, levels = unique_features))

  # ── 3. Threshold flags ────────────────────────────────────────────────────────
  plot_df$is_high <- if (!is.null(global_threshold)) {
    plot_df$Value >= global_threshold
  } else {
    FALSE
  }

  # ── 4. Pre-compute displaced x positions ─────────────────────────────────────
  # Both points AND labels will use x_jit directly via geom_point() /
  # geom_text_repel(), so they are guaranteed to share identical coordinates.
  # This avoids the mismatch that arises when geom_jitter uses its own internal
  # RNG independently of any manually computed offsets.
  plot_df <- plot_df %>%
    mutate(x_int = as.integer(Feature))

  if (use_beeswarm) {
    # beeswarm::swarmx() computes the same deterministic offsets that
    # geom_beeswarm uses internally, called once per feature group.
    plot_df <- plot_df %>%
      group_by(Feature) %>%
      mutate(
        x_jit = beeswarm::swarmx(
          x = rep(x_int[1], n()),
          y = Value,
          cex = beeswarm_cex,
          side = 0L # 0 = spread on both sides
        )$x
      ) %>%
      ungroup()
  } else {
    # Fixed-seed uniform jitter; stored in x_jit so points and labels both
    # read from the same column — no independent RNG calls.
    set.seed(42)
    plot_df <- plot_df %>%
      mutate(x_jit = x_int + runif(n(), -jitter_width, jitter_width))
  }

  # ── 5. Label candidates ───────────────────────────────────────────────────────
  top_hits <- plot_df %>%
    group_by(Feature) %>%
    slice_max(Value, n = top_hit, with_ties = FALSE) %>%
    ungroup()

  if (!is.null(label_threshold_quantile)) {
    q_val <- quantile(plot_df$Value, label_threshold_quantile, na.rm = TRUE)
    top_hits <- top_hits %>% filter(Value >= q_val)
  }

  if (!is.null(global_threshold)) {
    global_hits <- plot_df %>%
      filter(Value >= global_threshold) %>%
      group_by(Feature) %>%
      slice_max(Value, n = top_hit, with_ties = FALSE) %>%
      ungroup()
    top_hits <- bind_rows(top_hits, global_hits) %>% distinct()
  }

  # ── 6. Build plot ─────────────────────────────────────────────────────────────
  # Both layers use x_jit as the x aesthetic — no position objects needed.
  # The x axis is mapped to the original Feature factor so tick labels and
  # grid lines stay correct; x_jit values are numeric indices with small
  # offsets so they land on the right tick positions automatically.
  p <- ggplot(plot_df, aes(color = ColorGroup)) +
    # Point layer — x_jit is the actual rendered position
    geom_point(
      aes(
        x = x_jit,
        y = Value,
        alpha = ifelse(is_high, 1.0, 0.4),
        size = ifelse(is_high, 1.2, 0.7)
      )
    ) +
    # Label layer — starts from identical x_jit, so segment connects correctly
    geom_text_repel(
      data = top_hits,
      aes(x = x_jit, y = Value, label = GP),
      size = 3,
      direction = label_direction,
      nudge_y = label_nudge_y,
      hjust = 0,
      point.padding = 0.2,
      box.padding = 0.3,
      min.segment.length = 0,
      segment.size = 0.3,
      segment.alpha = 0.4,
      segment.color = "grey50",
      force = label_force,
      force_pull = label_force_pull,
      max.overlaps = Inf,
      seed = 42,
      show.legend = FALSE
    ) +
    # Restore proper categorical x axis using the original Feature factor
    scale_x_continuous(
      breaks = seq_len(n_features),
      labels = unique_features
    ) +
    scale_color_manual(values = color_map) +
    scale_size_identity() +
    scale_alpha_identity() +
    labs(title = title, x = "Protein Marker", y = y_label) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 7,
        color = feature_info$ActualColor
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  if (!is.null(global_threshold)) {
    p <- p +
      geom_hline(
        yintercept = global_threshold,
        linetype = "dotted",
        color = "red",
        alpha = 0.5
      )
  }

  return(p)
}
plot_gp_vs_protein <- function(
  loading_matrix,
  protein_matrix,
  gp,
  protein,
  n_cells = 1000,
  seed = 42,
  nonzero_protein = FALSE,
  nonzero_loading = FALSE,
  add_smoother = FALSE,
  smoother_method = "loess",
  smoother_se = TRUE,
  add_quantile_smoother = TRUE,
  n_quantile_bins = 10,
  quantile_smoother_color = "firebrick",
  quantile_smoother_size = 1.0,
  color_by_density = TRUE,
  point_alpha = 0.5,
  point_size = 0.8,
  point_color = "steelblue",
  add_rug = FALSE,
  xlim = NULL,
  ylim = NULL,
  vline_protein_threshold = NULL,
  vline_color = "dodgerblue",
  vline_label = NULL,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  highlight_cells = NULL,
  highlight_size = 2.0
) {
  # ── 1. Input checks ───────────────────────────────────────────────────────────
  if (!gp %in% colnames(loading_matrix)) {
    stop("GP '", gp, "' not found in loading_matrix columns.")
  }
  if (!protein %in% colnames(protein_matrix)) {
    stop("Protein '", protein, "' not found in protein_matrix columns.")
  }

  # ── 2. Shared cells ───────────────────────────────────────────────────────────
  shared_cells <- intersect(rownames(loading_matrix), rownames(protein_matrix))
  if (length(shared_cells) == 0) {
    stop("No shared row names between the two matrices.")
  }

  # ── 3. Build full data frame ──────────────────────────────────────────────────
  full_df <- data.frame(
    cell = shared_cells,
    x = as.numeric(protein_matrix[shared_cells, protein]), # protein on x
    y = as.numeric(loading_matrix[shared_cells, gp]), # loading on y
    row.names = NULL
  )

  # ── 4. Optional filters (applied before sampling) ─────────────────────────────
  if (nonzero_protein) {
    n_before <- nrow(full_df)
    full_df <- full_df[full_df$x != 0, ]
    message(sprintf(
      "nonzero_protein: removed %d cells with zero %s expression (%d remaining).",
      n_before - nrow(full_df),
      protein,
      nrow(full_df)
    ))
    if (nrow(full_df) == 0) {
      stop("No cells remain after nonzero_protein filter.")
    }
  }

  if (nonzero_loading) {
    n_before <- nrow(full_df)
    full_df <- full_df[full_df$y != 0, ]
    message(sprintf(
      "nonzero_loading: removed %d cells with zero %s loading (%d remaining).",
      n_before - nrow(full_df),
      gp,
      nrow(full_df)
    ))
    if (nrow(full_df) == 0) {
      stop("No cells remain after nonzero_loading filter.")
    }
  }

  # ── 5. Compute quantile smoother on full (post-filter) data ───────────────────
  # Using the full dataset (before subsampling) gives stable bin means even
  # when n_cells is small. Each bin contains an equal number of cells.
  quantile_df <- NULL
  if (add_quantile_smoother) {
    quantile_df <- full_df %>%
      mutate(
        bin = cut(
          x,
          breaks = quantile(
            x,
            probs = seq(0, 1, length.out = n_quantile_bins + 1),
            na.rm = TRUE
          ),
          include.lowest = TRUE,
          labels = FALSE # integer bin index
        )
      ) %>%
      group_by(bin) %>%
      summarise(
        x_mean = mean(x, na.rm = TRUE), # mean protein expression in bin
        y_mean = mean(y, na.rm = TRUE), # mean GP loading in bin
        n = n(),
        .groups = "drop"
      ) %>%
      arrange(x_mean)
  }

  # ── 6. Subsample for scatter points ───────────────────────────────────────────
  plot_df <- full_df
  if (!is.null(n_cells) && n_cells < nrow(plot_df)) {
    set.seed(seed)
    plot_df <- plot_df[sample(nrow(plot_df), size = n_cells, replace = FALSE), ]
  }

  # ── 7. Optional density colouring (on subsample) ─────────────────────────────
  if (color_by_density) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      warning(
        "Package 'MASS' needed for color_by_density. Falling back to fixed colour."
      )
      color_by_density <- FALSE
    } else {
      dens <- MASS::kde2d(plot_df$x, plot_df$y, n = 200)
      ix <- pmax(1L, pmin(findInterval(plot_df$x, dens$x), nrow(dens$z)))
      iy <- pmax(1L, pmin(findInterval(plot_df$y, dens$y), ncol(dens$z)))
      plot_df$density <- dens$z[cbind(ix, iy)]
      plot_df <- plot_df[order(plot_df$density), ]
    }
  }

  # ── 8. Highlight flag ─────────────────────────────────────────────────────────
  plot_df$highlighted <- plot_df$cell %in% highlight_cells

  # ── 9. Auto labels ────────────────────────────────────────────────────────────
  filter_tags <- c(
    if (nonzero_protein) "non-zero protein" else NULL,
    if (nonzero_loading) "non-zero loading" else NULL
  )
  filter_str <- if (length(filter_tags) > 0) {
    paste0(" [", paste(filter_tags, collapse = ", "), "]")
  } else {
    ""
  }

  x_lab <- if (!is.null(x_label)) x_label else paste0(protein, " expression")
  y_lab <- if (!is.null(y_label)) {
    y_label
  } else {
    paste0("Loading in ", gp)
  }
  ttl <- if (!is.null(title)) {
    title
  } else {
    paste0(
      gp,
      " loading vs ",
      protein,
      " (n = ",
      formatC(nrow(plot_df), format = "d", big.mark = ","),
      " shown / ",
      formatC(nrow(full_df), format = "d", big.mark = ","),
      " total)",
      filter_str
    )
  }

  # ── 10. Correlation annotation (on full post-filter data) ─────────────────────
  ct <- cor.test(full_df$x, full_df$y, method = "pearson")
  cor_label <- paste0(
    "r = ",
    round(ct$estimate, 3),
    ", p ",
    ifelse(ct$p.value < 0.001, "< 0.001", paste0("= ", round(ct$p.value, 3)))
  )

  # ── 11. Base plot ─────────────────────────────────────────────────────────────
  p <- ggplot(plot_df, aes(x = x, y = y))

  if (color_by_density) {
    p <- p +
      geom_point(
        data = subset(plot_df, !highlighted),
        aes(color = density),
        size = point_size,
        alpha = point_alpha
      ) +
      scale_color_viridis_c(
        option = "magma",
        name = "Density",
        guide = guide_colorbar(barwidth = 0.5, barheight = 5)
      )
  } else {
    p <- p +
      geom_point(
        data = subset(plot_df, !highlighted),
        color = point_color,
        size = point_size,
        alpha = point_alpha
      )
  }

  if (!is.null(highlight_cells) && any(plot_df$highlighted)) {
    p <- p +
      geom_point(
        data = subset(plot_df, highlighted),
        color = "red",
        size = highlight_size,
        alpha = 0.9
      )
  }

  # ── 12. Classical smoother ────────────────────────────────────────────────────
  if (add_smoother) {
    p <- p +
      geom_smooth(
        method = smoother_method,
        formula = y ~ x,
        se = smoother_se,
        color = "royalblue",
        linewidth = 0.8,
        alpha = 0.15
      )
  }

  # ── 13. Quantile smoother ─────────────────────────────────────────────────────
  # Line connecting mean loading within each equal-frequency protein bin,
  # with a dot at each bin centre. Drawn after points so it sits on top.
  if (add_quantile_smoother && !is.null(quantile_df)) {
    p <- p +
      geom_line(
        data = quantile_df,
        aes(x = x_mean, y = y_mean),
        color = quantile_smoother_color,
        linewidth = quantile_smoother_size,
        inherit.aes = FALSE
      ) +
      geom_point(
        data = quantile_df,
        aes(x = x_mean, y = y_mean),
        color = quantile_smoother_color,
        size = quantile_smoother_size * 2,
        shape = 21,
        fill = "white",
        stroke = 1,
        inherit.aes = FALSE
      )
  }

  # ── 14. Optional rug ──────────────────────────────────────────────────────────
  if (add_rug) {
    p <- p +
      geom_rug(
        sides = "bl",
        alpha = 0.15,
        length = unit(0.015, "npc")
      )
  }

  # ── 15. Vertical threshold line ───────────────────────────────────────────────
  if (!is.null(vline_protein_threshold)) {
    p <- p +
      geom_vline(
        xintercept = vline_protein_threshold,
        linetype = "dashed",
        color = vline_color,
        linewidth = 0.6,
        alpha = 0.8
      )
    if (!is.null(vline_label)) {
      p <- p +
        annotate(
          "text",
          x = vline_protein_threshold,
          y = Inf,
          label = vline_label,
          hjust = -0.1,
          vjust = 1.5,
          size = 3,
          color = vline_color
        )
    }
  }

  # ── 16. Axis limits ───────────────────────────────────────────────────────────
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }

  # ── 17. Annotations and theme ─────────────────────────────────────────────────
  p <- p +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = cor_label,
      hjust = 1.05,
      vjust = -0.5,
      size = 3,
      color = "grey40"
    ) +
    labs(title = ttl, x = x_lab, y = y_lab) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0),
      axis.line = element_line(linewidth = 0.4, color = "grey50"),
      axis.ticks = element_line(linewidth = 0.3, color = "grey50"),
      panel.grid.major = element_line(color = "grey94", linewidth = 0.3),
      legend.position = if (color_by_density) "right" else "none"
    )

  return(p)
}
plot_multi_gp_curves <- function(
  gp_vector,
  loading_matrix,
  protein_matrix,
  protein,
  n_cells = 6000,
  seed = 42,
  nonzero_protein = FALSE,
  nonzero_loading = FALSE,
  n_quantile_bins = 10,
  show_ribbon = TRUE,
  line_size = 0.9,
  point_size = 2.0,
  colors = NULL,
  xlim = NULL,
  ylim = NULL,
  vline_protein_threshold = NULL,
  vline_color = "grey40",
  vline_label = NULL,
  title = NULL,
  x_label = NULL,
  y_label = NULL
) {
  # ── 1. Validate inputs ────────────────────────────────────────────────────────
  missing_gps <- setdiff(gp_vector, colnames(loading_matrix))
  if (length(missing_gps) > 0) {
    stop(
      "The following GPs are not in loading_matrix: ",
      paste(missing_gps, collapse = ", ")
    )
  }
  if (!protein %in% colnames(protein_matrix)) {
    stop("Protein '", protein, "' not found in protein_matrix columns.")
  }

  # ── 2. Shared cells ───────────────────────────────────────────────────────────
  shared_cells <- intersect(rownames(loading_matrix), rownames(protein_matrix))
  if (length(shared_cells) == 0) {
    stop("No shared row names between matrices.")
  }

  # ── 3. Protein expression vector (shared across all GPs) ─────────────────────
  x_all <- as.numeric(protein_matrix[shared_cells, protein])

  # ── 4. Compute quantile bins for each GP ─────────────────────────────────────
  curve_list <- lapply(gp_vector, function(gp) {
    df <- data.frame(
      cell = shared_cells,
      x = x_all,
      y = as.numeric(loading_matrix[shared_cells, gp])
    )

    # Apply filters
    if (nonzero_protein) {
      df <- df[df$x != 0, ]
    }
    if (nonzero_loading) {
      df <- df[df$y != 0, ]
    }
    if (nrow(df) == 0) {
      warning("No cells remain for GP '", gp, "' after filtering. Skipping.")
      return(NULL)
    }

    # Subsample
    if (!is.null(n_cells) && n_cells < nrow(df)) {
      set.seed(seed)
      df <- df[sample(nrow(df), size = n_cells, replace = FALSE), ]
    }

    # Cut into equal-frequency bins and summarise
    df %>%
      mutate(
        bin = cut(
          x,
          breaks = quantile(
            x,
            probs = seq(0, 1, length.out = n_quantile_bins + 1),
            na.rm = TRUE
          ),
          include.lowest = TRUE,
          labels = FALSE
        )
      ) %>%
      group_by(bin) %>%
      summarise(
        x_mean = mean(x, na.rm = TRUE),
        y_mean = mean(y, na.rm = TRUE),
        y_se = sd(y, na.rm = TRUE) / sqrt(n()),
        n = n(),
        .groups = "drop"
      ) %>%
      arrange(x_mean) %>%
      mutate(GP = gp)
  })

  # Remove NULLs (skipped GPs)
  curve_list <- Filter(Negate(is.null), curve_list)
  if (length(curve_list) == 0) {
    stop("No valid curves to plot.")
  }

  all_curves <- bind_rows(curve_list) %>%
    mutate(GP = factor(GP, levels = gp_vector))

  # ── 5. Colour palette ─────────────────────────────────────────────────────────
  n_gps <- length(unique(all_curves$GP))
  if (is.null(colors)) {
    # Use a colorblind-friendly qualitative palette; fall back to hue if >12 GPs
    if (n_gps <= 12 && requireNamespace("RColorBrewer", quietly = TRUE)) {
      colors <- RColorBrewer::brewer.pal(max(3, n_gps), "Paired")[seq_len(
        n_gps
      )]
    } else {
      colors <- scales::hue_pal()(n_gps)
    }
  }

  # ── 6. Labels ─────────────────────────────────────────────────────────────────
  x_lab <- if (!is.null(x_label)) x_label else paste0(protein, " expression")
  y_lab <- if (!is.null(y_label)) y_label else "Mean GP loading"
  ttl <- if (!is.null(title)) {
    title
  } else {
    paste0("GP loading vs ", protein, " — top ", n_gps, " GPs")
  }

  # ── 7. Build plot ─────────────────────────────────────────────────────────────
  p <- ggplot(all_curves, aes(x = x_mean, y = y_mean, color = GP, fill = GP))

  # Optional SE ribbon
  if (show_ribbon) {
    p <- p +
      geom_ribbon(
        aes(ymin = y_mean - y_se, ymax = y_mean + y_se),
        alpha = 0.12,
        color = NA,
        show.legend = FALSE
      )
  }

  # Connecting lines
  p <- p + geom_line(linewidth = line_size)

  # Bin-centre dots
  p <- p +
    geom_point(
      shape = 21,
      size = point_size,
      fill = "white",
      stroke = 0.8
    )

  # ── 8. Optional vertical threshold line ───────────────────────────────────────
  if (!is.null(vline_protein_threshold)) {
    p <- p +
      geom_vline(
        xintercept = vline_protein_threshold,
        linetype = "dashed",
        color = vline_color,
        linewidth = 0.5,
        alpha = 0.7
      )
    if (!is.null(vline_label)) {
      p <- p +
        annotate(
          "text",
          x = vline_protein_threshold,
          y = Inf,
          label = vline_label,
          hjust = -0.1,
          vjust = 1.5,
          size = 3,
          color = vline_color
        )
    }
  }

  # ── 9. Axis limits ────────────────────────────────────────────────────────────
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }

  # ── 10. Scales and theme ──────────────────────────────────────────────────────
  p <- p +
    scale_color_manual(values = colors, name = "GP") +
    scale_fill_manual(values = colors, name = "GP") +
    labs(title = ttl, x = x_lab, y = y_lab) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0),
      axis.line = element_line(linewidth = 0.4, color = "grey50"),
      axis.ticks = element_line(linewidth = 0.3, color = "grey50"),
      panel.grid.major = element_line(color = "grey94", linewidth = 0.3),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = 9)
    )

  return(p)
}
plot_factor_heatmap <- function(
  F_matrix,
  gp_vector,
  n_top = 5,
  min_abs_loading = 0.5,
  transpose = FALSE,
  show_border = FALSE,
  border_color = "grey80",
  border_size = 0.3,
  title = "Factor loadings — top genes per GP",
  x_label = NULL,
  y_label = NULL,
  midpoint = 0,
  low_color = "steelblue",
  mid_color = "white",
  high_color = "firebrick",
  limit = NULL,
  font_size = 9
) {
  # ── 1. Validate ───────────────────────────────────────────────────────────────
  missing_gps <- setdiff(gp_vector, colnames(F_matrix))
  if (length(missing_gps) > 0) {
    stop("GPs not found in F_matrix: ", paste(missing_gps, collapse = ", "))
  }

  # ── 2. Subset to requested GPs ────────────────────────────────────────────────
  F_sub <- F_matrix[, gp_vector, drop = FALSE]

  # ── 3. Select top N positive and top N negative genes per GP ─────────────────
  selected_genes <- lapply(gp_vector, function(gp) {
    vals <- F_sub[, gp]
    top_pos <- names(sort(vals, decreasing = TRUE))[seq_len(min(
      n_top,
      sum(vals > 0)
    ))]
    top_neg <- names(sort(vals, decreasing = FALSE))[seq_len(min(
      n_top,
      sum(vals < 0)
    ))]
    c(top_pos, top_neg)
  })
  selected_genes <- unique(unlist(selected_genes))

  # ── 4. Filter: drop genes whose max abs loading is below threshold ────────────
  if (min_abs_loading > 0) {
    max_abs <- apply(F_sub[selected_genes, , drop = FALSE], 1, function(x) {
      max(abs(x), na.rm = TRUE)
    })
    kept <- names(max_abs[max_abs >= min_abs_loading])
    n_dropped <- length(selected_genes) - length(kept)
    if (n_dropped > 0) {
      message(sprintf(
        "min_abs_loading = %.2f: removed %d genes, %d remaining.",
        min_abs_loading,
        n_dropped,
        length(kept)
      ))
    }
    selected_genes <- kept
    if (length(selected_genes) == 0) {
      stop(
        "No genes remain after min_abs_loading filter. Try lowering the threshold."
      )
    }
  }

  # ── 5. Order genes by hierarchical clustering ─────────────────────────────────
  hc_genes <- hclust(dist(F_sub[selected_genes, , drop = FALSE]))
  gene_order <- rownames(F_sub[selected_genes, , drop = FALSE])[hc_genes$order]

  # ── 6. Build long data frame ──────────────────────────────────────────────────
  plot_df <- F_sub[selected_genes, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene") %>%
    pivot_longer(cols = -Gene, names_to = "GP", values_to = "Loading") %>%
    mutate(
      GP = factor(GP, levels = gp_vector),
      Gene = factor(Gene, levels = gene_order)
    )

  # ── 7. Colour scale limits ────────────────────────────────────────────────────
  if (is.null(limit)) {
    limit <- max(abs(plot_df$Loading), na.rm = TRUE)
  }

  # ── 8. Set x / y aesthetics depending on transpose ───────────────────────────
  if (!transpose) {
    # Default: x = GP, y = Gene
    x_aes <- "GP"
    y_aes <- "Gene"
    x_lab <- if (!is.null(x_label)) x_label else "GP"
    y_lab <- if (!is.null(y_label)) y_label else "Gene"
    x_angle <- 45
    x_hjust <- 1
  } else {
    # Transposed: x = Gene, y = GP
    x_aes <- "Gene"
    y_aes <- "GP"
    x_lab <- if (!is.null(x_label)) x_label else "Gene"
    y_lab <- if (!is.null(y_label)) y_label else "GP"
    x_angle <- 90
    x_hjust <- 1
  }

  # ── 9. Build plot ─────────────────────────────────────────────────────────────
  p <- ggplot(
    plot_df,
    aes(x = .data[[x_aes]], y = .data[[y_aes]], fill = Loading)
  ) +
    geom_tile(
      color = if (show_border) border_color else NA,
      linewidth = if (show_border) border_size else 0
    ) +
    scale_fill_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = midpoint,
      limits = c(-limit, limit),
      name = "Loading"
    ) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_minimal(base_size = font_size) +
    theme(
      plot.title = element_text(size = font_size + 2, face = "bold", hjust = 0),
      axis.text.x = element_text(
        angle = x_angle,
        hjust = x_hjust,
        vjust = 1,
        size = font_size
      ),
      axis.text.y = element_text(size = font_size),
      axis.title = element_text(size = font_size + 1),
      panel.grid = element_blank(),
      legend.key.height = unit(1.2, "cm"),
      legend.title = element_text(size = font_size),
      legend.text = element_text(size = font_size - 1)
    )

  return(p)
}


#####################################################
#####################################################
#####################################################
### Loading Data
#####################################################
#####################################################
# load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
# flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(L_pm_filtered) <- paste0("GP", 1:ncol(L_pm_filtered))
colnames(F_pm_filtered) <- paste0("GP", 1:ncol(F_pm_filtered))
seurat_meta <- readRDS(paste0(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
protein_mat_normalized_lognorm <- readRDS(paste0(
  data_path,
  "protein_mat_normalized_lognorm.rds"
))
umap_result <- readRDS(paste0(data_path, "umap_result.rds"))
colnames(umap_result) <- c("UMAP_1", "UMAP_2")
umap_result <- umap_result[rownames(L_pm_filtered), ]
df_umap <- as.data.frame(umap_result)
# Protein_F_pm <- readRDS(file = paste0(data_path, "protein_projection_OLS_lognorm.rds"))
Protein_flash_result <- readRDS(
  file = paste0(
    data_path,
    "protein_flash_selected_summary_lognorm_backfit200.rds"
  )
)
Protein_F_pm <- Protein_flash_result$F_pm
colnames(Protein_F_pm) <- paste0("GP", 1:ncol(Protein_F_pm))
saveRDS(Protein_F_pm, file = paste0(data_path, "U_full.rds"))
cells_citeseq <- seurat_meta_filtered$cellID[seurat_meta_filtered$cite_seq]
# focusing on cells_citeseq in this CITE-seq analysis, subset to those cells
L_pm_filtered <- L_pm_filtered[cells_citeseq, , drop = FALSE]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[
  cells_citeseq,
  ,
  drop = FALSE
]
mean_shifted_log_expr <- apply(protein_mat_normalized_lognorm, 2, function(x) {
  mean(x)
})
saveRDS(
  mean_shifted_log_expr,
  file = paste0(data_path, "mean_shifted_log_expr_protein_full.rds")
)
umap_result <- umap_result[cells_citeseq, , drop = FALSE]
seurat_meta_filtered <- seurat_meta_filtered[cells_citeseq, , drop = FALSE]
isotype_proteins <- grep("^Isotype", rownames(Protein_F_pm), value = TRUE)
proteins_quality <- read.csv(
  paste0(data_path, "TableS4_citeseq_qc_20250513.csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  skip = 1
)
good_proteins <- proteins_quality$protein[
  proteins_quality$classification == "good"
]
# put some proteins of interest into good_proteins
good_proteins <- c(good_proteins, "IL2RA.CD25", "ITB7", "CD69")
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% c(isotype_proteins), ]
Protein_F_pm <- Protein_F_pm[rownames(Protein_F_pm) %in% good_proteins, ]
exclude_proteins <- c(
  "CD19",
  "CD34",
  "CD45.1",
  "CD45.2",
  "CD138",
  "TCRVA2",
  "TER119"
)
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% exclude_proteins, ]
thy11_proteins <- grep("THY1.1", rownames(Protein_F_pm), value = TRUE)
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% thy11_proteins, ]
D_lognorm <- diag(1 / apply(Protein_F_pm, 2, function(x) max(abs(x))))
Protein_F_pm <- Protein_F_pm %*% D_lognorm
colnames(Protein_F_pm) <- paste0("GP", 1:ncol(Protein_F_pm))
# there might be NaN due to 0/0, replace them with 0
Protein_F_pm[is.na(Protein_F_pm)] <- 0
saveRDS(Protein_F_pm, file = paste0(data_path, "U_pm_filtered.rds"))
# read the defined threshold for selected proteins
threshold_results_subset <- read.csv(
  paste0(data_path, "Thresholds_Selected_Proteins.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)
threshold_results_subset <- threshold_results_subset[, c(
  "Protein",
  "Threshold",
  "Threshold_manual"
)]
threshold_results_subset <- rbind(
  threshold_results_subset,
  data.frame(Protein = "CD62L", Threshold = 3, Threshold_manual = 3)
)
threshold_results_subset_manual <- data.frame(
  Protein = threshold_results_subset$Protein,
  Threshold = threshold_results_subset$Threshold_manual
)


#####################################################
#####################################################
#####################################################
### sorted abs effect for each GP
#####################################################
#####################################################
Protein_F_pm_raw <- Protein_flash_result$F_pm
# normalize each column to have max abs value of 1, so that we can compare the effect size across GPs
D_lognorm_raw <- diag(1 / apply(Protein_F_pm_raw, 2, function(x) max(abs(x))))
Protein_F_pm_raw <- Protein_F_pm_raw %*% D_lognorm_raw
# NaN should be replaced with 0
Protein_F_pm_raw[is.na(Protein_F_pm_raw)] <- 0
colnames(Protein_F_pm_raw) <- paste0("GP", 1:ncol(Protein_F_pm_raw))
plot_sorted_protein_effects <- function(k, Protein_F_pm) {
  sorted_proteins <- names(sort(abs(Protein_F_pm[, k]), decreasing = TRUE))
  plot_df <- data.frame(
    Rank = seq_along(sorted_proteins),
    Abs_Effect = abs(Protein_F_pm[sorted_proteins, k])
  )
  ggplot(plot_df, aes(x = Rank, y = Abs_Effect)) +
    geom_line(color = "steelblue") +
    labs(
      title = paste("Proteins sorted by absolute effect in GP", k),
      x = "Rank",
      y = "Absolute Effect (|loading|)"
    ) +
    theme_minimal()
}
plot_all_sorted_protein_effects <- function(Protein_F_pm) {
  n_gps <- ncol(Protein_F_pm)
  active_gps <- which(apply(Protein_F_pm, 2, function(x) max(abs(x))) > 0)
  plot_df <- do.call(
    rbind,
    lapply(active_gps, function(k) {
      sorted_proteins <- names(sort(abs(Protein_F_pm[, k]), decreasing = TRUE))
      data.frame(
        GP = paste0("GP", k),
        Rank = seq_along(sorted_proteins),
        Abs_Effect = abs(Protein_F_pm[sorted_proteins, k])
      )
    })
  )
  ggplot(plot_df, aes(x = Rank, y = Abs_Effect, group = GP)) +
    geom_line(alpha = 0.2, color = "steelblue") +
    labs(
      title = "Proteins sorted by absolute effect across all GPs",
      x = "Rank",
      y = "Absolute Effect (|loading|)"
    ) +
    theme_minimal()
}
plot_all_sorted_protein_effects(Protein_F_pm = Protein_F_pm_raw)
ggsave(
  plot = last_plot(),
  filename = paste0(figure_path, "All_Proteins_Sorted_by_Abs_Effect.pdf"),
  width = 6,
  height = 4
)
plot_all_sorted_protein_effects(Protein_F_pm = Protein_F_pm)
ggsave(
  plot = last_plot(),
  filename = paste0(figure_path, "Selected_Proteins_Sorted_by_Abs_Effect.pdf"),
  width = 6,
  height = 4
)


#####################################################
#####################################################
#####################################################
### Giant heatmap of protein programs
#####################################################
#####################################################
bk <- seq(-1, 1, length.out = 101)
cols <- colorRampPalette(c("#4575B4", "white", "#D73027"))(100)
# this is a giant heatmap, we should make the size large
pdf(
  paste0(figure_path, "Protein_Programs_heatmap.pdf"),
  width = 20,
  height = 40
)
pheatmap::pheatmap(
  t(Protein_F_pm),
  main = "Protein programs (normalized)",
  # color scale, blue negative, white zero, red positive
  color = cols,
  breaks = bk,
  border_color = "black"
)
dev.off()

pdf(
  paste0(figure_path, "Protein_Programs_heatmap_all.pdf"),
  width = 20,
  height = 40
)
pheatmap::pheatmap(
  t(Protein_F_pm_raw),
  main = "Protein programs (normalized)",
  # color scale, blue negative, white zero, red positive
  color = cols,
  breaks = bk,
  # remove x/y axis labels for better visualization
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = "black"
)
dev.off()

#####################################################
#####################################################
#####################################################
### Giant heatmap of protein programs (simplified)
### - drop proteins/GPs with |loading| < 0 everywhere
#####################################################
#####################################################
threshold_simplified <- 0
keep_rows_simplified <- apply(Protein_F_pm, 1, function(v) {
  any(abs(v) > threshold_simplified, na.rm = TRUE)
})
Protein_F_pm_simplified <- Protein_F_pm[keep_rows_simplified, , drop = FALSE]
keep_cols_simplified <- apply(Protein_F_pm_simplified, 2, function(v) {
  any(abs(v) > threshold_simplified, na.rm = TRUE)
})
Protein_F_pm_simplified <- Protein_F_pm_simplified[,
  keep_cols_simplified,
  drop = FALSE
]

pdf(
  paste0(figure_path, "Protein_Programs_heatmap_all_simplified.pdf"),
  width = 20,
  height = 40
)
pheatmap::pheatmap(
  t(Protein_F_pm_simplified),
  main = sprintf(
    "Protein programs (simplified): %d proteins x %d GPs, |loading| >= %.2f in >=1 entry per row/col",
    nrow(Protein_F_pm_simplified),
    ncol(Protein_F_pm_simplified),
    threshold_simplified
  ),
  color = cols,
  breaks = bk,
  border_color = "black"
)
dev.off()


pdf(
  paste0(figure_path, "Protein_Programs_heatmap_all_simplified_sparse.pdf"),
  width = 20,
  height = 40
)
# sparse version: |loading| < 0.5 → white
sparse_cutoff <- 0.5
bk_sparse <- c(
  seq(-1, -sparse_cutoff, length.out = 26),
  seq(-sparse_cutoff, sparse_cutoff, length.out = 51),
  seq(sparse_cutoff, 1, length.out = 26)
)
bk_sparse <- unique(bk_sparse)
cols_sparse <- c(
  colorRampPalette(c("#4575B4", "white"))(25),
  rep("white", 50),
  colorRampPalette(c("white", "#D73027"))(25)
)
pheatmap::pheatmap(
  t(Protein_F_pm_simplified),
  main = sprintf(
    "Protein programs (sparse, |loading| < %.2f = white): %d proteins x %d GPs",
    sparse_cutoff,
    nrow(Protein_F_pm_simplified),
    ncol(Protein_F_pm_simplified)
  ),
  color = cols_sparse,
  breaks = bk_sparse,
  border_color = "black"
)
dev.off()

#### There are a few GPs that could represent contamination of other cell types.
#### We will take them out and redo the heatmap

GP_contamination <- c("GP40", "GP50", "GP55", "GP188")
Protein_F_pm_simplified_no_contamination <- Protein_F_pm_simplified[,
  !colnames(Protein_F_pm_simplified) %in% GP_contamination,
  drop = FALSE
]
pdf(
  paste0(
    figure_path,
    "Protein_Programs_heatmap_all_simplified_sparse_no_contamination.pdf"
  ),
  width = 20,
  height = 40
)
pheatmap::pheatmap(
  t(Protein_F_pm_simplified_no_contamination),
  main = sprintf(
    "Protein programs (sparse, no contamination GPs): %d proteins x %d GPs",
    nrow(Protein_F_pm_simplified_no_contamination),
    ncol(Protein_F_pm_simplified_no_contamination)
  ),
  color = cols_sparse,
  breaks = bk_sparse,
  border_color = "black"
)
dev.off()

# take a subset of well-aligned GPs
# > well_aligned_gps
#  [1] "GP10"  "GP26"  "GP68"  "GP171" "GP58"  "GP8"   "GP30"  "GP27"  "GP170" "GP80"  "GP35"  "GP12"  "GP3"   "GP29"  "GP77"  "GP22"
# [17] "GP25"  "GP41"  "GP181" "GP63"  "GP107" "GP23"  "GP126" "GP127" "GP192" "GP159" "GP153"
well_aligned_gps <- c(
  "GP10",
  "GP26",
  "GP68",
  "GP171",
  "GP58",
  "GP8",
  "GP30",
  "GP27",
  "GP170",
  "GP80",
  "GP35",
  "GP12",
  "GP3",
  "GP29",
  "GP77",
  "GP22",
  "GP25",
  "GP41",
  "GP181",
  "GP63",
  "GP107",
  "GP23",
  "GP126",
  "GP127",
  "GP192",
  "GP159",
  "GP153"
)
Protein_F_pm_simplified_well_aligned <- Protein_F_pm_simplified[,
  well_aligned_gps
]
quantile(colSums(abs(Protein_F_pm_simplified_well_aligned) > 0.5))
quantile(rowSums(abs(Protein_F_pm_simplified_well_aligned) > 0.5))

#####################################################
#####################################################
#####################################################
### Volcano plots of protein programs
#####################################################
#####################################################
# keep only those in Protein_F_pm
citeseq_F <- Protein_F_pm
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[, rownames(
  Protein_F_pm
)]
mean_shifted_log_expr <- apply(protein_mat_normalized_lognorm, 2, function(x) {
  mean(x)
})
saveRDS(
  mean_shifted_log_expr,
  file = paste0(data_path, "mean_shifted_log_expr_protein.rds")
)
pdf(paste0(figure_path, "CiteSeq_volcano_plots.pdf"), width = 6, height = 6)
for (i in 1:ncol(Protein_F_pm)) {
  p <- plot_scatterplot(i, n_higlight = 10)
  print(p)
}
dev.off()


#####################################################
#####################################################
#####################################################
### Table of top positive and negative protein markers for each GP
#####################################################
#####################################################
marker_threshold <- 0.5
gp_marker_list <- list()
gp_marker_sign_list <- list()
for (i in 1:ncol(Protein_F_pm)) {
  gp_name <- colnames(Protein_F_pm)[i]
  marker_values <- Protein_F_pm[, i]

  # 1. Identify all markers meeting the absolute threshold
  # No longer restricted to top 5
  marker_genes <- rownames(Protein_F_pm)[which(
    abs(marker_values) >= marker_threshold
  )]

  if (length(marker_genes) == 0) {
    gp_marker_sign_list[[gp_name]] <- list(
      pos = character(0),
      neg = character(0)
    )
    next
  }

  # 2. Sort markers by magnitude for better table readability
  marker_genes <- marker_genes[order(-abs(marker_values[marker_genes]))]
  gp_marker_list[[gp_name]] <- marker_genes

  # 3. Separate into positive and negative sets based on sign
  marker_genes_pos <- marker_genes[which(
    marker_values[marker_genes] >= marker_threshold
  )]
  marker_genes_neg <- marker_genes[which(
    marker_values[marker_genes] <= -marker_threshold
  )]

  gp_marker_sign_list[[gp_name]] <- list(
    pos = marker_genes_pos,
    neg = marker_genes_neg
  )
}
# 4. Build and Format the Table
df_markers <- marker_list_to_df(gp_marker_sign_list)
df_markers$GP_idx <- as.integer(sub("^GP", "", df_markers$Set))
df_markers <- df_markers[order(df_markers$GP_idx), ]
# Select only the marker columns for the final table
df_final <- df_markers[, c("Positive", "Negative")]
# 5. Generate LaTeX and HTML outputs
latex_code <- knitr::kable(
  df_final,
  format = "latex",
  booktabs = TRUE,
  longtable = TRUE,
  caption = "Comprehensive positive and negative protein markers (threshold >= 0.5)."
)
tab <- knitr::kable(df_final, format = "html") |>
  kableExtra::kable_styling(
    full_width = FALSE,
    bootstrap_options = c("striped", "hover")
  )
# Save results
html_file <- paste0(figure_path, "CITEseq_markers_full.html")
png_file <- paste0(figure_path, "CITEseq_markers_full.png")
kableExtra::save_kable(tab, file = html_file)
webshot2::webshot(html_file, file = png_file, zoom = 2)
# save df_markers into data_path for later use
saveRDS(df_markers, file = paste0(data_path, "CITEseq_markers_full.rds"))


#####################################################
#####################################################
#####################################################
### Table of top positive and negative protein markers for each GP
#####################################################
#####################################################
df_markers <- df_markers[order(df_markers$GP_idx), ]
# Select only the marker columns for the final table
df_final <- df_markers[, c("Positive", "Negative")]
# 5. Generate LaTeX and HTML outputs
latex_code <- knitr::kable(
  df_final,
  format = "latex",
  booktabs = TRUE,
  longtable = TRUE,
  caption = "Comprehensive positive and negative protein markers (threshold >= 0.5)."
)
tab <- knitr::kable(df_final, format = "html") |>
  kableExtra::kable_styling(
    full_width = FALSE,
    bootstrap_options = c("striped", "hover")
  )
# Save results
html_file <- paste0(figure_path, "CITEseq_markers_full.html")
png_file <- paste0(figure_path, "CITEseq_markers_full.png")


#####################################################
#####################################################
#####################################################
### KLRG1+/- for Treg and CD8
#####################################################
#####################################################
klrg1_threshold <- threshold_results_subset_manual$Threshold[
  threshold_results_subset_manual$Protein == "KLRG1"
]
cd8_split <- get_klrg1_split(
  "CD8",
  seurat_meta_filtered,
  protein_mat_normalized_lognorm,
  klrg1_threshold
)
treg_split <- get_klrg1_split(
  "Treg",
  seurat_meta_filtered,
  protein_mat_normalized_lognorm,
  klrg1_threshold
)
diff_factors_CD8 <- run_checked_dge(
  cd8_split,
  F_pm_filtered,
  L_pm_filtered,
  "CD8"
) %>%
  rename(mean_change_CD8 = mean_change_loadings, AveExpr_CD8 = AveExpr)
diff_factors_Treg <- run_checked_dge(
  treg_split,
  F_pm_filtered,
  L_pm_filtered,
  "Treg"
) %>%
  rename(mean_change_Treg = mean_change_loadings, AveExpr_Treg = AveExpr)
diff_factors_merged <- inner_join(
  diff_factors_Treg,
  diff_factors_CD8,
  by = "SYMBOL"
)
target_gps_to_show <- c("GP6", "GP10", "GP12", "GP27", "GP68", "GP58")
p_compare <- plot_target_gps(
  df = diff_factors_merged,
  x_var = mean_change_CD8,
  y_var = mean_change_Treg,
  label_var = SYMBOL,
  target_gps = target_gps_to_show,
  background_alpha = 0.8,
  x_limits = c(-0.2, 0.4),
  y_limits = c(-0.2, 0.4),
  highlight_color = c(
    "GP10" = "darkorange2",
    "GP27" = "deeppink",
    "GP6" = "deeppink",
    "GP68" = "deeppink",
    "GP12" = "deeppink",
    "GP58" = "darkorange2"
  ),
  title = "KLRG1 Modulation: CD8 vs Treg",
  xlab = "Effect Size in CD8 (KLRG1+ - KLRG1-)",
  ylab = "Effect Size in Treg (KLRG1+ - KLRG1-)"
) +
  theme_bw()
print(p_compare)
ggsave(
  paste0(figure_path, "KLRG1_CD8_vs_Treg_comparison.pdf"),
  p_compare,
  width = 7,
  height = 6
)


#####################################################
#####################################################
#####################################################
### KLRG1+/- for CD4 and CD8
#####################################################
#####################################################
klrg1_threshold <- threshold_results_subset_manual$Threshold[
  threshold_results_subset_manual$Protein == "KLRG1"
]
cd8_split <- get_klrg1_split(
  "CD8",
  seurat_meta_filtered,
  protein_mat_normalized_lognorm,
  klrg1_threshold
)
CD4_split <- get_klrg1_split(
  "CD4",
  seurat_meta_filtered,
  protein_mat_normalized_lognorm,
  klrg1_threshold
)
diff_factors_CD8 <- run_checked_dge(
  cd8_split,
  F_pm_filtered,
  L_pm_filtered,
  "CD8"
) %>%
  rename(mean_change_CD8 = mean_change_loadings, AveExpr_CD8 = AveExpr)
diff_factors_CD4 <- run_checked_dge(
  CD4_split,
  F_pm_filtered,
  L_pm_filtered,
  "CD4"
) %>%
  rename(mean_change_CD4 = mean_change_loadings, AveExpr_CD4 = AveExpr)
diff_factors_merged <- inner_join(
  diff_factors_CD4,
  diff_factors_CD8,
  by = "SYMBOL"
)
target_gps_to_show <- c("GP10", "GP58", "GP25", "GP26", "GP43")
p_compare <- plot_target_gps(
  df = diff_factors_merged,
  x_var = mean_change_CD8,
  y_var = mean_change_CD4,
  label_var = SYMBOL,
  target_gps = target_gps_to_show,
  background_alpha = 0.8,
  x_limits = c(-0.2, 0.4),
  y_limits = c(-0.2, 0.4),
  highlight_color = c(
    "GP10" = "darkorange2",
    "GP25" = "blue",
    "GP43" = "blue",
    "GP26" = "blue",
    "GP58" = "darkorange2"
  ),
  title = "KLRG1 Modulation: CD8 vs CD4",
  xlab = "Effect Size in CD8 (KLRG1+ - KLRG1-)",
  ylab = "Effect Size in CD4 (KLRG1+ - KLRG1-)"
) +
  theme_bw()
print(p_compare)
ggsave(
  paste0(figure_path, "KLRG1_CD8_vs_CD4_comparison.pdf"),
  p_compare,
  width = 7,
  height = 6
)


#####################################################
#####################################################
#####################################################
### Swarm plot for CD69
#####################################################
#####################################################
pve_results_matrix <- readRDS(paste0(
  "figures/Figures_Protein_GP/pve_results_matrix.rds"
))
colnames(pve_results_matrix) <- paste0("GP", 1:ncol(pve_results_matrix))
p_cd69 <- my_swarm_plot(
  input_matrix = pve_results_matrix["CD69", , drop = FALSE],
  title = "Predictive Power of CD69 across GPs",
  y_label = "PVE (R-squared)",
  label_threshold_quantile = 0.1,
  top_hit = 13,
  global_threshold = 0.03,
  manual_colors = c("steelblue", "darkorange"),
  n_colors = 2,
  use_beeswarm = F,
  jitter_width = 0.2,
  label_direction = "both",
  label_nudge_y = 0.005,
  label_force = 2,
  label_force_pull = 0.5
)
print(p_cd69)


#####################################################
#####################################################
#####################################################
### CD69 influence a few GPs
#####################################################
#####################################################
# find the col names of the top 13
cd69_top_gps <- names(sort(pve_results_matrix["CD69", ], decreasing = TRUE)[
  1:13
])
# for (gp in cd69_top_gps) {
#   p <- plot_gp_vs_protein(
#     loading_matrix = L_pm_filtered,
#     protein_matrix = protein_mat_normalized_lognorm,
#     ylim = c(0, 1),
#     gp = gp,
#     point_alpha = 0.3,
#     protein = "CD69",
#     vline_protein_threshold = 3.5,
#     n_cells = 6000,
#     add_smoother = FALSE,
#     smoother_se = FALSE,
#     smoother_method = "loess",
#     add_quantile_smoother = TRUE,
#     color_by_density = TRUE
#   )
#   print(p)
#   ggsave(
#     p,
#     filename = paste0(figure_path, "CD69_vs_", gp, ".pdf"),
#     width = 6,
#     height = 5
#   )
# }
p_multi <- plot_multi_gp_curves(
  gp_vector = cd69_top_gps,
  loading_matrix = L_pm_filtered,
  protein_matrix = protein_mat_normalized_lognorm,
  protein = "CD69",
  n_cells = 6000,
  nonzero_protein = FALSE,
  n_quantile_bins = 10,
  show_ribbon = TRUE,
  ylim = c(0, 0.5),
  vline_protein_threshold = 3.5,
  vline_label = "threshold"
)
print(p_multi)
ggsave(
  p_multi,
  filename = paste0(figure_path, "CD69_top_GPs_curves.pdf"),
  width = 8,
  height = 5
)


#####################################################
#####################################################
#####################################################
### Dot plot of top GPs associated with CD69
#####################################################
#####################################################
# scale each column to have max absolute value of 1
D <- diag(1 / apply(F_pm_filtered, 2, function(x) max(abs(x), na.rm = TRUE)))
F_pm_filtered_scaled <- F_pm_filtered %*% D
colnames(F_pm_filtered_scaled) <- paste0("GP", 1:ncol(F_pm_filtered_scaled))
annotation_heatmap_mat <- F_pm_filtered_scaled[, cd69_top_gps]
annotation_heatmap(
  annotation_heatmap_mat,
  n = 3,
  dims = cd69_top_gps,
  select_features = "largest",
  zero_value = 0.1,
  font_size = 9
)


#####################################################
#####################################################
#####################################################
### Heatmap of top GPs associated with CD69
#####################################################
#####################################################

cd69_top_gps_subset <- c(
  "GP35",
  "GP6",
  "GP170",
  "GP26",
  "GP58",
  "GP171",
  "GP63",
  "GP62",
  "GP3",
  "GP29"
)

# compute Spearman correlation between CD69 and each GP loading
shared_cells_cd69 <- intersect(
  rownames(L_pm_filtered),
  rownames(protein_mat_normalized_lognorm)
)
cd69_expr_vec <- protein_mat_normalized_lognorm[shared_cells_cd69, "CD69"]
cd69_corr <- sapply(cd69_top_gps_subset, function(gp) {
  cor(L_pm_filtered[shared_cells_cd69, gp], cd69_expr_vec, method = "spearman")
})

# sort: increasing so most-correlated GP ends up at top of y-axis
cd69_top_gps_sorted <- names(sort(cd69_corr, decreasing = FALSE))

p_heatmap <- plot_factor_heatmap(
  F_matrix = F_pm_filtered_scaled,
  gp_vector = cd69_top_gps_sorted,
  n_top = 5,
  font_size = 9,
  transpose = TRUE,
  low_color = "#4DAF4A",
  mid_color = "white",
  high_color = "#984EA3"
)

# correlation strip (one tile per GP, same y-factor order)
corr_strip_df <- data.frame(
  GP = factor(cd69_top_gps_sorted, levels = cd69_top_gps_sorted),
  Correlation = cd69_corr[cd69_top_gps_sorted],
  x = "Corr"
)
corr_limit <- max(abs(corr_strip_df$Correlation))
p_corr_strip <- ggplot(corr_strip_df, aes(x = x, y = GP, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "royalblue",
    mid = "white",
    high = "tomato",
    midpoint = 0,
    limits = c(-corr_limit, corr_limit),
    name = "Corr\n(CD69)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

p_combined <- p_corr_strip +
  p_heatmap +
  patchwork::plot_layout(widths = c(0.06, 1), guides = "collect")

ggsave(
  p_combined,
  filename = paste0(figure_path, "CD69_top_GPs_heatmap.pdf"),
  width = 11,
  height = 5
)

############### CD69_top_GP in level1 and in organ_simplified

cells_for_heatmap <- intersect(
  rownames(L_pm_filtered),
  rownames(seurat_meta_filtered)
)
L_cd69_sub <- L_pm_filtered[
  cells_for_heatmap,
  cd69_top_gps_sorted,
  drop = FALSE
]
meta_hm <- seurat_meta_filtered[
  cells_for_heatmap,
  c("annotation_level1", "organ_simplified")
]

mean_loading_long <- function(L_mat, group_vec, gp_levels) {
  as.data.frame(L_mat) %>%
    mutate(group = group_vec) %>%
    pivot_longer(cols = -group, names_to = "GP", values_to = "Loading") %>%
    group_by(group, GP) %>%
    summarise(mean_loading = mean(Loading, na.rm = TRUE), .groups = "drop") %>%
    mutate(GP = factor(GP, levels = gp_levels))
}

make_mean_loading_heatmap <- function(df, title) {
  fill_max <- max(df$mean_loading, na.rm = TRUE)
  ggplot(df, aes(x = group, y = GP, fill = mean_loading)) +
    geom_tile() +
    scale_fill_gradient(
      low = "white",
      high = "firebrick",
      limits = c(0, fill_max),
      name = "Mean\nloading"
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank()
    )
}

df_level1 <- mean_loading_long(
  L_cd69_sub,
  meta_hm$annotation_level1,
  cd69_top_gps_sorted
)
p_hm_level1 <- make_mean_loading_heatmap(
  df_level1,
  "Mean GP loading by cell type (level1)"
)
ggsave(
  p_hm_level1,
  filename = paste0(figure_path, "CD69_top_GPs_mean_loading_level1.pdf"),
  width = 7,
  height = 5
)

df_organ <- mean_loading_long(
  L_cd69_sub,
  meta_hm$organ_simplified,
  cd69_top_gps_sorted
)
p_hm_organ <- make_mean_loading_heatmap(
  df_organ,
  "Mean GP loading by tissue (organ_simplified)"
)
ggsave(
  p_hm_organ,
  filename = paste0(figure_path, "CD69_top_GPs_mean_loading_organ.pdf"),
  width = 9,
  height = 5
)

################################################################################
################################################################################
#### Quantify the agreement between protein-signature and high-loading cells.

# --- eligible cell exclusions (match gated_protein_loading_plot.R) ---
thymocyte_cells <- seurat_meta_filtered %>%
  filter(annotation_level1 == "thymocyte") %>%
  pull(cellID)
proliferating_cells <- seurat_meta_filtered %>%
  filter(annotation_level2_group == "proliferating") %>%
  pull(cellID)
miniverse_cells <- seurat_meta_filtered %>%
  filter(annotation_level2_group == "miniverse") %>%
  pull(cellID)
exclude_cells_all <- c(thymocyte_cells, proliferating_cells, miniverse_cells)

# --- curated proteins (match gated_protein_loading_plot.R) ---
select_proteins_for_alignment <- setdiff(
  c(
    proteins_quality$protein[proteins_quality$classification == "good"],
    "IL2RA.CD25",
    "ITB7",
    "CD69"
  ),
  c(
    c("CD19", "CD34", "CD45.1", "CD45.2", "CD138", "TCRVA2", "TER119"),
    grep("THY1.1", rownames(Protein_F_pm), value = TRUE),
    grep("^Isotype", rownames(Protein_F_pm), value = TRUE)
  )
)

# --- manually curated marker overrides (mirror gated_protein_loading_plot.R) ---
df_markers2 <- df_markers
df_markers2$Positive[8] <- "TCRVG3"
df_markers2$Negative[8] <- ""
df_markers2$Negative[22] <- "CD4, CD8B, TCRGD"
df_markers2$Positive[23] <- "ITB7, CD103, CD44"
df_markers2$Negative[23] <- "CD29, ITA4.CD49D"
df_markers2$Positive[29] <- "CD8A"
df_markers2$Negative[29] <- "CD4, CD8B"
df_markers2$Positive[30] <- "CD11A, CD49B, CD38"
df_markers2$Negative[30] <- "CD8B, CD8A, CD4"
df_markers2$Positive[
  41
] <- "CD55.DAF, CD4, CD44, CD45RB, CD62L, CD31, CD2, IL7RA.CD127, CD27, SCA1, CD5"
df_markers2$Positive[57] <- "CD44, ICOS.CD278"
df_markers2$Negative[57] <- "CD62L"
df_markers2$Positive[68] <- "IL2RA.CD25, FR4, GITR.CD357, NEUROPILIN1.CD304"
df_markers2$Negative[68] <- ""
df_markers2$Positive[80] <- "CD29, CD44, ITA4.CD49D"
df_markers2$Negative[80] <- "ITB7, CD103"
df_markers2$Positive[170] <- "ITB7, CD103, CD4, CD38"
df_markers2$Negative[170] <- "GITR.CD357, CD62L"
df_markers2$Positive[171] <- "CD62L"
df_markers2$Negative[171] <- "CD44"

# --- alignment metrics:
#     m = |protein-gated|, top_m = top-m cells by GP loading,
#     ratio = |gated ∩ top_m| / m
compute_alignment_scores <- function(
  df_m,
  protein_mat,
  loading_mat,
  threshold_df,
  selected_proteins = NULL,
  exclude_cells = NULL,
  missing_threshold_action = "skip",
  gp_pos_threshold = 0.1
) {
  eligible <- rownames(loading_mat)
  if (!is.null(exclude_cells)) {
    eligible <- setdiff(eligible, exclude_cells)
  }
  eligible <- intersect(eligible, rownames(protein_mat))
  sub_prot <- protein_mat[eligible, , drop = FALSE]
  sub_load <- loading_mat[eligible, , drop = FALSE]

  apply_gate <- function(markers, is_pos) {
    keep <- rep(TRUE, nrow(sub_prot))
    applied <- 0L
    for (p in markers) {
      if (!(p %in% colnames(sub_prot))) {
        next
      }
      thresh <- NULL
      if (p %in% threshold_df$Protein) {
        thresh <- threshold_df$Threshold[threshold_df$Protein == p]
      } else if (missing_threshold_action == "median") {
        thresh <- median(sub_prot[, p], na.rm = TRUE)
      }
      if (is.null(thresh)) {
        next
      }
      keep <- if (is_pos) {
        keep & (sub_prot[, p] > thresh)
      } else {
        keep & (sub_prot[, p] <= thresh)
      }
      applied <- applied + 1L
    }
    list(keep = keep, applied = applied)
  }

  gps <- rownames(df_m)
  N_eligible <- nrow(sub_load)
  res <- data.frame(
    GP = gps,
    n_pos = NA_integer_,
    n_neg = NA_integer_,
    n_applied = NA_integer_,
    m = NA_integer_,
    k = NA_integer_,
    ratio = NA_real_,
    mean_load_gated = NA_real_,
    mean_load_all = NA_real_,
    load_ratio = NA_real_,
    prop_pos_gated = NA_real_,
    prop_pos_all = NA_real_,
    z_enrich = NA_real_,
    p_enrich = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(gps)) {
    gp <- gps[i]
    if (!(gp %in% colnames(sub_load))) {
      next
    }

    pos_markers <- strsplit(df_m[i, "Positive"], ", ")[[1]]
    neg_markers <- strsplit(df_m[i, "Negative"], ", ")[[1]]
    if (!is.null(selected_proteins)) {
      pos_markers <- intersect(pos_markers, selected_proteins)
      neg_markers <- intersect(neg_markers, selected_proteins)
    }
    pos_markers <- pos_markers[pos_markers != ""]
    neg_markers <- neg_markers[neg_markers != ""]

    pg <- apply_gate(pos_markers, TRUE)
    ng <- apply_gate(neg_markers, FALSE)
    keep_vec <- pg$keep & ng$keep
    m <- sum(keep_vec)
    n_applied <- pg$applied + ng$applied

    loadings <- sub_load[, gp]
    mean_all <- mean(loadings, na.rm = TRUE)
    prop_pos_all <- mean(loadings > gp_pos_threshold, na.rm = TRUE)

    res$n_pos[i] <- length(pos_markers)
    res$n_neg[i] <- length(neg_markers)
    res$n_applied[i] <- n_applied
    res$m[i] <- m
    res$mean_load_all[i] <- mean_all
    res$prop_pos_all[i] <- prop_pos_all

    # vacuous gate or empty gate: ratio = 0 (no agreement),
    # loading-based metrics stay NA (no gated cells to average)
    if (n_applied == 0L || m == 0L) {
      res$k[i] <- 0L
      res$ratio[i] <- 0
      next
    }

    gated_cells <- rownames(sub_prot)[keep_vec]
    top_m_cells <- rownames(sub_load)[order(
      loadings,
      decreasing = TRUE
    )[seq_len(m)]]

    k <- length(intersect(gated_cells, top_m_cells))
    mean_gated <- mean(loadings[keep_vec], na.rm = TRUE)
    load_ratio <- if (is.finite(mean_all) && mean_all > 0) {
      mean_gated / mean_all
    } else {
      NA_real_
    }

    # one-sided z under H0: gated cells are a random size-m subsample of eligible
    # SE uses finite-population correction; valid for m < N
    sd_all <- sd(loadings, na.rm = TRUE)
    se_null <- if (
      is.finite(sd_all) && sd_all > 0 && N_eligible > 1 && m < N_eligible
    ) {
      sd_all / sqrt(m) * sqrt((N_eligible - m) / (N_eligible - 1))
    } else {
      NA_real_
    }
    z_enrich <- if (is.finite(se_null) && se_null > 0) {
      (mean_gated - mean_all) / se_null
    } else {
      NA_real_
    }
    p_enrich <- if (is.finite(z_enrich)) {
      pnorm(z_enrich, lower.tail = FALSE)
    } else {
      NA_real_
    }

    prop_pos_gated <- mean(loadings[keep_vec] > gp_pos_threshold, na.rm = TRUE)

    res$k[i] <- k
    res$ratio[i] <- k / m
    res$mean_load_gated[i] <- mean_gated
    res$load_ratio[i] <- load_ratio
    res$prop_pos_gated[i] <- prop_pos_gated
    res$z_enrich[i] <- z_enrich
    res$p_enrich[i] <- p_enrich
  }
  res
}

# df_markers rownames are "GP..."; ensure L_pm_filtered colnames match
# (gated_protein_loading_plot.R leaves them as "K...", so rename if needed)
if (!all(rownames(df_markers) %in% colnames(L_pm_filtered))) {
  colnames(L_pm_filtered) <- paste0("GP", seq_len(ncol(L_pm_filtered)))
}

scores_default <- compute_alignment_scores(
  df_m = df_markers,
  protein_mat = protein_mat_normalized_lognorm,
  loading_mat = L_pm_filtered,
  threshold_df = threshold_results_subset_manual,
  selected_proteins = select_proteins_for_alignment,
  exclude_cells = exclude_cells_all,
  missing_threshold_action = "skip"
)
scores_manual <- compute_alignment_scores(
  df_m = df_markers2,
  protein_mat = protein_mat_normalized_lognorm,
  loading_mat = L_pm_filtered,
  threshold_df = threshold_results_subset_manual,
  selected_proteins = select_proteins_for_alignment,
  exclude_cells = exclude_cells_all,
  missing_threshold_action = "skip"
)


# Which GPs have good alignment under default?
# choose those with prop_pos_gated > 0.25 and prop_pos_gated > prop_pos_all
scores_default %>%
  filter(prop_pos_gated > 0.25 & prop_pos_gated > prop_pos_all)

scores_manual %>%
  filter(prop_pos_gated > 0.25 & prop_pos_gated > prop_pos_all)

# Format scores_default into a nice csv table for presentation
format_scores_table <- function(scores, df_m) {
  tbl <- scores

  # attach the marker definitions used for gating (Positive / Negative)
  tbl$Positive <- df_m[tbl$GP, "Positive"]
  tbl$Negative <- df_m[tbl$GP, "Negative"]

  # flag well-aligned GPs using the same criterion as above
  tbl$well_aligned <- with(
    tbl,
    !is.na(prop_pos_gated) &
      prop_pos_gated >= 0.25 &
      prop_pos_gated > prop_pos_all
  )

  # round the numeric metrics for readability
  tbl$ratio <- round(tbl$ratio, 3)
  tbl$mean_load_gated <- round(tbl$mean_load_gated, 3)
  tbl$mean_load_all <- round(tbl$mean_load_all, 3)
  # log2 fold change of gated vs. overall mean loading
  tbl$log_fc <- round(log2(tbl$load_ratio), 2)
  tbl$prop_pos_gated <- round(tbl$prop_pos_gated, 3)
  tbl$prop_pos_all <- round(tbl$prop_pos_all, 3)

  # put well-aligned programs first, then by strength of positive-marker
  # enrichment in the gated cells (NA proportions sort to the bottom)
  tbl <- tbl[
    order(
      -tbl$well_aligned,
      -tbl$prop_pos_gated,
      na.last = TRUE
    ),
  ]

  # rename columns to presentation-friendly headers
  tbl <- tbl[, c(
    "GP",
    "Positive",
    "Negative",
    "n_pos",
    "n_neg",
    "n_applied",
    "m",
    "k",
    "ratio",
    "mean_load_gated",
    "mean_load_all",
    "log_fc",
    "prop_pos_gated",
    "prop_pos_all",
    "well_aligned"
  )]
  colnames(tbl) <- c(
    "GP",
    "Positive markers",
    "Negative markers",
    "# pos markers",
    "# neg markers",
    "# gates applied",
    "# gated cells (m)",
    "# top-m gated (k)",
    "Overlap ratio (k/m)",
    "Mean loading (gated)",
    "Mean loading (all)",
    "Log fold change",
    "Prop. positive (gated)",
    "Prop. positive (all)",
    "Well aligned"
  )
  rownames(tbl) <- NULL
  tbl
}

scores_default_table <- format_scores_table(scores_default, df_markers)
scores_manual_table <- format_scores_table(scores_manual, df_markers2)


write.csv(
  scores_default_table,
  file = paste0(figure_path, "CITEseq_alignment_scores_default.csv"),
  row.names = FALSE
)
write.csv(
  scores_manual_table,
  file = paste0(figure_path, "CITEseq_alignment_scores_manual.csv"),
  row.names = FALSE
)
