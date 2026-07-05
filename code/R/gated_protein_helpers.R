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

# ============================================================
# Protein-gate alignment scores (Extended Data "Protein Gating" table).
#
# For each GP: build a protein gate from its curated marker signature
# (positive markers strictly above threshold, negative markers at/below),
# count the gated cells m, compare them against the top-m cells by GP loading,
# and score how much the gate enriches for high GP loading. Ported verbatim
# from the pre-refactor Figure_CITEseq.R (compute_alignment_scores /
# format_scores_table); the manual-marker run produced the published
# data/CITEseq_alignment_scores_manual.csv.
#
#   m     = |protein-gated cells|
#   top_m = the m cells with highest GP loading
#   ratio = |gated ∩ top_m| / m
# ============================================================
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

# Turn compute_alignment_scores() output into the presentation table: attach the
# marker definitions, flag well-aligned GPs (prop_pos_gated >= 0.25 and greater
# than prop_pos_all), round metrics, order (well-aligned first, then by
# prop_pos_gated), and rename columns to friendly headers.
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
