# Plotting helpers shared by script-refactor/Figure2.R and FigureS2.R
# (T-cell lineage GPs). Ported from script/Figure_Lineage.R.

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
      dplyr::left_join(df_load, by = c("Lineage", "GP")) %>%
      dplyr::left_join(
        data.frame(GP = names(overall_loading_vec), Overall_Loading = unname(overall_loading_vec)),
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
    if (!has_loading) stop("filter_positive requires loading_mat and overall_loading_vec")
    df <- df %>% dplyr::filter(Loading >= Overall_Loading | .key %in% forced_keys)
  }

  if (!is.null(lineage_order)) {
    existing <- intersect(lineage_order, unique(as.character(df$Lineage)))
    df$Lineage <- factor(df$Lineage, levels = existing)
  }

  df <- df %>%
    dplyr::group_by(Lineage) %>%
    dplyr::arrange(dplyr::desc(Value)) %>%
    dplyr::mutate(
      Rank = dplyr::row_number(),
      is_forced = .key %in% forced_keys,
      is_top = (Rank <= top_k_labels) | is_forced,
      label_GP = ifelse(is_top, as.character(GP), "")
    ) %>%
    dplyr::ungroup()

  if (has_loading && shape_by_loading) {
    df <- df %>%
      dplyr::mutate(
        is_low_loading = Loading < Overall_Loading,
        plot_shape = dplyr::case_when(
          is_top & !is_low_loading ~ 21L, # hollow circle
          is_top & is_low_loading ~ 24L, # hollow triangle
          !is_top & !is_low_loading ~ 16L, # solid circle
          !is_top & is_low_loading ~ 17L # solid triangle
        )
      )
  } else if (has_loading && filter_positive) {
    df <- df %>%
      dplyr::mutate(
        is_low_loading = Loading < Overall_Loading,
        plot_shape = dplyr::case_when(
          is_top & is_low_loading ~ 24L, # hollow triangle (forced down-reg)
          is_top & !is_low_loading ~ 21L, # hollow circle (top up-reg)
          TRUE ~ 16L # solid background dot
        )
      )
  } else {
    df$plot_shape <- ifelse(df$is_top, 21L, 16L)
  }

  pos <- ggplot2::position_jitter(width = 0.2, height = 0, seed = seed)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Lineage, y = Value, color = Lineage))

  if (!is.null(threshold_line)) {
    p <- p +
      ggplot2::geom_hline(yintercept = threshold_line, linetype = "dashed", color = threshold_color, linewidth = 0.8)
  }

  p +
    ggplot2::geom_jitter(
      ggplot2::aes(
        shape = I(plot_shape),
        size = I(ifelse(is_top, 2.5, 1.5)),
        alpha = I(ifelse(is_top, 1, bg_alpha)),
        stroke = I(ifelse(is_top, 1.2, 0)),
        fill = I(ifelse(is_top, "white", "transparent"))
      ),
      position = pos
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label_GP),
      position = pos, size = 3.5, color = "black", fontface = "bold",
      box.padding = 0.5, point.padding = 0.3, min.segment.length = 0, segment.color = "grey50", max.overlaps = Inf
    ) +
    ggplot2::scale_color_manual(values = lineage_colors()) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Cell Lineage", y = y_label) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, face = if (italic_subtitle) "italic" else "plain"),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "none"
    )
}

# Colors an MDE embedding by a single GP's loading. Background cells
# (loading below `lower_bound`) are drawn faint/grey; the color scale is
# calibrated on active cells only, clipped at the `clip_q` quantile.
plot_loadings_on_mde <- function(
  mde,
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
  stopifnot(is.matrix(mde) || is.data.frame(mde))
  df <- as.data.frame(mde)
  stopifnot(all(c("MDE_1", "MDE_2") %in% colnames(df)))

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

  ggplot2::ggplot() +
    {
      if (show_bg) {
        scattermore::geom_scattermore(
          data = df, mapping = ggplot2::aes(x = MDE_1, y = MDE_2),
          pointsize = size, color = bg_color, alpha = bg_alpha
        )
      }
    } +
    scattermore::geom_scattermore(
      data = df[!is.na(df$loading_plot), ],
      mapping = ggplot2::aes(x = MDE_1, y = MDE_2, color = loading_plot),
      pointsize = size * 2, alpha = 1, na.rm = TRUE
    ) +
    ggplot2::scale_color_gradient(
      low = low_color, high = high_color, limits = c(0, q_hi),
      oob = scales::squish, breaks = scales::pretty_breaks(4), name = paste0("GP", factor_num)
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = paste0("Factor ", factor_num), x = "MDE 1", y = "MDE 2") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
}
