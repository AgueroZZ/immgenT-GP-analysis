# Bipartite TF-GP network layout + plotting.
# Used by script-refactor/Figure3.R (panel c) and FigureS3.R are both built
# from this shared machinery, ported from script/Figure_Activation.R (which
# itself absorbed this logic from script/Figure_TF.R and
# script/Figure_TF_and_Activation.R -- all three files had their own
# near-identical copy before this refactor).

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

# Draws a bipartite TF (left) -> GP (right) network from a normalized gene
# score matrix F (genes x GPs). Edges connect a TF to a GP when F[tf, gp]
# exceeds tf_gp_threshold. Each GP node is annotated with its top
# `top_genes_per_gp` positively-scoring genes.
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
    dplyr::left_join(tf_nodes %>% dplyr::rename(x0 = x, y0 = y), by = c("from" = "node")) %>%
    dplyr::left_join(gp_nodes %>% dplyr::rename(x1 = x, y1 = y), by = c("to" = "node")) %>%
    dplyr::filter(is.finite(x0), is.finite(y0), is.finite(x1), is.finite(y1))

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

  ggplot2::ggplot() +
    ggplot2::geom_curve(
      data = edge_df,
      ggplot2::aes(
        x = x0,
        y = y0,
        xend = x1,
        yend = y1,
        color = from,
        linewidth = weight
      ),
      curvature = edge_curvature,
      alpha = edge_alpha,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.07, "inches"), type = "closed")
    ) +
    ggplot2::geom_point(
      data = gp_nodes,
      ggplot2::aes(x = x, y = y, fill = fill_color),
      shape = 22,
      color = "black",
      size = node_size_gp
    ) +
    ggplot2::geom_point(
      data = tf_nodes,
      ggplot2::aes(x = x, y = y, fill = fill_color),
      shape = 21,
      color = "black",
      size = node_size_tf,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_text(
      data = tf_nodes,
      ggplot2::aes(x = x, y = y, label = node),
      hjust = 1.3,
      size = label_size_tf,
      fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = gp_nodes,
      ggplot2::aes(x = x, y = y, label = node),
      hjust = 1.3,
      size = label_size_gp,
      color = "grey25"
    ) +
    ggplot2::geom_text(
      data = gene_label_df,
      ggplot2::aes(x = x, y = y, label = label),
      hjust = 0,
      size = label_size_gene,
      color = "grey40",
      fontface = "italic"
    ) +
    ggplot2::scale_color_manual(values = tf_colors, name = "TF") +
    ggplot2::scale_linewidth_continuous(range = c(0.3, 1.8), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(-0.25, 4.5), clip = "off") +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(20, 20, 20, 80)
    )
}
