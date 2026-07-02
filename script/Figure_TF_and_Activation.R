# BiocManager::install("org.Mm.eg.db")

library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(scales)
library(pheatmap)
library(RColorBrewer)

mm <- org.Mm.eg.db
go2eg <- as.list(org.Mm.egGO2ALLEGS)
symbols <- AnnotationDbi::select(
  mm,
  keys = unique(unlist(go2eg)),
  columns = "SYMBOL",
  keytype = "ENTREZID"
)
# DNA-binding TF activity (GO:0003700) plus the Tox family.
tf <- c(
  sort(symbols$SYMBOL[symbols$ENTREZID %in% unique(go2eg[["GO:0003700"]])]),
  "Tox",
  "Tox2",
  "Tox3",
  "Tox4"
) %>%
  sort() %>%
  unique()


#####################################################
#####################################################
#####################################################
### Defining directory and loading functions
#####################################################
#####################################################
data_path <- "data/"
code_path <- "code/"
figure_path <- "figures/Figure_ActivationTF/"
dir.create(figure_path, showWarnings = FALSE, recursive = TRUE)
source(paste0(code_path, "filtering_membership.R")) # for scale_cols()


#####################################################
#####################################################
#####################################################
### Loading data and normalizing F
#####################################################
#####################################################
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F_pm_filtered) <- paste0("GP", seq_len(ncol(F_pm_filtered)))
# Per-GP normalization so each column has max(|F|) = 1, matching Figure_TF.R.
Df <- apply(F_pm_filtered, 2, function(x) max(abs(x), na.rm = TRUE))
F_pm_norm_col <- scale_cols(F_pm_filtered, 1 / Df)
colnames(F_pm_norm_col) <- paste0("GP", seq_len(ncol(F_pm_norm_col)))


#####################################################
#####################################################
#####################################################
### 21 highlighted GPs from Figure_Activation.R
### (GPs_of_interest 17 + GPs_extra_z 4)
#####################################################
#####################################################
highlighted_GPs <- paste0(
  "GP",
  c(
    25, 26, 10, 43, 12, 58, 171, 9, 79, 35, 72, 177,
    161, 152, 162, 36, 56,        # GPs_of_interest
    181, 32, 80, 163               # GPs_extra_z
  )
)
F_sub <- F_pm_norm_col[, highlighted_GPs, drop = FALSE]

# Activation-figure colors (Figure_Activation.R highlight_colors_z, with
# GP35 overwritten to darkred and the 4 extras_z set to darkred).
activation_gp_colors <- c(
  GP56 = "blue", GP162 = "blue", GP36 = "blue", GP152 = "blue",
  GP161 = "blue", GP177 = "blue", GP72 = "blue", GP79 = "blue", GP12 = "blue",
  GP43 = "darkorange2", GP10 = "darkorange2", GP58 = "darkorange2",
  GP25 = "darkred", GP26 = "darkred", GP35 = "darkred",
  GP181 = "darkred", GP32 = "darkred", GP80 = "darkred", GP163 = "darkred",
  GP9 = "darkgreen", GP171 = "darkgreen"
)
# Order in which color groups appear top-to-bottom in the network.
gp_color_group_order <- c("darkred", "darkorange2", "darkgreen", "blue")


#####################################################
#####################################################
#####################################################
### Filter the TF list: keep TFs with max normalized score > 0.25
### in at least one of the 21 highlighted GPs.
#####################################################
#####################################################
tf_gp_threshold <- 0.25
tf_in_F <- intersect(tf, rownames(F_sub))
tf_max_score <- apply(F_sub[tf_in_F, , drop = FALSE], 1, max, na.rm = TRUE)
selected_tfs <- sort(names(tf_max_score)[tf_max_score > tf_gp_threshold])
message(sprintf(
  "TF filter: %d / %d TFs (from %d in F) clear score > %.2f in any of the %d GPs",
  length(selected_tfs),
  length(tf),
  length(tf_in_F),
  tf_gp_threshold,
  length(highlighted_GPs)
))


#####################################################
#####################################################
#####################################################
### TF-GP network plot (inlined from Figure_TF.R::plot_tf_gp_network_v2)
### Genes are shown as italic text next to each GP, not as nodes.
#####################################################
#####################################################
# Iterative barycenter ordering for a bipartite TF-GP graph.
# `adj`           : TF x GP weighted matrix.
# `gp_groups`     : named vector (color group label per GP) or NULL.
# `gp_group_order`: preferred top-to-bottom order of color groups (NULL = data order).
# GP order is constrained to keep each color group contiguous; within each
# group, both GPs and TFs are reordered by weighted barycenter until stable.
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
    if (length(missing_g)) {
      gp_groups[missing_g] <- "_unassigned"
    }
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

  # Seed: GPs grouped by color, then iterate.
  adj <- adj[, reorder_gps_within_groups(adj), drop = FALSE]
  for (i in seq_len(n_iter)) {
    new_tf <- reorder_tfs(adj)
    adj <- adj[new_tf, , drop = FALSE]
    new_gp <- reorder_gps_within_groups(adj)
    if (identical(new_gp, colnames(adj)) &&
        identical(new_tf, rownames(adj))) {
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
  # Initial GP order: numeric (will be overwritten by optimization).
  connected_gps <- connected_gps[order(as.numeric(sub(
    "GP", "", connected_gps
  )))]
  # Drop TFs that ended up with no edges.
  selected_tfs <- intersect(selected_tfs, unique(tf_gp_edges$from))

  # Weighted adjacency matrix for layout optimization.
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

  # Top of plot -> first GP in optimized order.
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
        x = x0, y = y0, xend = x1, yend = y1,
        color = from, linewidth = weight
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


#####################################################
#####################################################
#####################################################
### Build and save the network plot
#####################################################
#####################################################
# F_sub is restricted to the 21 highlighted GPs, so the network is naturally
# scoped — `connected_gps` inside the function can only land on these 21.
network_plot <- plot_tf_gp_network_v2(
  F                = F_sub,
  selected_tfs     = selected_tfs,
  tf_gp_threshold  = tf_gp_threshold,
  top_genes_per_gp = 5,
  gp_colors        = activation_gp_colors,
  gp_group_order   = gp_color_group_order,
  optimize_layout  = TRUE,
  barycenter_iter  = 12,
  gp_spacing       = 1.5,
  node_size_tf     = 6,
  node_size_gp     = 5,
  label_size_tf    = 4.5,
  label_size_gp    = 4,
  label_size_gene  = 3.4
)

# Size the page so many TFs / GPs do not get crushed.
plot_height <- min(
  60,
  max(12, length(selected_tfs) * 0.35, length(highlighted_GPs) * 1.5 * 0.55 + 2)
)
ggsave(
  filename = file.path(figure_path, "TF_GP_network_activation21.pdf"),
  plot = network_plot,
  width = 18,
  height = plot_height,
  limitsize = FALSE
)


#####################################################
#####################################################
#####################################################
### Companion heatmap of the same TF x 21-GP submatrix
#####################################################
#####################################################
F_tf_sub <- F_sub[selected_tfs, , drop = FALSE]
pheatmap(
  F_tf_sub,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(
    100
  ),
  breaks = seq(-1, 1, length.out = 101),
  border_color = "grey80",
  fontsize_row = 6,
  fontsize_col = 9,
  main = sprintf(
    "TF (n=%d) x 21 activation GPs (normalized F, score > %.2f filter)",
    length(selected_tfs), tf_gp_threshold
  ),
  filename = paste0(figure_path, "TF_GP_heatmap_activation21.pdf"),
  width = 8,
  height = max(8, length(selected_tfs) * 0.15)
)
