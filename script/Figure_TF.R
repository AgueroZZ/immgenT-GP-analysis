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
figure_path <- "figures/Figure5_TF/"
source(paste0(code_path, "filtering_membership.R"))


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
seurat_meta <- readRDS(paste0(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
# protein_mat_normalized_lognorm <- readRDS(paste0(data_path, "protein_mat_normalized_lognorm.rds"))
# protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[rownames(L_pm_filtered) , "CD44"]
# umap_result <- readRDS(paste0(data_path, "umap_result.rds"))
# colnames(umap_result) <- c("UMAP_1", "UMAP_2")
# umap_result <- umap_result[rownames(L_pm_filtered),]
# df_umap <- as.data.frame(umap_result)
CD4_TF_list <- read.table(
  paste0(data_path, "/TF-list/CD4_tf_level2_AUC.txt"),
  header = TRUE
)
CD8_TF_list <- read.table(
  paste0(data_path, "/TF-list/CD8_tf_level2_AUC.txt"),
  header = TRUE
)
CD8aa_TF_list <- read.table(
  paste0(data_path, "/TF-list/CD8aa_tf_level2_AUC.txt"),
  header = TRUE
)
gdT_TF_list <- read.table(
  paste0(data_path, "/TF-list/gdT_tf_level2_AUC.txt"),
  header = TRUE
)
DN_TF_list <- read.table(
  paste0(data_path, "/TF-list/DN_tf_level2_AUC.txt"),
  header = TRUE
)
DP_TF_list <- read.table(
  paste0(data_path, "/TF-list/DP_tf_level2_AUC.txt"),
  header = TRUE
)
Treg_TF_list <- read.table(
  paste0(data_path, "/TF-list/Treg_tf_level2_AUC.txt"),
  header = TRUE
)
non_conv_TF_list <- read.table(
  paste0(data_path, "/TF-list/nonconv_tf_level2_AUC.txt"),
  header = TRUE
)


#####################################################
#####################################################
#####################################################
### Defining additional functions
#####################################################
#####################################################
compute_tf_rank_matrix <- function(
  F,
  TFs,
  mode = c("positive", "abs", "raw"),
  ties.method = "average",
  verbose = TRUE
) {
  mode <- match.arg(mode)

  if (is.null(rownames(F))) {
    stop("F must have rownames (gene names).")
  }
  if (is.null(colnames(F))) {
    colnames(F) <- paste0("GP", seq_len(ncol(F)))
  }

  TFs_in <- intersect(TFs, rownames(F))
  TFs_drop <- setdiff(TFs, TFs_in)

  if (verbose && length(TFs_drop) > 0) {
    message(sprintf(
      "Dropping %d TFs not found in rownames(F). Example: %s",
      length(TFs_drop),
      paste(head(TFs_drop, 5), collapse = ", ")
    ))
  }
  if (length(TFs_in) == 0) {
    stop("None of the requested TFs are present in rownames(F).")
  }

  score_mat <- switch(
    mode,
    positive = pmax(F, 0),
    abs = abs(F),
    raw = F
  )

  # rank within each column (GP): higher score => smaller rank (1 = best)
  rank_all <- apply(score_mat, 2, function(x) {
    rank(-x, ties.method = ties.method, na.last = "keep")
  })

  # subset rows to TFs
  rank_tf <- rank_all[TFs_in, , drop = FALSE]

  list(
    rank = rank_tf, # matrix: TF x GP
    TFs_used = TFs_in,
    TFs_dropped = TFs_drop
  )
}
.compute_tf_rank_score <- function(F, tf, mode = c("positive", "abs", "raw")) {
  mode <- match.arg(mode)

  if (is.null(rownames(F))) {
    stop("F must have rownames (gene names).")
  }
  if (!tf %in% rownames(F)) {
    stop(sprintf("TF '%s' not found in rownames(F).", tf))
  }

  G <- ncol(F)
  cols <- colnames(F)
  if (is.null(cols)) {
    cols <- paste0("GP", seq_len(G))
  }

  score_mat <- switch(
    mode,
    positive = pmax(F, 0),
    abs = abs(F),
    raw = F
  )

  # Rank within each column: higher score -> smaller rank (rank=1 is best)
  rank_mat <- apply(score_mat, 2, function(x) {
    rank(-x, ties.method = "average", na.last = "keep")
  })
  tf_rank <- as.numeric(rank_mat[tf, ])

  data.frame(
    GEP = factor(cols, levels = cols),
    rank = tf_rank,
    tf = tf,
    stringsAsFactors = FALSE
  )
}
plot_tf_rank_bars <- function(
  F,
  TFs,
  mode = c("positive", "abs", "raw"),
  rank_max = 100,
  custom_x_breaks = NULL,
  x_label_every = 10,
  facet = TRUE,
  alpha = 0.9,
  n_label_best = 1,
  # 3. New: Take a custom color vector (named or unnamed)
  tf_colors = NULL
) {
  mode <- match.arg(mode)

  # .compute_tf_rank_score helper assumed to exist
  df <- do.call(
    rbind,
    lapply(TFs, function(tf) .compute_tf_rank_score(F, tf, mode = mode))
  )

  # Ensure numeric rank; drop missing ranks
  df$rank <- as.numeric(df$rank)
  df <- df[is.finite(df$rank), , drop = FALSE]

  # Update column name to GP
  colnames(df)[colnames(df) == "GEP"] <- "GP"
  levs <- levels(df$GP)
  df$x <- as.integer(df$GP)

  # "Skyscraper" logic
  df$rank_clip <- pmin(df$rank, rank_max)
  df$height <- rank_max - df$rank_clip

  # --- X-axis Handling ---
  if (!is.null(custom_x_breaks)) {
    x_breaks <- which(levs %in% custom_x_breaks)
    x_labels <- levs[x_breaks]
  } else {
    keep_idx <- seq(1, length(levs), by = x_label_every)
    x_breaks <- keep_idx
    x_labels <- levs[keep_idx]
  }

  # --- Y-axis Handling ---
  y_breaks <- seq(0, rank_max, by = 25)
  y_labels <- rank_max - y_breaks
  y_labels[y_labels == 0] <- 1

  p <- ggplot(df, aes(x = x, y = height, fill = tf)) +
    geom_col(alpha = alpha, width = 0.9) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = c(0, rank_max * 1.2)
    ) +
    labs(
      title = sprintf("TF rank skyline across GPs (top %d)", rank_max),
      x = "GP",
      y = "Rank within each GP (1 = highest)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = if (facet) "none" else "right",
      strip.text = element_text(face = "bold", size = 14, color = "black"),
      strip.background = element_rect(
        fill = "grey95",
        color = "black",
        linewidth = 0.5
      ),
      panel.spacing = unit(1, "lines")
    )

  # --- Apply Custom Colors if provided ---
  if (!is.null(tf_colors)) {
    p <- p + scale_fill_manual(values = tf_colors)
  }

  if (facet) {
    p <- p + facet_wrap(~tf, ncol = 1, scales = "fixed")
  }

  # --- Labeling logic with ggrepel ---
  if (n_label_best > 0 && nrow(df) > 0) {
    best_df <- df %>%
      filter(rank <= rank_max) %>%
      group_by(tf) %>%
      slice_min(order_by = rank, n = n_label_best, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(label = as.character(GP))

    p <- p +
      ggrepel::geom_text_repel(
        data = best_df,
        aes(x = x, y = height, label = label),
        inherit.aes = FALSE,
        vjust = -0.5,
        size = 3.5,
        fontface = "bold",
        box.padding = 0.3,
        point.padding = 0.2,
        direction = "y",
        segment.color = "grey50"
      )
  }

  return(p)
}
plot_tf_rank_bars_ordered <- function(
  F,
  TFs,
  mode = c("positive", "abs", "raw"),
  rank_max = 100,
  gp_order = NULL, # New: Vector of GP names in desired order
  x_label_every = 10,
  facet = TRUE,
  alpha = 0.9,
  n_label_best = 1,
  tf_colors = NULL
) {
  mode <- match.arg(mode)

  # Compute ranks (using your existing helper)
  df <- do.call(
    rbind,
    lapply(TFs, function(tf) .compute_tf_rank_score(F, tf, mode = mode))
  )
  df$rank <- as.numeric(df$rank)
  df <- df[is.finite(df$rank), , drop = FALSE]
  colnames(df)[colnames(df) == "GEP"] <- "GP"

  # --- HANDLE ORDERING ---
  if (!is.null(gp_order)) {
    # Filter df to only include GPs in the order vector (if necessary)
    df <- df[df$GP %in% gp_order, ]
    # Force the factor levels to match the heatmap order
    df$GP <- factor(df$GP, levels = gp_order)
  } else {
    # Default behavior if no order is provided
    df$GP <- factor(df$GP, levels = paste0("GP", seq_len(ncol(F))))
  }

  # Use the factor level as the x coordinate
  df$x <- as.numeric(df$GP)

  # "Skyscraper" logic
  df$rank_clip <- pmin(df$rank, rank_max)
  df$height <- rank_max - df$rank_clip

  # Define X-axis breaks based on the NEW order
  # Label every N-th GP in the clustered sequence
  total_gps <- length(unique(df$GP))
  x_breaks <- seq(1, total_gps, by = x_label_every)
  x_labels <- levels(df$GP)[x_breaks]

  # Y-axis Handling
  y_breaks <- seq(0, rank_max, by = 25)
  y_labels <- rank_max - y_breaks
  y_labels[y_labels == 0] <- 1

  p <- ggplot(df, aes(x = x, y = height, fill = tf)) +
    geom_col(alpha = alpha, width = 0.9) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = c(0, rank_max * 1.2)
    ) +
    labs(
      title = "TF rank skyline (Clustered by Heatmap Order)",
      x = "Gene Programs (Ordered by Similarity)",
      y = "Rank (1 = highest)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      legend.position = if (facet) "none" else "right",
      strip.text = element_text(face = "bold")
    )

  if (!is.null(tf_colors)) {
    p <- p + scale_fill_manual(values = tf_colors)
  }
  if (facet) {
    p <- p + facet_wrap(~tf, ncol = 1, scales = "fixed")
  }

  # Add labels for top GPs
  if (n_label_best > 0) {
    best_df <- df %>%
      filter(rank <= rank_max) %>%
      group_by(tf) %>%
      slice_min(order_by = rank, n = n_label_best, with_ties = FALSE) %>%
      ungroup()

    p <- p +
      ggrepel::geom_text_repel(
        data = best_df,
        aes(x = x, y = height, label = GP),
        inherit.aes = FALSE,
        size = 3,
        fontface = "bold"
      )
  }

  return(p)
}
plot_tf_rank_bars_free_y <- function(
  F,
  TFs,
  mode = c("positive", "abs", "raw"),
  rank_max = 100, # Depth of the window shown (e.g., 200)
  gp_order = NULL, # Vector of GP names in clustered order
  x_label_every = 1,
  facet = TRUE,
  alpha = 0.9,
  n_label_best = 1,
  tf_colors = NULL,
  free_y = FALSE, # Adapt y-axis per TF
  label_internal = FALSE # Move TF name inside plot area
) {
  mode <- match.arg(mode)

  # 1. Compute ranks
  df <- do.call(
    rbind,
    lapply(TFs, function(tf) .compute_tf_rank_score(F, tf, mode = mode))
  )
  df$rank <- as.numeric(df$rank)
  df <- df[is.finite(df$rank), , drop = FALSE]
  colnames(df)[colnames(df) == "GEP"] <- "GP"

  # 2. Handle Ordering (Crucial: Keep as Factor for Discrete X-axis)
  if (!is.null(gp_order)) {
    # Keep only the GPs that exist in the clustered order provided
    df <- df[df$GP %in% gp_order, ]
    df$GP <- factor(df$GP, levels = gp_order)
  } else {
    df$GP <- factor(df$GP, levels = paste0("GP", seq_len(ncol(F))))
  }

  # 3. Logic for "Skyscraper" Hanging from Top
  if (free_y) {
    df <- df %>%
      group_by(tf) %>%
      mutate(
        true_min = min(rank, na.rm = TRUE),
        # Height is rank_max at the best rank, 0 at the window floor
        height = rank_max - (rank - true_min)
      ) %>%
      # if height < 0, force it to be 0 (floor of the window)
      mutate(height = pmax(height, 0)) %>%
      ungroup()
  } else {
    # Global mode: ceiling is 1, floor is rank_max
    df$height <- rank_max - pmin(df$rank, rank_max)
  }

  # 4. Axis Setup
  all_gps <- levels(df$GP)
  x_breaks <- all_gps[seq(1, length(all_gps), by = x_label_every)]

  # 5. Build Plot (Using x = GP for discrete spacing)
  p <- ggplot(df, aes(x = GP, y = height, fill = tf)) +
    # --- ADDED: Horizontal line at the bottom (y=0) ---
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +

    geom_col(alpha = alpha, width = 0.9) +
    scale_x_discrete(breaks = x_breaks) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      legend.position = "none"
    )

  # 6. Y-Axis Labelling
  if (free_y) {
    p <- p +
      scale_y_continuous(
        breaks = seq(0, rank_max, by = 25),
        labels = function(h) paste0("+", rank_max - h)
      ) +
      labs(y = "Rank (Relative to TF Peak)")
  } else {
    y_breaks <- seq(0, rank_max, by = 25)
    p <- p +
      scale_y_continuous(
        breaks = y_breaks,
        labels = rev(y_breaks + 1)
      ) +
      labs(y = "Rank (1 = highest)")
  }

  if (!is.null(tf_colors)) {
    p <- p + scale_fill_manual(values = tf_colors)
  }

  # 7. Faceting and Internal Labels
  facet_scales <- if (free_y) "free_y" else "fixed"
  if (facet) {
    p <- p + facet_wrap(~tf, ncol = 1, scales = facet_scales)
    if (label_internal) {
      p <- p +
        theme(
          strip.background = element_blank(),
          strip.text = element_blank(),
          # Removing axis line since we added the manual hline for better control per facet
          axis.line.x = element_blank()
        ) +
        geom_text(
          data = data.frame(tf = TFs, label = TFs),
          aes(x = -Inf, y = Inf, label = label),
          inherit.aes = FALSE,
          hjust = -0.1,
          vjust = 1.2,
          fontface = "bold",
          size = 5
        )
    }
  }

  # 8. Best GP Labeling (Annotates with Absolute Rank)
  if (n_label_best > 0) {
    best_df <- df %>%
      group_by(tf) %>%
      slice_min(order_by = rank, n = n_label_best, with_ties = FALSE) %>%
      ungroup()

    p <- p +
      ggrepel::geom_text_repel(
        data = best_df,
        aes(
          x = GP,
          y = height,
          label = paste0(GP, "\n(rank ", round(rank), ")")
        ),
        inherit.aes = FALSE,
        size = 2.5,
        fontface = "bold",
        lineheight = 0.8
      )
  }

  return(p)
}


#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
### Network helpers for TF-GP-gene visualization
#####################################################
#####################################################
validate_tf_network_edges <- function(
  F,
  edge_csv = NULL,
  selected_tfs = NULL,
  tf_gp_threshold = 0.25,
  tol = 1e-6
) {
  if (is.null(edge_csv) || !file.exists(edge_csv)) {
    stop("edge_csv must be a valid file path.")
  }

  edge_df <- read.csv(edge_csv, stringsAsFactors = FALSE)
  if (!all(c("from", "to") %in% colnames(edge_df))) {
    stop("edge_csv must contain columns: from, to.")
  }

  if (!is.null(selected_tfs)) {
    edge_df <- edge_df[edge_df$from %in% selected_tfs, , drop = FALSE]
  }

  if (is.null(colnames(F))) {
    colnames(F) <- paste0("GP", seq_len(ncol(F)))
  }

  edge_df$to <- gsub("^K", "GP", edge_df$to)
  edge_df <- edge_df[
    edge_df$to %in% colnames(F) & edge_df$from %in% rownames(F),
  ]

  edge_keys <- paste(edge_df$from, edge_df$to, sep = "->")
  tfs_in <- unique(edge_df$from)

  candidate_list <- lapply(tfs_in, function(tf) {
    scores <- as.numeric(F[tf, ])
    names(scores) <- colnames(F)
    sel <- names(scores)[is.finite(scores) & scores > tf_gp_threshold]
    if (length(sel) == 0) {
      return(NULL)
    }
    data.frame(
      from = tf,
      to = sel,
      weight = scores[sel],
      stringsAsFactors = FALSE
    )
  })

  candidate_edges <- do.call(rbind, candidate_list)
  if (is.null(candidate_edges)) {
    candidate_edges <- data.frame(
      from = character(0),
      to = character(0),
      weight = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  candidate_keys <- paste(candidate_edges$from, candidate_edges$to, sep = "->")

  missing_in_file <- setdiff(candidate_keys, edge_keys)
  extra_in_file <- setdiff(edge_keys, candidate_keys)
  matched <- setdiff(edge_keys, extra_in_file)

  merged <- merge(
    edge_df,
    candidate_edges,
    by = c("from", "to"),
    suffixes = c("_file", "_candidate")
  )
  if (nrow(merged) > 0 && "weight" %in% colnames(merged)) {
    merged$weight_file <- as.numeric(merged$weight_file)
    merged$weight_candidate <- as.numeric(merged$weight_candidate)
    merged$weight_diff <- abs(merged$weight_file - merged$weight_candidate)
    max_diff <- max(merged$weight_diff, na.rm = TRUE)
  } else {
    max_diff <- NA_real_
  }

  list(
    expected_edges = nrow(edge_df),
    candidate_edges = nrow(candidate_edges),
    matched_edges = length(matched),
    missing_in_file = missing_in_file,
    extra_in_file = extra_in_file,
    max_weight_diff = max_diff,
    is_consistent = (length(missing_in_file) == 0 &&
      length(extra_in_file) == 0 &&
      (is.na(max_diff) || max_diff <= tol))
  )
}

build_tf_gp_gene_network <- function(
  F,
  selected_tfs,
  tf_gp_threshold = 0.25,
  top_genes_per_gp = 3,
  gene_threshold = 0.25
) {
  if (is.null(rownames(F))) {
    stop("F must have rownames.")
  }
  if (is.null(colnames(F))) {
    stop("F must have column names.")
  }

  selected_tfs <- intersect(selected_tfs, rownames(F))
  if (length(selected_tfs) == 0) {
    stop("No requested TFs were found in F rownames.")
  }

  tf_gp_edges <- do.call(
    rbind,
    lapply(selected_tfs, function(tf) {
      vals <- as.numeric(F[tf, ])
      names(vals) <- colnames(F)
      keep <- which(is.finite(vals) & vals > tf_gp_threshold)
      if (length(keep) == 0) {
        return(NULL)
      }
      data.frame(
        from = tf,
        to = names(vals)[keep],
        weight = vals[keep],
        edge_type = "TF-GP",
        stringsAsFactors = FALSE
      )
    })
  )
  if (is.null(tf_gp_edges)) {
    tf_gp_edges <- data.frame(
      from = character(0),
      to = character(0),
      weight = numeric(0),
      edge_type = character(0),
      stringsAsFactors = FALSE
    )
  }

  connected_gps <- unique(tf_gp_edges$to)
  gp_gene_edges <- do.call(
    rbind,
    lapply(connected_gps, function(gp) {
      vals <- as.numeric(F[, gp])
      names(vals) <- rownames(F)
      keep <- which(is.finite(vals) & vals > gene_threshold)
      if (length(keep) == 0) {
        return(NULL)
      }
      ranked_genes <- keep[order(vals[keep], decreasing = TRUE)]
      top_k <- ranked_genes[1:min(top_genes_per_gp, length(ranked_genes))]
      data.frame(
        from = gp,
        gene = names(vals)[top_k],
        weight = vals[top_k],
        edge_type = "GP-GENE",
        stringsAsFactors = FALSE
      )
    })
  )
  if (is.null(gp_gene_edges)) {
    gp_gene_edges <- data.frame(
      from = character(0),
      gene = character(0),
      weight = numeric(0),
      edge_type = character(0),
      stringsAsFactors = FALSE
    )
  }
  gp_gene_edges$to <- if (nrow(gp_gene_edges) > 0) {
    gp_gene_edges$gene
  } else {
    character(0)
  }
  gp_gene_edges$node <- if (nrow(gp_gene_edges) > 0) {
    paste0(gp_gene_edges$from, ":", gp_gene_edges$gene)
  } else {
    character(0)
  }

  # Build nodes in three layers.
  tf_nodes <- data.frame(
    node = unique(tf_gp_edges$from),
    type = "TF",
    label = unique(tf_gp_edges$from),
    x = 1,
    y = seq_len(length(unique(tf_gp_edges$from))),
    stringsAsFactors = FALSE
  )

  gp_nodes <- data.frame(
    node = connected_gps,
    type = "GP",
    label = connected_gps,
    x = 2,
    y = seq_len(length(connected_gps)),
    stringsAsFactors = FALSE
  )

  gene_node_df <- gp_gene_edges %>%
    group_by(from) %>%
    arrange(desc(weight), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  if (nrow(gene_node_df) > 0) {
    gp_y <- setNames(gp_nodes$y, gp_nodes$node)
    gene_node_df <- gene_node_df %>%
      group_by(from) %>%
      mutate(
        n_in_gp = n(),
        rank_centered = rank - (n_in_gp + 1) / 2,
        y = gp_y[from] + rank_centered * 0.15
      ) %>%
      ungroup()
  }

  gene_nodes <- if (nrow(gene_node_df) > 0) {
    data.frame(
      node = gene_node_df$node,
      type = "Gene",
      label = gene_node_df$gene,
      x = 3,
      y = gene_node_df$y,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      node = character(0),
      type = character(0),
      label = character(0),
      x = numeric(0),
      y = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  list(
    tf_gp_edges = tf_gp_edges,
    gp_gene_edges = gp_gene_edges,
    tf_nodes = tf_nodes,
    gp_nodes = gp_nodes,
    gene_nodes = gene_nodes,
    selected_tfs = selected_tfs
  )
}

plot_tf_gp_gene_network <- function(
  F,
  selected_tfs,
  tf_gp_threshold = 0.25,
  top_genes_per_gp = 3,
  gene_threshold = 0.25,
  gp_order = NULL,
  node_size = 4,
  edge_arrow_size = 0.1,
  tf_spacing = 1.8,
  gp_spacing = 2.0,
  gene_spacing = 0.55,
  panel_ratio = 0.75
) {
  net <- build_tf_gp_gene_network(
    F = F,
    selected_tfs = selected_tfs,
    tf_gp_threshold = tf_gp_threshold,
    top_genes_per_gp = top_genes_per_gp,
    gene_threshold = gene_threshold
  )

  tf_nodes <- net$tf_nodes
  gp_nodes <- net$gp_nodes
  gene_nodes <- net$gene_nodes
  tf_gp_edges <- net$tf_gp_edges
  gp_gene_edges <- net$gp_gene_edges

  if (!is.null(gp_order) && length(gp_order) > 0) {
    gp_nodes <- gp_nodes[gp_nodes$node %in% gp_order, , drop = FALSE]
    gp_nodes <- gp_nodes[order(match(gp_nodes$node, gp_order), na.last = NA), ]
  }

  if (nrow(gp_nodes) == 0) {
    stop("No GP nodes were generated from selected TFs.")
  }

  gp_nodes <- gp_nodes %>%
    arrange(as.numeric(sub("GP", "", node))) %>%
    mutate(y = seq_len(n()) * gp_spacing)

  if (nrow(tf_nodes) > 0) {
    tf_center <- mean(range(gp_nodes$y))
    if (nrow(tf_nodes) == 1) {
      tf_nodes$y <- tf_center
    } else {
      tf_nodes <- tf_nodes %>%
        arrange(desc(node)) %>%
        mutate(y = tf_center + (seq_len(n()) - (n() + 1) / 2) * tf_spacing)
    }
  }

  if (nrow(gene_nodes) > 0) {
    gp_y <- setNames(gp_nodes$y, gp_nodes$node)
    gene_nodes <- gene_nodes %>%
      mutate(gp = sub(":.*$", "", node)) %>%
      group_by(gp) %>%
      mutate(
        n_in_gp = n(),
        rank = row_number(),
        y = gp_y[gp] + (rank - (n_in_gp + 1) / 2) * gene_spacing
      ) %>%
      ungroup() %>%
      select(-n_in_gp, -rank)
  }

  tf_nodes$x <- 1
  gp_nodes$x <- 2
  gene_nodes$x <- 3

  tf_gp_edges <- tf_gp_edges %>%
    filter(from %in% tf_nodes$node, to %in% gp_nodes$node)
  gp_gene_edges <- gp_gene_edges %>%
    filter(from %in% gp_nodes$node, node %in% gene_nodes$node)

  nodes <- bind_rows(tf_nodes, gp_nodes, gene_nodes)
  if (nrow(nodes) == 0) {
    stop("No nodes were generated for the network plot.")
  }
  coords <- nodes[, c("node", "x", "y")]

  tf_gp_edges_plot <- tf_gp_edges %>%
    left_join(
      coords %>% rename(from_x = x, from_y = y),
      by = c("from" = "node")
    ) %>%
    left_join(coords %>% rename(to_x = x, to_y = y), by = c("to" = "node")) %>%
    mutate(
      size = 0.5 + 1.0 * (weight / max(weight, na.rm = TRUE))
    ) %>%
    filter(
      is.finite(from_x),
      is.finite(from_y),
      is.finite(to_x),
      is.finite(to_y)
    )

  gp_gene_edges_plot <- if (nrow(gp_gene_edges) > 0) {
    gp_gene_edges %>%
      left_join(
        coords %>% rename(from_x = x, from_y = y),
        by = c("from" = "node")
      ) %>%
      left_join(
        coords %>% rename(to_x = x, to_y = y),
        by = c("node" = "node")
      ) %>%
      mutate(
        size = 0.4 + 0.8 * (weight / max(weight, na.rm = TRUE))
      ) %>%
      filter(
        is.finite(from_x),
        is.finite(from_y),
        is.finite(to_x),
        is.finite(to_y)
      )
  } else {
    data.frame(
      from_x = numeric(0),
      from_y = numeric(0),
      to_x = numeric(0),
      to_y = numeric(0),
      size = numeric(0),
      weight = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  p <- ggplot() +
    geom_segment(
      data = tf_gp_edges_plot,
      aes(
        x = from_x,
        y = from_y,
        xend = to_x,
        yend = to_y,
        size = size,
        color = weight,
        alpha = 0.9
      ),
      arrow = arrow(length = unit(edge_arrow_size, "inch"), type = "closed")
    ) +
    geom_segment(
      data = gp_gene_edges_plot,
      aes(
        x = from_x,
        y = from_y,
        xend = to_x,
        yend = to_y,
        size = size,
        color = weight
      ),
      lineend = "round"
    ) +
    geom_point(
      data = nodes,
      aes(x = x, y = y, fill = type, shape = type),
      color = "black",
      size = node_size
    ) +
    geom_text(
      data = nodes %>% filter(type != "Gene"),
      aes(x = x, y = y, label = label),
      hjust = -0.05,
      vjust = 0.4,
      size = 3.2,
      fontface = "bold"
    ) +
    ggrepel::geom_text_repel(
      data = nodes %>% filter(type == "Gene"),
      aes(x = x, y = y, label = label),
      direction = "y",
      nudge_x = 0.08,
      size = 2.8,
      box.padding = 0.2,
      point.padding = 0.1,
      segment.size = 0.15,
      segment.color = "grey80",
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c("TF" = "#1f77b4", "GP" = "#2ca02c", "Gene" = "#ff7f0e"),
      limits = c("TF", "GP", "Gene"),
      name = "Node type"
    ) +
    scale_shape_manual(
      values = c("TF" = 21, "GP" = 22, "Gene" = 24),
      name = "Node type"
    ) +
    scale_color_gradient(low = "#bdbdbd", high = "#1f78b4") +
    scale_size_identity() +
    guides(
      size = "none",
      alpha = "none",
      fill = guide_legend(),
      shape = guide_legend()
    ) +
    labs(
      title = "TF to GP to top genes network",
      x = NULL,
      y = NULL,
      color = "Weight"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    ) +
    coord_fixed(
      ratio = panel_ratio,
      xlim = c(0.4, 3.6),
      ylim = c(min(nodes$y) - 1, max(nodes$y) + 1),
      clip = "off"
    )

  return(list(
    plot = p,
    network = net,
    tf_gp_edges_plot = tf_gp_edges_plot,
    gp_gene_edges_plot = gp_gene_edges_plot,
    nodes = nodes
  ))
}

#####################################################
### Clean TF → GP network with inline gene labels
### (no gene nodes; genes shown as text next to each GP)
#####################################################
plot_tf_gp_network_v2 <- function(
  F,
  selected_tfs,
  tf_gp_threshold = 0.25,
  top_genes_per_gp = 5,
  tf_colors = NULL,
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

  # TF → GP edges above threshold
  tf_gp_edges <- do.call(
    rbind,
    lapply(selected_tfs, function(tf) {
      vals <- setNames(as.numeric(F[tf, ]), colnames(F))
      idx <- which(is.finite(vals) & vals > tf_gp_threshold)
      if (!length(idx)) {
        return(NULL)
      }
      data.frame(
        from = tf,
        to = names(vals)[idx],
        weight = vals[idx],
        stringsAsFactors = FALSE
      )
    })
  )
  if (is.null(tf_gp_edges)) {
    stop("No TF-GP edges found above threshold.")
  }

  # GPs in numeric order
  connected_gps <- unique(tf_gp_edges$to)
  connected_gps <- connected_gps[order(as.numeric(sub(
    "GP",
    "",
    connected_gps
  )))]
  n_gps <- length(connected_gps)
  n_tfs <- length(selected_tfs)

  # Top genes per GP (text only, no nodes)
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

  # Node positions: GPs on right (x=1), TFs on left (x=0)
  gp_y <- rev(seq_len(n_gps)) * gp_spacing
  gp_nodes <- data.frame(
    node = connected_gps,
    x = 1,
    y = gp_y,
    stringsAsFactors = FALSE
  )
  tf_nodes <- data.frame(
    node = selected_tfs,
    x = 0,
    y = seq(max(gp_y), min(gp_y), length.out = n_tfs),
    stringsAsFactors = FALSE
  )

  # Default qualitative colors for TFs
  if (is.null(tf_colors)) {
    tf_colors <- setNames(scales::hue_pal()(n_tfs), selected_tfs)
  }

  # Edge coordinates
  edge_df <- tf_gp_edges %>%
    left_join(tf_nodes %>% rename(x0 = x, y0 = y), by = c("from" = "node")) %>%
    left_join(gp_nodes %>% rename(x1 = x, y1 = y), by = c("to" = "node")) %>%
    filter(is.finite(x0), is.finite(y0), is.finite(x1), is.finite(y1))

  # Gene labels: comma-separated string placed right of GP node
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
    # Curved edges colored by TF, thickness ~ weight
    geom_curve(
      data = edge_df,
      aes(
        x = x0,
        y = y0,
        xend = x1,
        yend = y1,
        color = from,
        linewidth = weight
      ),
      curvature = edge_curvature,
      alpha = edge_alpha,
      arrow = arrow(length = unit(0.07, "inches"), type = "closed")
    ) +
    # GP squares
    geom_point(
      data = gp_nodes,
      aes(x = x, y = y),
      shape = 22,
      fill = "#4daf4a",
      color = "black",
      size = node_size_gp
    ) +
    # TF circles (colored)
    geom_point(
      data = tf_nodes,
      aes(x = x, y = y, fill = node),
      shape = 21,
      color = "black",
      size = node_size_tf,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = tf_colors) +
    # TF labels (left of node)
    geom_text(
      data = tf_nodes,
      aes(x = x, y = y, label = node),
      hjust = 1.3,
      size = label_size_tf,
      fontface = "bold"
    ) +
    # GP name (left of square)
    geom_text(
      data = gp_nodes,
      aes(x = x, y = y, label = node),
      hjust = 1.3,
      size = label_size_gp,
      color = "grey25"
    ) +
    # Top genes (right of GP square, italic)
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

### Heatmap of TF loadings across GPs
#####################################################
#####################################################
threshold <- 0 # include all TFs
TF_highlights <- c(
  CD4_TF_list$Gene[CD4_TF_list$auc >= threshold],
  CD8_TF_list$Gene[CD8_TF_list$auc >= threshold],
  CD8aa_TF_list$Gene[CD8aa_TF_list$auc >= threshold],
  gdT_TF_list$Gene[gdT_TF_list$auc >= threshold],
  DN_TF_list$Gene[DN_TF_list$auc >= threshold],
  DP_TF_list$Gene[DP_TF_list$auc >= threshold],
  non_conv_TF_list$Gene[non_conv_TF_list$auc >= threshold],
  Treg_TF_list$Gene[Treg_TF_list$auc >= threshold]
)
TF_highlights <- unique(TF_highlights)
# find abs max in each column of F_pm_filtered
Df <- apply(F_pm_filtered, 2, function(x) max(abs(x), na.rm = TRUE))
F_pm_norm_col <- scale_cols(F_pm_filtered, 1 / Df)
colnames(F_pm_norm_col) <- paste0("GP", 1:ncol(F_pm_norm_col))
F_pm_TF <- F_pm_norm_col[rownames(F_pm_norm_col) %in% TF_highlights, ]
# only include TFs with at least one factor > 0.4
keep <- apply((F_pm_TF), 1, function(v) any(v > 0.4, na.rm = TRUE))
F_pm_TF <- F_pm_TF[keep, , drop = FALSE]
# only include GPs with at least one TF factor > 0.4
keep_gp <- apply((F_pm_TF), 2, function(v) any(v > 0.4, na.rm = TRUE))
F_pm_TF <- F_pm_TF[, keep_gp, drop = FALSE]
pdf(paste0(figure_path, "TF_GP_heatmap.pdf"), width = 12, height = 10)
pheatmap(
  (F_pm_TF),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(
    100
  ),
  breaks = seq(-1, 1, length.out = 101),
  # add grey grid lines
  border_color = "grey10",
  fontsize_row = 6,
  fontsize_col = 8,
  main = "TF × GP loading heatmap (normalized; TFs only)",
  legend = TRUE
)
dev.off()


#####################################################
#####################################################
#####################################################
### Heatmap of TF ranks across GPs
#####################################################
#####################################################
rk_obj <- compute_tf_rank_matrix(F_pm_norm_col, TF_highlights, mode = "pos")
rank_tf <- rk_obj$rank
rank_max <- 100
rank_clip <- pmin(rank_tf, rank_max)
keep <- apply(rank_clip, 1, function(v) any(v < rank_max, na.rm = TRUE))
rank_clip <- rank_clip[keep, , drop = FALSE]
pal <- colorRampPalette(c("red", "white"))(10)
bk <- c(1, seq(10, 100, by = 10))
rank_filter <- 20
# keep only GPs that have at least one TF with rank < rank_filter
keep_gp <- apply(rank_clip, 2, function(v) any(v < rank_filter, na.rm = TRUE))
rank_clip <- rank_clip[, keep_gp, drop = FALSE]
# keep only TFs that have at least one GP with rank < rank_filter
keep_tf <- apply(rank_clip, 1, function(v) any(v < rank_filter, na.rm = TRUE))
rank_clip <- rank_clip[keep_tf, , drop = FALSE]
pheatmap(
  mat = rank_clip,
  color = pal,
  breaks = bk,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontface_row = "bold",
  fontsize_col = 7,
  main = sprintf("TF × GP rank heatmap (pos; ranks 1–%d)", rank_max),
  border_color = "grey80", # Use a lighter grey for better visibility
  legend_breaks = seq(1, rank_max, by = 10),
  legend_labels = as.character(seq(1, rank_max, by = 10)),
  filename = paste0(figure_path, "TF_GP_rank_heatmap.pdf"),
  width = 14,
  height = 12
)


#####################################################
#####################################################
#####################################################
### Skyscraper plot of TF ranks across GPs
#####################################################
#####################################################
x_breaks <- paste0("GP", seq(0, ncol(F_pm_norm_col), by = 10))
x_breaks[1] <- "GP1" # Ensure the first GP is labeled as GP1
p1 <- plot_tf_rank_bars(
  F = F_pm_norm_col,
  TFs = c("Tcf7", "Foxp3", "Sox13"),
  mode = "pos",
  rank_max = 100,
  custom_x_breaks = x_breaks,
  facet = TRUE,
  n_label_best = 5,
  tf_colors = c("#1b9e57", "#1b9e57", "#1b9e57")
)
p2 <- plot_tf_rank_bars(
  F = F_pm_norm_col,
  TFs = c("Runx3", "Gata3", "Foxn3"),
  mode = "pos",
  rank_max = 100,
  custom_x_breaks = x_breaks,
  facet = TRUE,
  n_label_best = 5,
  tf_colors = c("#1b9e57", "#1b9e57", "#1b9e57")
)
ggsave(
  p1,
  filename = paste0(figure_path, "TF_rank_skyline_pos.pdf"),
  width = 8,
  height = 10
)
ggsave(
  p2,
  filename = paste0(figure_path, "TF_rank_skyline_pos_2.pdf"),
  width = 8,
  height = 10
)


#####################################################
#####################################################
#####################################################
### Skyscraper plot of TF ranks across GPs (ordered)
#####################################################
#####################################################
JJ <- pheatmap(
  mat = rank_clip,
  color = pal,
  breaks = bk,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)
gp_clustered_order <- colnames(rank_clip)[JJ$tree_col$order]
p_clustered <- plot_tf_rank_bars_free_y(
  F = F_pm_norm_col,
  TFs = c(
    "Foxp3",
    "Zbtb16",
    "Sox13",
    "Rorc",
    "Pparg",
    "Hes1",
    "Gata3",
    "Tcf7",
    "Tbx21",
    "Tox",
    "Tox2"
  ),
  mode = "pos",
  rank_max = 200,
  gp_order = gp_clustered_order, # This is the magic line
  x_label_every = 1,
  n_label_best = 1,
  free_y = TRUE,
  label_internal = TRUE
)
ggsave(
  p_clustered,
  filename = paste0(figure_path, "TF_skyline_clustered.pdf"),
  width = 20,
  height = 20
)


#####################################################
#####################################################
#####################################################
### Skyscraper plot of TF values across GPs (ordered)
#####################################################
#####################################################
F_pm_TF <- F_pm_norm_col[rownames(F_pm_norm_col) %in% TF_highlights, ]
JJJ <- pheatmap(
  (F_pm_TF),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(
    100
  ),
  breaks = seq(-1, 1, length.out = 101),
  # add grey grid lines
  border_color = "grey10",
  fontsize_row = 6,
  fontsize_col = 8,
  main = "TF × GP loading heatmap (normalized; TFs only)",
  legend = TRUE,
  width = 12,
  height = 10
)
gp_clustered_order <- colnames(rank_clip)[JJ$tree_col$order]
p_clustered2 <- plot_tf_rank_bars_free_y(
  F = F_pm_norm_col,
  TFs = c(
    "Foxp3",
    "Zbtb16",
    "Sox13",
    "Rorc",
    "Pparg",
    "Hes1",
    "Gata3",
    "Tcf7",
    "Tbx21",
    "Tox",
    "Tox2"
  ),
  mode = "pos",
  rank_max = 200,
  gp_order = gp_clustered_order, # This is the magic line
  x_label_every = 1,
  n_label_best = 1,
  free_y = TRUE,
  label_internal = TRUE
)


#####################################################
#####################################################
#####################################################
### TF-GP-gene network plot
#####################################################
#####################################################
# TFs selected based on prior biological knowledge / AUC analysis.
# GPs are connected automatically: any GP where the TF loading > tf_gp_threshold.
# Top genes per GP are the highest-loading genes in F (no separate threshold).
selected_tfs <- c(
  "Tox",
  "Tcf7",
  "Foxp3",
  "Hes1",
  "Rorc",
  "Sox13",
  "Tbx21",
  "Tox2",
  "Zbtb16"
)

network_plot <- plot_tf_gp_network_v2(
  F                = F_pm_norm_col,
  selected_tfs     = selected_tfs,
  tf_gp_threshold  = 0.25,
  top_genes_per_gp = 5,
  gp_spacing       = 1.5,
  node_size_tf     = 7,
  node_size_gp     = 5,
  label_size_tf    = 5.5,
  label_size_gp    = 4.5,
  label_size_gene  = 3.8
)

n_gps_network <- length(unique(unlist(lapply(selected_tfs, function(tf) {
  vals <- setNames(as.numeric(F_pm_norm_col[tf, ]), colnames(F_pm_norm_col))
  names(vals)[is.finite(vals) & vals > 0.25]
}))))

ggsave(
  filename = file.path(figure_path, "TF_GP_gene_network.pdf"),
  plot = network_plot,
  width = 16,
  height = min(48, max(10, n_gps_network * 1.5 * 0.55 + 2)),
  limitsize = FALSE
)
