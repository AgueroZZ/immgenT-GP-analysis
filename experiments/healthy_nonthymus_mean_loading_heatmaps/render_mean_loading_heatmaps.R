# Render GP-by-group mean-loading heatmaps for the healthy non-thymocyte reference.
#
# The hierarchical and triangular-first PDFs are intentionally kept in this
# experiment directory. They do not modify any formal figure, workflowr page,
# or published output.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ZemmourLib)
})

if (!file.exists("code/R/setup_data.R")) {
  stop("Run this script from the immgenT-GP-analysis repository root.")
}

source("code/R/setup_data.R")

output_dir <- "experiments/healthy_nonthymus_mean_loading_heatmaps"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mean_loading_by_group <- function(L_mat, labels) {
  if (length(labels) != nrow(L_mat)) {
    stop("`labels` must contain one value for each row of `L_mat`.")
  }
  if (anyNA(labels) || any(labels == "")) {
    stop("Group labels must be present for every retained cell.")
  }

  labels <- droplevels(factor(as.character(labels)))
  group_sums <- rowsum(L_mat, group = labels, reorder = TRUE)
  group_counts <- as.integer(table(labels)[rownames(group_sums)])
  mean_matrix <- t(sweep(group_sums, 1L, group_counts, "/"))

  list(
    matrix = mean_matrix,
    counts = data.frame(group = rownames(group_sums), n_cells = group_counts)
  )
}

normalize_by_gp_max <- function(mean_matrix) {
  row_max <- apply(mean_matrix, 1L, max)
  if (any(!is.finite(row_max)) || any(row_max <= 0)) {
    stop("Every GP must have a finite positive maximum mean loading.")
  }
  sweep(mean_matrix, 1L, row_max, "/")
}

triangular_first_order <- function(normalized_matrix, support_cutoff = 0.50) {
  visible_mask <- normalized_matrix >= support_cutoff
  row_visible_count <- rowSums(visible_mask)
  column_visible_count <- colSums(visible_mask)

  if (any(row_visible_count == 0L) || any(column_visible_count == 0L)) {
    stop("The triangular support cutoff leaves an empty GP row or group column.")
  }

  gp_number <- suppressWarnings(as.integer(sub("^GP", "", rownames(normalized_matrix))))
  if (anyNA(gp_number)) {
    stop("GP row names must have the form `GP<number>` for deterministic ordering.")
  }

  column_order <- order(-column_visible_count, colnames(normalized_matrix))
  visible_mask_ordered_columns <- visible_mask[, column_order, drop = FALSE]
  rightmost_visible_column <- apply(
    visible_mask_ordered_columns,
    1L,
    function(values) max(which(values))
  )
  column_rarity_weights <- seq_len(ncol(visible_mask_ordered_columns))^2
  row_rarity_score <- as.numeric(
    visible_mask_ordered_columns %*% column_rarity_weights
  )
  row_order <- order(
    -rightmost_visible_column,
    -row_visible_count,
    -row_rarity_score,
    gp_number
  )

  ordered_right_boundary <- rightmost_visible_column[row_order]
  ordered_column_support <- column_visible_count[column_order]
  if (any(diff(ordered_right_boundary) > 0L) || any(diff(ordered_column_support) > 0L)) {
    stop("Triangular ordering invariants failed.")
  }

  list(
    support_cutoff = support_cutoff,
    row_order = row_order,
    column_order = column_order,
    row_visible_count = row_visible_count,
    column_visible_count = column_visible_count,
    rightmost_visible_column = rightmost_visible_column,
    row_rarity_score = row_rarity_score,
    row_right_boundary_increases = sum(diff(ordered_right_boundary) > 0L),
    column_support_increases = sum(diff(ordered_column_support) > 0L)
  )
}

dominant_group_order <- function(raw_mean_matrix) {
  if (ncol(raw_mean_matrix) < 2L) {
    stop("Dominant-group ordering requires at least two groups.")
  }

  gp_number <- suppressWarnings(as.integer(sub("^GP", "", rownames(raw_mean_matrix))))
  if (anyNA(gp_number)) {
    stop("GP row names must have the form `GP<number>` for deterministic ordering.")
  }

  dominant_index <- max.col(raw_mean_matrix, ties.method = "first")
  dominant_mean <- raw_mean_matrix[cbind(seq_len(nrow(raw_mean_matrix)), dominant_index)]
  second_mean <- apply(
    raw_mean_matrix,
    1L,
    function(values) sort(values, decreasing = TRUE)[2L]
  )
  dominance_gap <- dominant_mean - second_mean
  dominant_gp_count <- tabulate(dominant_index, nbins = ncol(raw_mean_matrix))

  column_order <- order(
    -dominant_gp_count,
    -colMeans(raw_mean_matrix),
    colnames(raw_mean_matrix)
  )
  dominant_group_position <- match(dominant_index, column_order)
  row_order <- order(
    dominant_group_position,
    -dominance_gap,
    -dominant_mean,
    gp_number
  )

  if (any(diff(dominant_group_position[row_order]) < 0L)) {
    stop("Dominant-group blocks are not monotone after ordering.")
  }

  list(
    row_order = row_order,
    column_order = column_order,
    dominant_index = dominant_index,
    dominant_mean = dominant_mean,
    second_mean = second_mean,
    dominance_gap = dominance_gap,
    dominant_gp_count = dominant_gp_count,
    dominant_group_position = dominant_group_position
  )
}

palette_for_groups <- function(groups, palette, label) {
  missing <- setdiff(groups, names(palette))
  if (length(missing) > 0L) {
    stop(
      "The canonical ", label, " palette lacks: ",
      paste(missing, collapse = ", ")
    )
  }
  palette[groups]
}

write_matrix_csv <- function(matrix, filename) {
  write.csv(
    cbind(GP = rownames(matrix), matrix),
    filename,
    row.names = FALSE,
    quote = FALSE
  )
}

filter_raw_mean_matrix <- function(raw_matrix, raw_mean_cutoff) {
  keep_gp <- rowSums(raw_matrix >= raw_mean_cutoff) > 0L
  filtered_raw <- raw_matrix[keep_gp, , drop = FALSE]
  keep_group <- colSums(filtered_raw >= raw_mean_cutoff) > 0L
  filtered_raw <- filtered_raw[, keep_group, drop = FALSE]

  if (nrow(filtered_raw) == 0L || ncol(filtered_raw) < 2L) {
    stop("The raw-mean filter must retain at least one GP and two groups.")
  }

  filtered_raw
}

subset_group_counts <- function(group_counts, groups) {
  group_index <- match(groups, group_counts$group)
  if (anyNA(group_index)) {
    stop("A retained matrix group is absent from the group-count table.")
  }
  group_counts[group_index, , drop = FALSE]
}

write_raw_filter_summary <- function(full_matrix, filtered_matrix, raw_mean_cutoff, filename) {
  write.csv(
    data.frame(
      raw_mean_cutoff = raw_mean_cutoff,
      full_gp_count = nrow(full_matrix),
      retained_gp_count = nrow(filtered_matrix),
      full_group_count = ncol(full_matrix),
      retained_group_count = ncol(filtered_matrix)
    ),
    filename,
    row.names = FALSE,
    quote = FALSE
  )
}

render_heatmap <- function(
    matrix,
    group_counts,
    group_palette,
    group_label,
    scale_label,
    filename,
    raw_limit = NULL,
    row_order = NULL,
    column_order = NULL,
    order_description = NULL
) {
  if (!identical(colnames(matrix), group_counts$group)) {
    stop("Matrix columns and group-count labels are not aligned.")
  }
  using_fixed_order <- !is.null(row_order) || !is.null(column_order)
  if (xor(is.null(row_order), is.null(column_order))) {
    stop("Supply both `row_order` and `column_order`, or neither.")
  }
  if (using_fixed_order) {
    if (
      length(row_order) != nrow(matrix) || length(column_order) != ncol(matrix) ||
      !identical(sort(row_order), seq_len(nrow(matrix))) ||
      !identical(sort(column_order), seq_len(ncol(matrix)))
    ) {
      stop("Fixed row and column orders must be complete permutations.")
    }
  }

  if (identical(scale_label, "Raw mean loading")) {
    color_fun <- circlize::colorRamp2(
      c(0, raw_limit * 0.15, raw_limit),
      c("#FFFFFF", "#FCAE91", "#99000D")
    )
    legend_at <- c(0, raw_limit / 2, raw_limit)
  } else {
    color_fun <- circlize::colorRamp2(
      c(0, 0.5, 1),
      c("#FFFFFF", "#FCAE91", "#99000D")
    )
    legend_at <- c(0, 0.5, 1)
  }

  heatmap_width_mm <- max(180, ncol(matrix) * 4.2)
  heatmap_height_mm <- max(480, nrow(matrix) * 3.0)
  pdf_width_in <- (heatmap_width_mm + 130) / 25.4
  pdf_height_in <- (heatmap_height_mm + 120) / 25.4

  column_annotation <- ComplexHeatmap::HeatmapAnnotation(
    group = factor(colnames(matrix), levels = colnames(matrix)),
    col = list(group = group_palette),
    show_legend = FALSE,
    annotation_name_side = "left",
    annotation_height = grid::unit(4, "mm")
  )
  column_title <- if (is.null(order_description)) {
    paste0(scale_label, ": all GPs by healthy non-thymocyte ", group_label)
  } else {
    short_scale_label <- if (identical(scale_label, "Raw mean loading")) {
      "Raw GP mean loading"
    } else {
      "Normalized GP mean loading"
    }
    paste0(
      short_scale_label, ": ", group_label,
      "\n", order_description
    )
  }

  heatmap <- ComplexHeatmap::Heatmap(
    matrix,
    name = scale_label,
    col = color_fun,
    cluster_rows = !using_fixed_order,
    cluster_columns = !using_fixed_order,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    row_order = row_order,
    column_order = column_order,
    top_annotation = column_annotation,
    column_title = column_title,
    column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
    row_title = "GP",
    row_title_gp = grid::gpar(fontsize = 10),
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 5),
    column_names_rot = 45,
    heatmap_legend_param = list(
      title = scale_label,
      at = legend_at,
      labels = format(legend_at, trim = TRUE, scientific = FALSE)
    ),
    width = grid::unit(heatmap_width_mm, "mm"),
    height = grid::unit(heatmap_height_mm, "mm"),
    use_raster = TRUE,
    raster_quality = 2
  )

  grDevices::pdf(filename, width = pdf_width_in, height = pdf_height_in)
  drawn_heatmap <- ComplexHeatmap::draw(
    heatmap,
    heatmap_legend_side = "right",
    padding = grid::unit(c(8, 8, 8, 8), "mm")
  )
  grDevices::dev.off()

  list(
    row_order = ComplexHeatmap::row_order(drawn_heatmap),
    column_order = ComplexHeatmap::column_order(drawn_heatmap),
    width_in = pdf_width_in,
    height_in = pdf_height_in
  )
}

write_triangular_order_csv <- function(matrix, triangular_order, filename) {
  row_table <- data.frame(
    dimension = "GP row",
    order = seq_len(nrow(matrix)),
    label = rownames(matrix)[triangular_order$row_order],
    visible_count = unname(triangular_order$row_visible_count[triangular_order$row_order]),
    total_entries = ncol(matrix),
    rightmost_visible_column = unname(triangular_order$rightmost_visible_column[triangular_order$row_order]),
    rarity_score = unname(triangular_order$row_rarity_score[triangular_order$row_order]),
    support_cutoff = triangular_order$support_cutoff
  )
  column_table <- data.frame(
    dimension = "Group column",
    order = seq_len(ncol(matrix)),
    label = colnames(matrix)[triangular_order$column_order],
    visible_count = unname(triangular_order$column_visible_count[triangular_order$column_order]),
    total_entries = nrow(matrix),
    rightmost_visible_column = NA_integer_,
    rarity_score = NA_real_,
    support_cutoff = triangular_order$support_cutoff
  )
  write.csv(rbind(row_table, column_table), filename, row.names = FALSE, quote = FALSE)
}

write_dominant_group_order_csv <- function(matrix, dominant_order, filename) {
  row_table <- data.frame(
    dimension = "GP row",
    order = seq_len(nrow(matrix)),
    label = rownames(matrix)[dominant_order$row_order],
    dominant_group = colnames(matrix)[dominant_order$dominant_index[dominant_order$row_order]],
    dominant_group_position = unname(dominant_order$dominant_group_position[dominant_order$row_order]),
    dominant_mean = unname(dominant_order$dominant_mean[dominant_order$row_order]),
    second_mean = unname(dominant_order$second_mean[dominant_order$row_order]),
    dominance_gap = unname(dominant_order$dominance_gap[dominant_order$row_order]),
    dominant_gp_count = NA_integer_
  )
  column_table <- data.frame(
    dimension = "Group column",
    order = seq_len(ncol(matrix)),
    label = colnames(matrix)[dominant_order$column_order],
    dominant_group = NA_character_,
    dominant_group_position = seq_len(ncol(matrix)),
    dominant_mean = NA_real_,
    second_mean = NA_real_,
    dominance_gap = NA_real_,
    dominant_gp_count = unname(dominant_order$dominant_gp_count[dominant_order$column_order])
  )
  write.csv(rbind(row_table, column_table), filename, row.names = FALSE, quote = FALSE)
}

write_order_csv <- function(matrix, draw_result, prefix) {
  write.csv(
    data.frame(order = seq_len(nrow(matrix)), GP = rownames(matrix)[draw_result$row_order]),
    paste0(prefix, "_row_order.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  write.csv(
    data.frame(
      order = seq_len(ncol(matrix)),
      group = colnames(matrix)[draw_result$column_order]
    ),
    paste0(prefix, "_column_order.csv"),
    row.names = FALSE,
    quote = FALSE
  )
}

gp_data <- load_gp_data()
meta <- gp_data$seurat_meta_filtered

healthy_nonthymus <-
  meta$condition_broad == "healthy" &
  meta$annotation_level1 != "thymocyte"

if (anyNA(healthy_nonthymus)) {
  stop("Healthy non-thymocyte selection contains missing values.")
}

L_reference <- gp_data$L_pm_filtered[healthy_nonthymus, , drop = FALSE]
meta_reference <- meta[healthy_nonthymus, , drop = FALSE]

if (ncol(L_reference) != 200L || nrow(L_reference) != nrow(meta_reference)) {
  stop("The healthy non-thymocyte reference matrix has unexpected dimensions.")
}
if (anyNA(L_reference)) {
  stop("The loading matrix contains missing values.")
}

organ_result <- mean_loading_by_group(L_reference, meta_reference$organ_simplified)
level2_result <- mean_loading_by_group(L_reference, meta_reference$annotation_level2)

organ_raw <- organ_result$matrix
level2_raw <- level2_result$matrix
organ_normalized <- normalize_by_gp_max(organ_raw)
level2_normalized <- normalize_by_gp_max(level2_raw)
support_cutoff <- 0.50
organ_triangular <- triangular_first_order(organ_normalized, support_cutoff)
level2_triangular <- triangular_first_order(level2_normalized, support_cutoff)
organ_dominant <- dominant_group_order(organ_raw)
level2_dominant <- dominant_group_order(level2_raw)
raw_dominant_cutoff <- 0.1
organ_raw_dominant_filtered <- filter_raw_mean_matrix(organ_raw, raw_dominant_cutoff)
level2_raw_dominant_filtered <- filter_raw_mean_matrix(level2_raw, raw_dominant_cutoff)
organ_normalized_dominant_filtered <- organ_normalized[
  rownames(organ_raw_dominant_filtered),
  colnames(organ_raw_dominant_filtered),
  drop = FALSE
]
level2_normalized_dominant_filtered <- level2_normalized[
  rownames(level2_raw_dominant_filtered),
  colnames(level2_raw_dominant_filtered),
  drop = FALSE
]
organ_raw_dominant_filtered_order <- dominant_group_order(organ_raw_dominant_filtered)
level2_raw_dominant_filtered_order <- dominant_group_order(level2_raw_dominant_filtered)

stopifnot(
  nrow(organ_raw) == 200L,
  nrow(level2_raw) == 200L,
  all(abs(apply(organ_normalized, 1L, max) - 1) < 1e-10),
  all(abs(apply(level2_normalized, 1L, max) - 1) < 1e-10),
  identical(dimnames(organ_raw_dominant_filtered), dimnames(organ_normalized_dominant_filtered)),
  identical(dimnames(level2_raw_dominant_filtered), dimnames(level2_normalized_dominant_filtered)),
  all(apply(organ_raw_dominant_filtered, 1L, max) >= raw_dominant_cutoff),
  all(apply(level2_raw_dominant_filtered, 1L, max) >= raw_dominant_cutoff)
)

organ_palette <- palette_for_groups(
  colnames(organ_raw),
  ZemmourLib::immgent_colors$organ_simplified,
  "organ_simplified"
)
level2_palette <- palette_for_groups(
  colnames(level2_raw),
  ZemmourLib::immgent_colors$level2,
  "level2"
)
organ_raw_dominant_filtered_counts <- subset_group_counts(
  organ_result$counts,
  colnames(organ_raw_dominant_filtered)
)
level2_raw_dominant_filtered_counts <- subset_group_counts(
  level2_result$counts,
  colnames(level2_raw_dominant_filtered)
)
organ_raw_dominant_filtered_palette <- organ_palette[colnames(organ_raw_dominant_filtered)]
level2_raw_dominant_filtered_palette <- level2_palette[colnames(level2_raw_dominant_filtered)]

write_matrix_csv(organ_raw, file.path(output_dir, "organ_simplified_raw_mean_loading_matrix.csv"))
write_matrix_csv(organ_normalized, file.path(output_dir, "organ_simplified_row_normalized_mean_loading_matrix.csv"))
write_matrix_csv(level2_raw, file.path(output_dir, "annotation_level2_raw_mean_loading_matrix.csv"))
write_matrix_csv(level2_normalized, file.path(output_dir, "annotation_level2_row_normalized_mean_loading_matrix.csv"))
write_matrix_csv(
  organ_normalized_dominant_filtered,
  file.path(output_dir, "organ_simplified_row_normalized_mean_loading_raw_mean_ge_0.1_dominant_group_matrix.csv")
)
write_matrix_csv(
  level2_normalized_dominant_filtered,
  file.path(output_dir, "annotation_level2_row_normalized_mean_loading_raw_mean_ge_0.1_dominant_group_matrix.csv")
)
write.csv(organ_result$counts, file.path(output_dir, "organ_simplified_group_cell_counts.csv"), row.names = FALSE, quote = FALSE)
write.csv(level2_result$counts, file.path(output_dir, "annotation_level2_group_cell_counts.csv"), row.names = FALSE, quote = FALSE)
write_triangular_order_csv(
  organ_normalized,
  organ_triangular,
  file.path(output_dir, "organ_simplified_triangular_order.csv")
)
write_triangular_order_csv(
  level2_normalized,
  level2_triangular,
  file.path(output_dir, "annotation_level2_triangular_order.csv")
)
write_dominant_group_order_csv(
  organ_raw,
  organ_dominant,
  file.path(output_dir, "organ_simplified_dominant_group_order.csv")
)
write_dominant_group_order_csv(
  level2_raw,
  level2_dominant,
  file.path(output_dir, "annotation_level2_dominant_group_order.csv")
)
write_dominant_group_order_csv(
  organ_raw_dominant_filtered,
  organ_raw_dominant_filtered_order,
  file.path(output_dir, "organ_simplified_raw_mean_ge_0.1_dominant_group_order.csv")
)
write_dominant_group_order_csv(
  level2_raw_dominant_filtered,
  level2_raw_dominant_filtered_order,
  file.path(output_dir, "annotation_level2_raw_mean_ge_0.1_dominant_group_order.csv")
)
write_raw_filter_summary(
  organ_raw,
  organ_raw_dominant_filtered,
  raw_dominant_cutoff,
  file.path(output_dir, "organ_simplified_raw_mean_ge_0.1_dominant_group_filter_summary.csv")
)
write_raw_filter_summary(
  level2_raw,
  level2_raw_dominant_filtered,
  raw_dominant_cutoff,
  file.path(output_dir, "annotation_level2_raw_mean_ge_0.1_dominant_group_filter_summary.csv")
)

raw_limit <- max(c(organ_raw, level2_raw))

organ_raw_draw <- render_heatmap(
  organ_raw,
  organ_result$counts,
  organ_palette,
  "tissue (organ_simplified)",
  "Raw mean loading",
  file.path(output_dir, "organ_simplified_raw_mean_loading_heatmap.pdf"),
  raw_limit
)
write_order_csv(organ_raw, organ_raw_draw, file.path(output_dir, "organ_simplified_raw_mean_loading"))

organ_normalized_draw <- render_heatmap(
  organ_normalized,
  organ_result$counts,
  organ_palette,
  "tissue (organ_simplified)",
  "Within-GP normalized mean loading",
  file.path(output_dir, "organ_simplified_row_normalized_mean_loading_heatmap.pdf")
)
write_order_csv(organ_normalized, organ_normalized_draw, file.path(output_dir, "organ_simplified_row_normalized_mean_loading"))

level2_raw_draw <- render_heatmap(
  level2_raw,
  level2_result$counts,
  level2_palette,
  "cluster (annotation_level2)",
  "Raw mean loading",
  file.path(output_dir, "annotation_level2_raw_mean_loading_heatmap.pdf"),
  raw_limit
)
write_order_csv(level2_raw, level2_raw_draw, file.path(output_dir, "annotation_level2_raw_mean_loading"))

level2_normalized_draw <- render_heatmap(
  level2_normalized,
  level2_result$counts,
  level2_palette,
  "cluster (annotation_level2)",
  "Within-GP normalized mean loading",
  file.path(output_dir, "annotation_level2_row_normalized_mean_loading_heatmap.pdf")
)
write_order_csv(level2_normalized, level2_normalized_draw, file.path(output_dir, "annotation_level2_row_normalized_mean_loading"))

triangular_description <- paste0(
  "triangular-first ordering (relative support >= ", support_cutoff, ")"
)

organ_raw_triangular_draw <- render_heatmap(
  organ_raw,
  organ_result$counts,
  organ_palette,
  "tissue (organ_simplified)",
  "Raw mean loading",
  file.path(output_dir, "organ_simplified_raw_mean_loading_triangular_heatmap.pdf"),
  raw_limit,
  organ_triangular$row_order,
  organ_triangular$column_order,
  triangular_description
)

organ_normalized_triangular_draw <- render_heatmap(
  organ_normalized,
  organ_result$counts,
  organ_palette,
  "tissue (organ_simplified)",
  "Within-GP normalized mean loading",
  file.path(output_dir, "organ_simplified_row_normalized_mean_loading_triangular_heatmap.pdf"),
  row_order = organ_triangular$row_order,
  column_order = organ_triangular$column_order,
  order_description = triangular_description
)

level2_raw_triangular_draw <- render_heatmap(
  level2_raw,
  level2_result$counts,
  level2_palette,
  "cluster (annotation_level2)",
  "Raw mean loading",
  file.path(output_dir, "annotation_level2_raw_mean_loading_triangular_heatmap.pdf"),
  raw_limit,
  level2_triangular$row_order,
  level2_triangular$column_order,
  triangular_description
)

level2_normalized_triangular_draw <- render_heatmap(
  level2_normalized,
  level2_result$counts,
  level2_palette,
  "cluster (annotation_level2)",
  "Within-GP normalized mean loading",
  file.path(output_dir, "annotation_level2_row_normalized_mean_loading_triangular_heatmap.pdf"),
  row_order = level2_triangular$row_order,
  column_order = level2_triangular$column_order,
  order_description = triangular_description
)

dominant_group_description <- "dominant-group blocks (within block: dominance gap)"
filtered_dominant_group_description <- paste0(
  dominant_group_description,
  "; focus: raw mean >= ", raw_dominant_cutoff, " filter"
)

organ_raw_dominant_draw <- render_heatmap(
  organ_raw_dominant_filtered,
  organ_raw_dominant_filtered_counts,
  organ_raw_dominant_filtered_palette,
  "tissue (organ_simplified)",
  "Raw mean loading",
  file.path(output_dir, "organ_simplified_raw_mean_loading_dominant_group_heatmap.pdf"),
  raw_limit,
  organ_raw_dominant_filtered_order$row_order,
  organ_raw_dominant_filtered_order$column_order,
  filtered_dominant_group_description
)

organ_normalized_dominant_draw <- render_heatmap(
  organ_normalized_dominant_filtered,
  organ_raw_dominant_filtered_counts,
  organ_raw_dominant_filtered_palette,
  "tissue (organ_simplified)",
  "Within-GP normalized mean loading",
  file.path(output_dir, "organ_simplified_row_normalized_mean_loading_dominant_group_heatmap.pdf"),
  row_order = organ_raw_dominant_filtered_order$row_order,
  column_order = organ_raw_dominant_filtered_order$column_order,
  order_description = filtered_dominant_group_description
)

level2_raw_dominant_draw <- render_heatmap(
  level2_raw_dominant_filtered,
  level2_raw_dominant_filtered_counts,
  level2_raw_dominant_filtered_palette,
  "cluster (annotation_level2)",
  "Raw mean loading",
  file.path(output_dir, "annotation_level2_raw_mean_loading_dominant_group_heatmap.pdf"),
  raw_limit,
  level2_raw_dominant_filtered_order$row_order,
  level2_raw_dominant_filtered_order$column_order,
  filtered_dominant_group_description
)

level2_normalized_dominant_draw <- render_heatmap(
  level2_normalized_dominant_filtered,
  level2_raw_dominant_filtered_counts,
  level2_raw_dominant_filtered_palette,
  "cluster (annotation_level2)",
  "Within-GP normalized mean loading",
  file.path(output_dir, "annotation_level2_row_normalized_mean_loading_dominant_group_heatmap.pdf"),
  row_order = level2_raw_dominant_filtered_order$row_order,
  column_order = level2_raw_dominant_filtered_order$column_order,
  order_description = filtered_dominant_group_description
)

summary_lines <- c(
  "Healthy non-thymocyte GP mean-loading heatmaps",
  paste("Reference cells:", nrow(L_reference)),
  paste("GPs:", ncol(L_reference)),
  paste("Observed organ_simplified groups:", ncol(organ_raw)),
  paste("Observed annotation_level2 groups:", ncol(level2_raw)),
  paste("Shared raw-mean color maximum:", format(raw_limit, digits = 8)),
  "Clustering: Euclidean distance with complete linkage on both rows and columns.",
  "Normalized matrices: each GP row divided by its maximum group mean.",
  paste("Triangular support cutoff:", support_cutoff),
  "Triangular ordering: columns by decreasing support; rows by rightmost support boundary, support count, and rarity score.",
  paste("Organ triangular right-boundary increases:", organ_triangular$row_right_boundary_increases),
  paste("Level2 triangular right-boundary increases:", level2_triangular$row_right_boundary_increases),
  "Dominant-group ordering: columns by number of dominant GPs; rows grouped by dominant group and ordered by decreasing dominance gap.",
  paste("Organ dominant groups represented:", sum(organ_dominant$dominant_gp_count > 0L)),
  paste("Level2 dominant groups represented:", sum(level2_dominant$dominant_gp_count > 0L)),
  paste("Final raw dominant heatmap cutoff:", raw_dominant_cutoff),
  paste("Filtered tissue raw dominant heatmap:", nrow(organ_raw_dominant_filtered), "GPs x", ncol(organ_raw_dominant_filtered), "groups"),
  paste("Filtered level2 raw dominant heatmap:", nrow(level2_raw_dominant_filtered), "GPs x", ncol(level2_raw_dominant_filtered), "groups"),
  "Final raw and normalized dominant-group heatmaps share the raw-filtered rows, columns, and order for each grouping.",
  paste("Organ PDF dimensions (inches):", round(organ_raw_draw$width_in, 1), "x", round(organ_raw_draw$height_in, 1)),
  paste("Level2 PDF dimensions (inches):", round(level2_raw_draw$width_in, 1), "x", round(level2_raw_draw$height_in, 1))
)
writeLines(summary_lines, file.path(output_dir, "heatmap_summary.txt"))

message("Wrote healthy non-thymocyte GP mean-loading heatmaps to ", normalizePath(output_dir))
