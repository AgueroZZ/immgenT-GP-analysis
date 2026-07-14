# Render GP-by-group mean-loading heatmaps for the healthy non-thymocyte reference.
#
# The experiment retains only the current filtered raw, normalized, and
# centered heatmaps plus unfiltered all-200-GP centered internal views.

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

center_by_gp_mean <- function(mean_matrix) {
  sweep(mean_matrix, 1L, rowMeans(mean_matrix), "-")
}

dominant_group_order <- function(raw_mean_matrix, fixed_column_order = NULL) {
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
  if (is.null(fixed_column_order)) {
    column_order <- order(
      -dominant_gp_count,
      -colMeans(raw_mean_matrix),
      colnames(raw_mean_matrix)
    )
  } else {
    if (
      length(fixed_column_order) != ncol(raw_mean_matrix) ||
      !identical(sort(fixed_column_order), seq_len(ncol(raw_mean_matrix)))
    ) {
      stop("The fixed column order must be a complete permutation.")
    }
    column_order <- fixed_column_order
  }
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

level2_to_level1_map <- function(meta, groups, level1_order) {
  mapping <- unique(data.frame(
    group = as.character(meta$annotation_level2),
    level1 = as.character(meta$annotation_level1),
    stringsAsFactors = FALSE
  ))
  if (anyDuplicated(mapping$group)) {
    stop("Each annotation_level2 label must map to exactly one annotation_level1 label.")
  }

  group_level1 <- mapping$level1[match(groups, mapping$group)]
  names(group_level1) <- groups
  if (anyNA(group_level1) || any(!group_level1 %in% level1_order)) {
    stop("Every displayed level2 group must map to the Figure 1 level1 order.")
  }
  group_level1
}

level2_column_order <- function(groups, group_level1, level1_order) {
  order(match(group_level1[groups], level1_order), groups)
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

filter_centered_mean_matrix <- function(centered_matrix, centered_mean_cutoff) {
  supported_entries <- centered_matrix >= centered_mean_cutoff
  keep_gp <- rowSums(supported_entries) > 0L
  keep_group <- colSums(supported_entries[keep_gp, , drop = FALSE]) > 0L
  filtered_centered <- centered_matrix[keep_gp, keep_group, drop = FALSE]

  if (nrow(filtered_centered) == 0L || ncol(filtered_centered) < 2L) {
    stop("The centered-mean filter must retain at least one GP and two groups.")
  }

  filtered_centered
}

subset_group_counts <- function(group_counts, groups) {
  group_index <- match(groups, group_counts$group)
  if (anyNA(group_index)) {
    stop("A retained matrix group is absent from the group-count table.")
  }
  group_counts[group_index, , drop = FALSE]
}

write_filter_summary <- function(
    full_matrix,
    filtered_matrix,
    filter_basis,
    filter_value,
    filename
) {
  write.csv(
    data.frame(
      filter_basis = filter_basis,
      filter_value = filter_value,
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
    order_description = NULL,
    centered_limit = NULL,
    group_level1 = NULL,
    level1_palette = NULL
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
  } else if (identical(scale_label, "Within-GP normalized mean loading")) {
    color_fun <- circlize::colorRamp2(
      c(0, 0.5, 1),
      c("#FFFFFF", "#FCAE91", "#99000D")
    )
    legend_at <- c(0, 0.5, 1)
  } else if (identical(scale_label, "Row-centered mean loading")) {
    color_fun <- circlize::colorRamp2(
      c(-centered_limit, 0, centered_limit),
      c("#2166AC", "#FFFFFF", "#B2182B")
    )
    legend_at <- c(-centered_limit, 0, centered_limit)
  } else {
    stop("Unsupported scale label: ", scale_label)
  }

  heatmap_width_mm <- max(180, ncol(matrix) * 4.2)
  heatmap_height_mm <- max(160, nrow(matrix) * 3.5)
  pdf_width_in <- (heatmap_width_mm + 130) / 25.4
  pdf_height_in <- (heatmap_height_mm + 90) / 25.4
  cell_width_mm <- heatmap_width_mm / ncol(matrix)
  cell_height_mm <- heatmap_height_mm / nrow(matrix)
  row_label_fontsize <- min(14, max(9, floor(cell_height_mm * 2.8)))
  column_label_fontsize <- min(14, max(9, floor(cell_width_mm * 2.8)))

  if (is.null(group_level1)) {
    column_annotation <- ComplexHeatmap::HeatmapAnnotation(
      group = factor(colnames(matrix), levels = colnames(matrix)),
      col = list(group = group_palette),
      show_legend = FALSE,
      annotation_name_side = "left",
      annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold"),
      annotation_height = grid::unit(4, "mm")
    )
  } else {
    group_level1 <- group_level1[colnames(matrix)]
    if (anyNA(group_level1) || is.null(level1_palette)) {
      stop("Level2 heatmaps require complete level1 annotations and a palette.")
    }
    column_annotation <- ComplexHeatmap::HeatmapAnnotation(
      level1 = factor(group_level1, levels = names(level1_palette)),
      group = factor(colnames(matrix), levels = colnames(matrix)),
      col = list(level1 = level1_palette, group = group_palette),
      show_legend = FALSE,
      annotation_name_side = "left",
      annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold"),
      annotation_height = grid::unit(c(4, 4), "mm")
    )
  }
  column_title <- if (is.null(order_description)) {
    paste0(scale_label, ": all GPs by healthy non-thymocyte ", group_label)
  } else {
    short_scale_label <- switch(
      scale_label,
      "Raw mean loading" = "Raw GP mean loading",
      "Within-GP normalized mean loading" = "Normalized GP mean loading",
      "Row-centered mean loading" = "Row-centered GP mean loading"
    )
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
    column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
    row_title = "GP",
    row_title_gp = grid::gpar(fontsize = 12),
    row_names_gp = grid::gpar(fontsize = row_label_fontsize),
    column_names_gp = grid::gpar(fontsize = column_label_fontsize),
    column_names_rot = 90,
    heatmap_legend_param = list(
      title = scale_label,
      at = legend_at,
      labels = format(legend_at, trim = TRUE, scientific = FALSE),
      title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 10)
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
organ_centered <- center_by_gp_mean(organ_raw)
level2_centered <- center_by_gp_mean(level2_raw)
level1_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN", "DP")
level2_group_level1 <- level2_to_level1_map(
  meta_reference, colnames(level2_raw), level1_order
)
organ_dominant <- dominant_group_order(organ_raw)
raw_dominant_cutoff <- 0.1
centered_dominant_cutoff <- 0.01
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
organ_centered_dominant_filtered <- filter_centered_mean_matrix(
  organ_centered, centered_dominant_cutoff
)
level2_centered_dominant_filtered <- filter_centered_mean_matrix(
  level2_centered, centered_dominant_cutoff
)
organ_raw_dominant_filtered_order <- dominant_group_order(organ_raw_dominant_filtered)
level2_raw_dominant_filtered_order <- dominant_group_order(
  level2_raw_dominant_filtered,
  level2_column_order(
    colnames(level2_raw_dominant_filtered), level2_group_level1, level1_order
  )
)
level2_centered_full_order <- dominant_group_order(
  level2_raw,
  level2_column_order(colnames(level2_raw), level2_group_level1, level1_order)
)
organ_centered_dominant_filtered_order <- dominant_group_order(
  organ_centered_dominant_filtered
)
level2_centered_dominant_filtered_order <- dominant_group_order(
  level2_centered_dominant_filtered,
  level2_column_order(
    colnames(level2_centered_dominant_filtered), level2_group_level1, level1_order
  )
)

organ_centered_supported <- organ_centered_dominant_filtered >= centered_dominant_cutoff
level2_centered_supported <- level2_centered_dominant_filtered >= centered_dominant_cutoff

stopifnot(
  nrow(organ_raw) == 200L,
  nrow(level2_raw) == 200L,
  all(abs(apply(organ_normalized, 1L, max) - 1) < 1e-10),
  all(abs(apply(level2_normalized, 1L, max) - 1) < 1e-10),
  identical(dimnames(organ_raw_dominant_filtered), dimnames(organ_normalized_dominant_filtered)),
  identical(dimnames(level2_raw_dominant_filtered), dimnames(level2_normalized_dominant_filtered)),
  all(apply(organ_raw_dominant_filtered, 1L, max) >= raw_dominant_cutoff),
  all(apply(level2_raw_dominant_filtered, 1L, max) >= raw_dominant_cutoff),
  all(rowSums(organ_centered_supported) > 0L),
  all(rowSums(level2_centered_supported) > 0L),
  all(colSums(organ_centered_supported) > 0L),
  all(colSums(level2_centered_supported) > 0L),
  "GP37" %in% rownames(organ_centered_dominant_filtered),
  "GP37" %in% rownames(level2_centered_dominant_filtered),
  max(abs(rowMeans(organ_centered))) < 1e-12,
  max(abs(rowMeans(level2_centered))) < 1e-12
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
organ_centered_dominant_filtered_counts <- subset_group_counts(
  organ_result$counts,
  colnames(organ_centered_dominant_filtered)
)
level2_centered_dominant_filtered_counts <- subset_group_counts(
  level2_result$counts,
  colnames(level2_centered_dominant_filtered)
)
organ_centered_dominant_filtered_palette <- organ_palette[colnames(organ_centered_dominant_filtered)]
level2_centered_dominant_filtered_palette <- level2_palette[colnames(level2_centered_dominant_filtered)]
level1_palette <- ZemmourLib::immgent_colors$level1[level1_order]

write_matrix_csv(organ_centered, file.path(output_dir, "organ_simplified_row_centered_mean_loading_matrix.csv"))
write_matrix_csv(level2_centered, file.path(output_dir, "annotation_level2_row_centered_mean_loading_matrix.csv"))
write_matrix_csv(
  organ_raw_dominant_filtered,
  file.path(output_dir, "organ_simplified_raw_mean_loading_raw_mean_ge_0.1_dominant_group_matrix.csv")
)
write_matrix_csv(
  level2_raw_dominant_filtered,
  file.path(output_dir, "annotation_level2_raw_mean_loading_raw_mean_ge_0.1_dominant_group_matrix.csv")
)
write_matrix_csv(
  organ_normalized_dominant_filtered,
  file.path(output_dir, "organ_simplified_row_normalized_mean_loading_raw_mean_ge_0.1_dominant_group_matrix.csv")
)
write_matrix_csv(
  level2_normalized_dominant_filtered,
  file.path(output_dir, "annotation_level2_row_normalized_mean_loading_raw_mean_ge_0.1_dominant_group_matrix.csv")
)
write_matrix_csv(
  organ_centered_dominant_filtered,
  file.path(output_dir, "organ_simplified_row_centered_mean_loading_centered_mean_ge_0.01_dominant_group_matrix.csv")
)
write_matrix_csv(
  level2_centered_dominant_filtered,
  file.path(output_dir, "annotation_level2_row_centered_mean_loading_centered_mean_ge_0.01_dominant_group_matrix.csv")
)
write.csv(organ_result$counts, file.path(output_dir, "organ_simplified_group_cell_counts.csv"), row.names = FALSE, quote = FALSE)
write.csv(level2_result$counts, file.path(output_dir, "annotation_level2_group_cell_counts.csv"), row.names = FALSE, quote = FALSE)
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
write_dominant_group_order_csv(
  organ_centered_dominant_filtered,
  organ_centered_dominant_filtered_order,
  file.path(output_dir, "organ_simplified_centered_mean_ge_0.01_dominant_group_order.csv")
)
write_dominant_group_order_csv(
  level2_centered_dominant_filtered,
  level2_centered_dominant_filtered_order,
  file.path(output_dir, "annotation_level2_centered_mean_ge_0.01_dominant_group_order.csv")
)
write_filter_summary(
  organ_raw,
  organ_raw_dominant_filtered,
  "raw mean cutoff",
  raw_dominant_cutoff,
  file.path(output_dir, "organ_simplified_raw_mean_ge_0.1_dominant_group_filter_summary.csv")
)
write_filter_summary(
  level2_raw,
  level2_raw_dominant_filtered,
  "raw mean cutoff",
  raw_dominant_cutoff,
  file.path(output_dir, "annotation_level2_raw_mean_ge_0.1_dominant_group_filter_summary.csv")
)
write_filter_summary(
  organ_centered,
  organ_centered_dominant_filtered,
  "centered mean cutoff",
  centered_dominant_cutoff,
  file.path(output_dir, "organ_simplified_centered_mean_ge_0.01_dominant_group_filter_summary.csv")
)
write_filter_summary(
  level2_centered,
  level2_centered_dominant_filtered,
  "centered mean cutoff",
  centered_dominant_cutoff,
  file.path(output_dir, "annotation_level2_centered_mean_ge_0.01_dominant_group_filter_summary.csv")
)

raw_limit <- max(c(organ_raw, level2_raw))
centered_limit <- max(abs(c(organ_centered, level2_centered)))

centered_full_description <- "all 200 GPs; dominant-group blocks (within block: dominance gap)"
organ_centered_full_draw <- render_heatmap(
  organ_centered,
  organ_result$counts,
  organ_palette,
  "tissue (organ_simplified)",
  "Row-centered mean loading",
  file.path(output_dir, "organ_simplified_row_centered_mean_loading_full_dominant_group_heatmap.pdf"),
  row_order = organ_dominant$row_order,
  column_order = organ_dominant$column_order,
  order_description = centered_full_description,
  centered_limit = centered_limit
)
write_order_csv(
  organ_centered,
  organ_centered_full_draw,
  file.path(output_dir, "organ_simplified_row_centered_mean_loading_full_dominant_group")
)

level2_centered_full_description <- paste0(
  "all 200 GPs; level2 columns: Figure 1 level1 order (CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP); ",
  "alphabetical within level1; GP rows: dominant-group blocks"
)
level2_centered_full_draw <- render_heatmap(
  level2_centered,
  level2_result$counts,
  level2_palette,
  "cluster (annotation_level2)",
  "Row-centered mean loading",
  file.path(output_dir, "annotation_level2_row_centered_mean_loading_full_level1_order_heatmap.pdf"),
  row_order = level2_centered_full_order$row_order,
  column_order = level2_centered_full_order$column_order,
  order_description = level2_centered_full_description,
  centered_limit = centered_limit,
  group_level1 = level2_group_level1,
  level1_palette = level1_palette
)
write_order_csv(
  level2_centered,
  level2_centered_full_draw,
  file.path(output_dir, "annotation_level2_row_centered_mean_loading_full_level1_order")
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
  paste0(
    "level2 columns: Figure 1 level1 order (CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP); ",
    "alphabetical within level1; GP rows: dominant-group blocks; focus: raw mean >= ",
    raw_dominant_cutoff, " filter"
  ),
  group_level1 = level2_group_level1,
  level1_palette = level1_palette
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
  order_description = paste0(
    "level2 columns: Figure 1 level1 order (CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP); ",
    "alphabetical within level1; GP rows: dominant-group blocks; focus: raw mean >= ",
    raw_dominant_cutoff, " filter"
  ),
  group_level1 = level2_group_level1,
  level1_palette = level1_palette
)

organ_centered_dominant_draw <- render_heatmap(
  organ_centered_dominant_filtered,
  organ_centered_dominant_filtered_counts,
  organ_centered_dominant_filtered_palette,
  "tissue (organ_simplified)",
  "Row-centered mean loading",
  file.path(output_dir, "organ_simplified_row_centered_mean_loading_dominant_group_heatmap.pdf"),
  row_order = organ_centered_dominant_filtered_order$row_order,
  column_order = organ_centered_dominant_filtered_order$column_order,
  order_description = paste0(
    dominant_group_description,
    "; focus: centered mean >= ", centered_dominant_cutoff, " filter"
  ),
  centered_limit = centered_limit
)

level2_centered_dominant_draw <- render_heatmap(
  level2_centered_dominant_filtered,
  level2_centered_dominant_filtered_counts,
  level2_centered_dominant_filtered_palette,
  "cluster (annotation_level2)",
  "Row-centered mean loading",
  file.path(output_dir, "annotation_level2_row_centered_mean_loading_dominant_group_heatmap.pdf"),
  row_order = level2_centered_dominant_filtered_order$row_order,
  column_order = level2_centered_dominant_filtered_order$column_order,
  order_description = paste0(
    "level2 columns: Figure 1 level1 order (CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP); ",
    "alphabetical within level1; GP rows: dominant-group blocks; focus: centered mean >= ",
    centered_dominant_cutoff, " filter"
  ),
  centered_limit = centered_limit,
  group_level1 = level2_group_level1,
  level1_palette = level1_palette
)

summary_lines <- c(
  "Healthy non-thymocyte GP mean-loading heatmaps",
  paste("Reference cells:", nrow(L_reference)),
  paste("GPs:", ncol(L_reference)),
  paste("Observed organ_simplified groups:", ncol(organ_raw)),
  paste("Observed annotation_level2 groups:", ncol(level2_raw)),
  paste("Shared raw-mean color maximum:", format(raw_limit, digits = 8)),
  paste("Shared centered color maximum:", format(centered_limit, digits = 8)),
  "Normalized matrices: each GP row divided by its maximum group mean.",
  "Centered matrices: each GP row has its mean across groups subtracted.",
  "Tissue dominant-group ordering: columns by number of dominant GPs; rows grouped by dominant group and ordered by decreasing dominance gap.",
  "Final level2 ordering: Figure 1 level1 order (CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP), alphabetized within level1; rows grouped by dominant level2.",
  paste("Final raw dominant heatmap cutoff:", raw_dominant_cutoff),
  paste("Final centered dominant heatmap cutoff:", centered_dominant_cutoff),
  paste("Filtered tissue raw dominant heatmap:", nrow(organ_raw_dominant_filtered), "GPs x", ncol(organ_raw_dominant_filtered), "groups"),
  paste("Filtered level2 raw dominant heatmap:", nrow(level2_raw_dominant_filtered), "GPs x", ncol(level2_raw_dominant_filtered), "groups"),
  paste("Filtered tissue centered dominant heatmap:", nrow(organ_centered_dominant_filtered), "GPs x", ncol(organ_centered_dominant_filtered), "groups"),
  paste("Filtered level2 centered dominant heatmap:", nrow(level2_centered_dominant_filtered), "GPs x", ncol(level2_centered_dominant_filtered), "groups"),
  "Final raw and normalized dominant-group heatmaps share their raw-filtered rows, columns, and order; centered heatmaps use independent positive centered mean >= 0.01 filtering and ordering.",
  paste("Full centered tissue heatmap:", nrow(organ_centered), "GPs x", ncol(organ_centered), "groups"),
  paste("Full centered level2 heatmap:", nrow(level2_centered), "GPs x", ncol(level2_centered), "groups"),
  paste("Filtered tissue PDF dimensions (inches):", round(organ_raw_dominant_draw$width_in, 1), "x", round(organ_raw_dominant_draw$height_in, 1)),
  paste("Filtered level2 PDF dimensions (inches):", round(level2_raw_dominant_draw$width_in, 1), "x", round(level2_raw_dominant_draw$height_in, 1))
)
writeLines(summary_lines, file.path(output_dir, "heatmap_summary.txt"))

message("Wrote healthy non-thymocyte GP mean-loading heatmaps to ", normalizePath(output_dir))
