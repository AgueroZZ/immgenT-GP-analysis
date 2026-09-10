# Row-centered group mean-loading heatmaps.
#
# Shared by script/FigureS2d.R (GP activity across the level-2 clusters) and
# script/Figure4.R (panel 4b). Both show the same quantity: for each GP, the
# mean loading within a group of cells minus that GP's mean across all the
# groups shown, so a row says where a program is more or less active than its
# own average, not how large its loading is.
#
# The column annotation comes in three shapes: one colour bar naming the
# columns; that bar under a level1 bar; or -- pass `group_annotation` --
# Figure 1d's arrangement, a level1 bar and a categorical bar, both with a
# legend and no per-column bar at all.
#
# --- internal ---
# Extracted verbatim (bar the added `palette`/`legend_title` arguments and the
# `keep_gps` selector) from a retired activation-era FigureS4.R, which drew the
# all-200-GP tissue and cluster heatmaps as one two-panel Extended Data figure
# before the 2026-09-09 reorganisation split them apart. The `group_annotation`
# branch was added on 2026-09-10 for Extended Data Figure 2d.
# --- end internal ---

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# The two diverging ramps the figures use, each written low -> mid -> high, so
# the first colour is the negative end. (a)-style panels keep the published
# blue-white-red; green-white-purple is kept for the tissue/lineage pair of the
# internal experiments/fig_n5/ draft, so the two halves of one figure cannot be
# mistaken for each other. Purple is the positive end there by Ziang's call on
# 2026-09-09.
heatmap_palettes <- list(
  blue_red = c("#2166AC", "#FFFFFF", "#B2182B"),
  green_purple = c("#1B7837", "#FFFFFF", "#762A83")
)

# Healthy non-thymocyte reference: the cells every panel here is computed on.
healthy_nonthymocyte_reference <- function(gp_data) {
  meta <- gp_data$seurat_meta_filtered
  keep <- meta$condition_broad == "healthy" & meta$annotation_level1 != "thymocyte"

  if (anyNA(keep)) {
    stop("Healthy non-thymocyte selection contains missing values.")
  }

  L <- gp_data$L_pm_filtered[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  if (ncol(L) != 200L || nrow(L) != nrow(meta) || anyNA(L)) {
    stop("The healthy non-thymocyte loading matrix has unexpected dimensions or missing values.")
  }

  list(L = L, meta = meta)
}

mean_loading_by_group <- function(L_mat, labels) {
  if (length(labels) != nrow(L_mat) || anyNA(labels) || any(labels == "")) {
    stop("Group labels must be present for every retained cell.")
  }

  labels <- droplevels(factor(as.character(labels)))
  group_sums <- rowsum(L_mat, group = labels, reorder = TRUE)
  group_counts <- as.integer(table(labels)[rownames(group_sums)])

  list(
    matrix = t(sweep(group_sums, 1L, group_counts, "/")),
    counts = data.frame(group = rownames(group_sums), n_cells = group_counts)
  )
}

center_by_gp_mean <- function(mean_matrix) {
  sweep(mean_matrix, 1L, rowMeans(mean_matrix), "-")
}

dominant_group_order <- function(raw_mean_matrix, fixed_column_order = NULL) {
  gp_number <- suppressWarnings(as.integer(sub("^GP", "", rownames(raw_mean_matrix))))
  if (ncol(raw_mean_matrix) < 2L || anyNA(gp_number)) {
    stop("Dominant-group ordering requires at least two groups and GP<number> row names.")
  }

  dominant_index <- max.col(raw_mean_matrix, ties.method = "first")
  dominant_mean <- raw_mean_matrix[cbind(seq_len(nrow(raw_mean_matrix)), dominant_index)]
  second_mean <- apply(raw_mean_matrix, 1L, function(values) sort(values, decreasing = TRUE)[2L])
  dominance_gap <- dominant_mean - second_mean

  if (is.null(fixed_column_order)) {
    dominant_gp_count <- tabulate(dominant_index, nbins = ncol(raw_mean_matrix))
    column_order <- order(-dominant_gp_count, -colMeans(raw_mean_matrix), colnames(raw_mean_matrix))
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
  row_order <- order(dominant_group_position, -dominance_gap, -dominant_mean, gp_number)

  if (any(diff(dominant_group_position[row_order]) < 0L)) {
    stop("Dominant-group blocks are not monotone after ordering.")
  }

  list(row_order = row_order, column_order = column_order)
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
    stop("The canonical ", label, " palette lacks: ", paste(missing, collapse = ", "))
  }
  palette[groups]
}

render_centered_heatmap <- function(
    matrix,
    group_palette,
    group_label,
    filename,
    row_order,
    column_order,
    centered_color_limit,
    order_description,
    group_level1 = NULL,
    level1_palette = NULL,
    palette = heatmap_palettes$blue_red,
    group_annotation = NULL,
    group_annotation_palette = NULL
) {
  if (
    length(row_order) != nrow(matrix) || length(column_order) != ncol(matrix) ||
    !identical(sort(row_order), seq_len(nrow(matrix))) ||
    !identical(sort(column_order), seq_len(ncol(matrix)))
  ) {
    stop("Fixed row and column orders must be complete permutations.")
  }

  color_fun <- circlize::colorRamp2(
    c(-centered_color_limit, 0, centered_color_limit),
    palette
  )
  legend_at <- c(-centered_color_limit, 0, centered_color_limit)

  heatmap_width_mm <- max(180, ncol(matrix) * 4.2)
  heatmap_height_mm <- max(160, nrow(matrix) * 3.5)
  pdf_width_in <- (heatmap_width_mm + 130) / 25.4
  pdf_height_in <- (heatmap_height_mm + 90) / 25.4
  cell_width_mm <- heatmap_width_mm / ncol(matrix)
  cell_height_mm <- heatmap_height_mm / nrow(matrix)
  row_label_fontsize <- min(14, max(9, floor(cell_height_mm * 2.8)))
  column_label_fontsize <- min(14, max(9, floor(cell_width_mm * 2.8)))

  if (!is.null(group_annotation)) {
    # Figure 1d's style: no per-column colour bar (the column labels already
    # name every group), a level1 bar and a categorical bar above it, and a
    # legend for each. `group_palette` is unused in this branch.
    if (is.null(group_level1) || is.null(level1_palette) || is.null(group_annotation_palette)) {
      stop("A grouped column annotation needs level1 annotations and both palettes.")
    }
    group_level1 <- group_level1[colnames(matrix)]
    group_annotation <- group_annotation[colnames(matrix)]
    if (anyNA(group_level1) || anyNA(group_annotation)) {
      stop("Every column needs a level1 and a group annotation.")
    }
    level1_present <- intersect(names(level1_palette), unique(group_level1))
    group_present <- intersect(names(group_annotation_palette), unique(group_annotation))
    column_annotation <- ComplexHeatmap::HeatmapAnnotation(
      `Cell Type` = factor(group_level1, levels = level1_present),
      `Level2 Group` = factor(group_annotation, levels = group_present),
      col = list(
        `Cell Type` = level1_palette[level1_present],
        `Level2 Group` = group_annotation_palette[group_present]
      ),
      show_legend = TRUE,
      annotation_name_side = "left",
      annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold"),
      simple_anno_size = grid::unit(4, "mm"),
      annotation_legend_param = list(
        `Cell Type` = list(
          title = "Cell Type",
          title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
          labels_gp = grid::gpar(fontsize = 10)
        ),
        `Level2 Group` = list(
          title = "Level2 Group",
          title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
          labels_gp = grid::gpar(fontsize = 10)
        )
      )
    )
  } else if (is.null(group_level1)) {
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

  heatmap <- ComplexHeatmap::Heatmap(
    matrix,
    name = "Row-centered mean loading",
    col = color_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = row_order,
    column_order = column_order,
    top_annotation = column_annotation,
    column_title = paste0(
      "Row-centered GP mean loading: ", group_label, "\n", order_description
    ),
    column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
    row_title = "GP",
    row_title_gp = grid::gpar(fontsize = 12),
    row_names_gp = grid::gpar(fontsize = row_label_fontsize),
    column_names_gp = grid::gpar(fontsize = column_label_fontsize),
    column_names_rot = 90,
    heatmap_legend_param = list(
      title = "Row-centered mean loading",
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
  ComplexHeatmap::draw(
    heatmap,
    heatmap_legend_side = "right",
    merge_legend = TRUE,
    padding = grid::unit(c(8, 8, 8, 8), "mm")
  )
  grDevices::dev.off()
}
