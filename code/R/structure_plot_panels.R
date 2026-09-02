# Shared definitions for the per-lineage GP structure plots.
#
# Extended Data Figure 5 is one stacked figure: seven rows, one structure plot
# per T cell lineage, over the GPs that reach AUC > 0.9 for at least one of that
# lineage's sub-lineage clusters. Two scripts need the same row map, the same
# selection rule and the same palette:
#
#   script/FigureS5.R                    draws the seven rows, stacks them, and
#                                        records what it drew
#   script/verify_structure_plot_gps.R   re-derives the GP sets from the
#                                        published Extended Data Table 6 and
#                                        fails if they disagree with that
#                                        record or with the page's caption
#
# Both read them from here rather than keeping a copy, so the figure and the
# check cannot drift apart, and the caption's per-row GP counts have exactly
# one source.

# --- doc:panels ---
# Row letter per lineage, in stacking order -- the level-1 lineage order of
# Figure 1D and Figure 3, with thymocytes and DP cells excluded as they are
# there. The letters are the panel labels cowplot::plot_grid() draws on the
# assembled figure. Both consumers iterate over the names of this map, so a
# lineage cannot be drawn, or checked, under another lineage's letter.
structure_plot_panels <- c(
  "CD8" = "a",
  "CD4" = "b",
  "Treg" = "c",
  "gdT" = "d",
  "CD8aa" = "e",
  "Tz" = "f",
  "DN" = "g"
)

# A GP is shown in a lineage's panel when its one-vs-rest AUC for predicting
# membership in at least one of that lineage's annotation_level2 clusters
# exceeds this threshold. Those AUCs are the values published as Extended Data
# Table 6 -- script/ExtendedDataTable6_GP_AUC_cluster.R reads the same file --
# which is what makes the two checkable against each other.
structure_plot_auc_threshold <- 0.9

# Display filters. These apply to the plotted cells only, never to the
# AUC-based GP selection above, which always uses every cluster of the lineage.
structure_plot_min_cluster_cells <- 100      # smaller clusters are not drawn
structure_plot_max_cells_per_cluster <- 2000 # larger clusters are subsampled

# Assembled geometry: each row is drawn wide and short, and the seven are
# stacked into one figure. A row carries up to 21 cluster blocks and their
# labels, so it needs the width; the height then follows from keeping seven
# rows on one page.
structure_plot_width <- 16     # inches, the whole figure
structure_plot_row_height <- 3 # inches per lineage row
# --- doc:end-panels ---

# annotation_level2 -> annotation_level1, one entry per cluster.
#
# Stops if a cluster maps to more than one lineage: the panels are built by
# splitting clusters between lineages, so an ambiguous cluster would land in
# one panel and quietly inflate that panel's GP set.
cluster_lineage_map <- function(level2, level1) {
  if (length(level2) != length(level1)) {
    stop("cluster_lineage_map(): level2 and level1 must be the same length.")
  }
  pairs <- unique(data.frame(
    cluster = as.character(level2),
    lineage = as.character(level1),
    stringsAsFactors = FALSE
  ))
  pairs <- pairs[!is.na(pairs$cluster) & nzchar(pairs$cluster), ]
  if (anyDuplicated(pairs$cluster)) {
    ambiguous <- unique(pairs$cluster[duplicated(pairs$cluster)])
    stop(sprintf(
      "cluster_lineage_map(): these clusters map to more than one lineage: %s",
      paste(ambiguous, collapse = ", ")
    ))
  }
  stats::setNames(pairs$lineage, pairs$cluster)
}

# GPs that pass the AUC threshold in each lineage.
#
# `auc` is a clusters x GPs matrix (rows named by annotation_level2, columns by
# GP), `cluster_lineage` the map above. Returns one character vector of GP
# names per lineage, in the column order of `auc` -- ascending GP number, which
# is the order the rows' legends use.
gps_above_auc_by_lineage <- function(auc,
                                     cluster_lineage,
                                     lineages = names(structure_plot_panels),
                                     threshold = structure_plot_auc_threshold) {
  auc <- as.matrix(auc)
  if (is.null(rownames(auc)) || is.null(colnames(auc))) {
    stop("gps_above_auc_by_lineage(): `auc` needs cluster row names and GP column names.")
  }
  unmapped <- setdiff(rownames(auc), names(cluster_lineage))
  if (length(unmapped)) {
    stop(sprintf(
      "gps_above_auc_by_lineage(): no lineage for these clusters: %s",
      paste(unmapped, collapse = ", ")
    ))
  }
  row_lineage <- cluster_lineage[rownames(auc)]
  missing_lineages <- setdiff(lineages, row_lineage)
  if (length(missing_lineages)) {
    stop(sprintf(
      "gps_above_auc_by_lineage(): no clusters for these lineages: %s",
      paste(missing_lineages, collapse = ", ")
    ))
  }

  out <- lapply(lineages, function(lineage) {
    auc_lineage <- auc[row_lineage == lineage, , drop = FALSE]
    passes <- apply(auc_lineage, 2, function(x) any(x > threshold, na.rm = TRUE))
    colnames(auc_lineage)[passes]
  })
  stats::setNames(out, lineages)
}

# One color per GP, held fixed across rows.
#
# Assigning per row would give the same GP a different color in each lineage it
# appears in; this takes the union up front so GP2, shown in both the CD8aa and
# DN rows, reads as the same program in both.
#
# With up to 44 GPs in one row, no exported qualitative palette has enough
# distinct colors: pals::glasbey() stops at 32 and pals::polychrome() at 36.
# fastTopics's own 256-color glasbey -- the fallback its structure_plot() uses
# for >= 22 topics -- is large enough but unexported, so it is fetched by name:
# a future version that renames it fails here, loudly, instead of silently
# recoloring the figure. The assigned subset is then re-sorted by hue, and
# assigned in ascending GP order, so neighboring legend entries step through the
# hue wheel rather than jumping around it.
structure_plot_gp_colors <- function(gps) {
  gps <- unique(as.character(gps))
  gp_number <- suppressWarnings(as.integer(sub("^GP", "", gps)))
  if (anyNA(gp_number)) {
    stop("structure_plot_gp_colors(): every GP must be named GP<number>.")
  }
  gps <- gps[order(gp_number)]

  if (!requireNamespace("fastTopics", quietly = TRUE)) {
    stop("Package 'fastTopics' is required for the GP palette.")
  }
  glasbey <- utils::getFromNamespace("glasbey", "fastTopics")
  palette <- glasbey()[-1] # entry 1 is white: unusable as a bar fill
  if (length(gps) > length(palette)) {
    stop(sprintf(
      "structure_plot_gp_colors(): %d GPs but only %d palette colors.",
      length(gps), length(palette)
    ))
  }

  assigned <- palette[seq_along(gps)]
  hue <- grDevices::rgb2hsv(grDevices::col2rgb(assigned))["h", ]
  stats::setNames(assigned[order(hue)], gps)
}

# structure_plot()'s ggplot_call, with the bars rasterized.
#
# Same call signature and the same layers as fastTopics's default; the only
# change is ggrastr::rasterise() around geom_col(). A panel draws one bar per
# cell, so a vector PDF of the gdT panel is tens of thousands of rectangles --
# slow to open and larger than the rest of the figure put together. The axes,
# labels and legend stay vector.
rasterized_structure_plot_call <- function(dat, colors, ticks = NULL,
                                           font.size = 9, linewidth = 0) {
  ggplot2::ggplot(dat, ggplot2::aes(x = sample, y = prop, fill = topic)) +
    ggrastr::rasterise(
      ggplot2::geom_col(linewidth = linewidth, width = 1),
      dpi = 300
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, max(dat$sample) + 1),
      breaks = ticks,
      labels = names(ticks)
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = "", y = "topic proportion") +
    cowplot::theme_cowplot(font.size) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
