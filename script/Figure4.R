# Figure 4. GP activity across level-2 clusters.
#
# Panels produced:
#   4a  Max AUC (level-1 lineage) vs max AUC (level-2 cluster) scatter, per GP.
#   4b  The 69 cluster-associated GPs across the level-2 clusters: row-centered
#       mean-loading heatmap, columns grouped by level 1.
#   4c  The Treg row of the per-lineage structure plots as stacked bars: one
#       bar per cluster, holding that cluster's mean loading per GP.
#   4d  The Treg structure plot itself, one bar per cell.
#
# 4b, 4c and 4d show the GPs and the cells of the per-lineage structure plots,
# so the row map, the AUC > 0.9 rule, the display filters and the palette come
# from code/R/structure_plot_panels.R, and 4b's rendering from
# code/R/centered_mean_heatmap.R -- the module Extended Data Figure 2d uses --
# rather than being restated here. The checks at the bottom compare what this
# script drew against output/structure_plot_record/, written by
# script/structure_plot_record.R.
#
# --- internal ---
# New figure, inserted at 4 on 2026-09-10 (the activation, tissue and CITE-seq
# figures moved to 5, 6, 7 and the RQVI figure to 8). It has no published
# counterpart, so verify_panels.sh has no pair for it and the captions on
# analysis/Figure4.Rmd are ours, not the published legends every other page
# carries.
#
# Ported from experiments/GP-level2/, which is retired: this script is the only
# copy. Decisions made there and kept here: 4a is labelled with GP names, not
# Figure 6a's top-three-categories blocks (53 of 200 GPs clear 0.9 here, 52 of
# them above the diagonal); 4c's bars are raw mean loadings, not renormalized,
# because that is what makes them the average of a structure plot's bars; 4c is
# drawn tall and narrow rather than at the structure plot's shape; and every
# lineage's miniverse (.wM) cluster is dropped, as in Figure 1d.
# --- end internal ---
#
# Required inputs (data/) -- see code/README.md's "Data provenance" table:
#   igt1_96_..._ADTonly.Rds                           [primary input Seurat object]
#   L_pm_filtered.rds                                 [code/pipeline/01b_filter_cells.R]
#   level_1_AUC_list_figure_no_thymocytes_healthy.rds [code/pipeline/02_compute_auc.R]
#   level_2_AUC_list_figure_no_thymocytes_healthy.rds [code/pipeline/02_compute_auc.R]

# --- doc:setup ---
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggrastr)
  library(cowplot)
  library(dplyr)
  library(fastTopics)      # structure_plot()
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ZemmourLib)      # immgent_colors
})

if (!file.exists("code/R/structure_plot_panels.R")) {
  stop("Run this script from the immgenT-GP-analysis repository root.")
}
source("code/R/structure_plot_panels.R")   # the per-lineage structure plots' rows, GP rule, palette
source("code/R/centered_mean_heatmap.R")   # Extended Data Figure 2d's heatmap rendering
source("code/R/level2_group_palette.R")    # EXCLUDE_LEVEL2_GROUPS, as Figure 1d and ED 2d use

data_path   <- "data/"
figure_path <- "figures/final-selected/Figure 4/"
record_path <- "output/Figure4/"
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)
dir.create(record_path, recursive = TRUE, showWarnings = FALSE)

# Figure scripts here have been seen to exit 0 having written nothing, so every
# PDF is checked against this timestamp at the end rather than against the exit
# code (script/README.md, "A re-run can silently not write").
run_started_at <- Sys.time()

# ============================================================
# Load data
# ============================================================
seurat_meta <- readRDS(paste0(
  data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
colnames(L_pm_filtered) <- gsub("^K", "GP", colnames(L_pm_filtered))
if (!all(grepl("^GP[0-9]+$", colnames(L_pm_filtered)))) {
  stop("The loading matrix does not carry K##/GP## column names.")
}

level_1_AUC_list <- readRDS(paste0(
  data_path, "level_1_AUC_list_figure_no_thymocytes_healthy.rds"
))
level_2_AUC_list <- readRDS(paste0(
  data_path, "level_2_AUC_list_figure_no_thymocytes_healthy.rds"
))
colnames(level_1_AUC_list$auc) <- gsub("^K", "GP", colnames(level_1_AUC_list$auc))
colnames(level_2_AUC_list$auc) <- gsub("^K", "GP", colnames(level_2_AUC_list$auc))

# Both AUC matrices index GPs by names taken from the loading matrix, so they
# have to be talking about the same GPs in the same order.
if (!identical(colnames(level_1_AUC_list$auc), colnames(L_pm_filtered)) ||
      !identical(colnames(level_2_AUC_list$auc), colnames(L_pm_filtered))) {
  stop("The AUC matrices and L_pm_filtered do not carry the same GP columns in the same order.")
}

level1_all <- seurat_meta_filtered$annotation_level1
level2_all <- seurat_meta_filtered$annotation_level2
healthy_non_thymocyte <- which(
  seurat_meta_filtered$condition_broad == "healthy" & level1_all != "thymocyte"
)
meta_reference <- seurat_meta_filtered[healthy_non_thymocyte, , drop = FALSE]
L_healthy <- L_pm_filtered[healthy_non_thymocyte, , drop = FALSE]
overall_mean <- colMeans(L_healthy, na.rm = TRUE)

message(sprintf(
  "%d healthy non-thymocyte cells, %d GPs, %d level-2 clusters",
  nrow(L_healthy), ncol(L_healthy), dplyr::n_distinct(meta_reference$annotation_level2)
))

# ============================================================
# 4a: Max AUC level-1 vs max AUC level-2
# ============================================================
# Figure 6a and 6d's construction, with their category-count filters kept as
# they are there: level-1 lineages need 1000 cells, level-2 clusters 100.
# --- internal ---
# Ported from script/Figure6.R's "6a" and "6d prep" sections. The only change
# is which two of the three maxima are plotted against each other.
# --- end internal ---
level_1_small <- table(meta_reference$annotation_level1)
level_1_small <- names(level_1_small[level_1_small < 1000])
level_1_AUC <- level_1_AUC_list$auc[
  !rownames(level_1_AUC_list$auc) %in% level_1_small, , drop = FALSE
]

level_2_small <- table(meta_reference$annotation_level2)
level_2_small <- names(level_2_small[level_2_small < 100])
level_2_AUC <- level_2_AUC_list$auc[
  !rownames(level_2_AUC_list$auc) %in% level_2_small, , drop = FALSE
]

# A high AUC only means "GP predicts this category" when the category's mean
# loading is above the overall mean; the other direction is a GP the category
# lacks. Same masking as Figure 6.
category_mean <- function(labels, categories) {
  t(vapply(categories, function(cat) {
    colMeans(L_healthy[labels == cat, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(L_healthy))))
}
level_1_positive <- sweep(
  category_mean(meta_reference$annotation_level1, rownames(level_1_AUC)),
  2, overall_mean, "-"
) > 0
level_2_positive <- sweep(
  category_mean(meta_reference$annotation_level2, rownames(level_2_AUC)),
  2, overall_mean, "-"
) > 0

masked_max <- function(auc, positive) {
  auc[!positive] <- NA
  list(
    value = apply(auc, 2, max, na.rm = TRUE),
    name = apply(auc, 2, function(x) rownames(auc)[which.max(x)])
  )
}
l1_max <- masked_max(level_1_AUC, level_1_positive)
l2_max <- masked_max(level_2_AUC, level_2_positive)

df_a <- data.frame(
  Factor = colnames(level_1_AUC),
  annotation_Level1 = l1_max$name,
  annotation_Level2 = l2_max$name[colnames(level_1_AUC)],
  Max_AUC_Level1 = l1_max$value,
  Max_AUC_Level2 = l2_max$value[colnames(level_1_AUC)],
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(residual = Max_AUC_Level2 - Max_AUC_Level1)

# Figure 6a's highlight rule verbatim: AUC > 0.9 on either axis, coloured by
# which of the two maxima is the larger. Only the labelling differs -- GP names,
# not each point's top categories.
highlighted <- df_a |>
  dplyr::filter(is.finite(residual), Max_AUC_Level1 > 0.9 | Max_AUC_Level2 > 0.9) |>
  dplyr::pull(Factor)

label_above <- df_a |>
  dplyr::filter(Factor %in% highlighted, residual > 0) |>
  dplyr::mutate(
    label_text = Factor
  )
label_below <- df_a |>
  dplyr::filter(Factor %in% highlighted, residual <= 0) |>
  dplyr::mutate(
    label_text = Factor
  )
message(sprintf(
  "a: %d of %d GPs highlighted (%d above the diagonal, %d on or below)",
  length(highlighted), nrow(df_a), nrow(label_above), nrow(label_below)
))

axis_limits <- function(x, pad = 0.04) {
  x <- x[is.finite(x)]
  c(min(x) - pad, max(x) + pad)
}

# --- internal ---
# Figure 6a labels each highlighted point with its top three categories. Here 53
# of 200 GPs clear 0.9 and 52 do it above the diagonal, which at that label size
# is unreadable and costs the points ggrepel cannot place, so Ziang's call on
# 2026-09-09 was GP names only. Which categories each maximum is attained in is
# in records/a_max_auc_level1_level2.csv.
# --- end internal ---
p_4a <- ggplot(df_a, aes(Max_AUC_Level1, Max_AUC_Level2)) +
  geom_point(alpha = 0.3, size = 1.8) +
  geom_point(data = label_above, color = "#1f78b4", alpha = 0.8, size = 1.8) +
  geom_point(data = label_below, color = "#e31a1c", alpha = 0.8, size = 1.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(
    xlim = axis_limits(df_a$Max_AUC_Level1),
    ylim = axis_limits(df_a$Max_AUC_Level2),
    expand = FALSE, clip = "off"
  ) +
  labs(
    x = "Max AUC (Level-1)", y = "Max AUC (Level-2)",
    title = "Max AUC: Level-1 vs Level-2"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.margin = margin(10, 40, 10, 40)) +
  geom_text_repel(
    seed = 42, data = label_above, aes(label = label_text), color = "#1f78b4",
    size = 2.5, lineheight = 0.85, direction = "y", nudge_x = -0.035,
    segment.color = "#1f78b4",
    arrow = arrow(length = unit(0.008, "npc"), type = "closed", angle = 20),
    force = 3, force_pull = 0.1, box.padding = 0.4, point.padding = 0.15,
    max.time = 10, max.iter = 2e4, max.overlaps = 20,
    min.segment.length = 0.01, segment.alpha = 0.7
  ) +
  geom_text_repel(
    seed = 42, data = label_below, aes(label = label_text), color = "#e31a1c",
    size = 2.5, lineheight = 0.85, direction = "y", nudge_x = 0.035,
    segment.color = "#e31a1c",
    arrow = arrow(length = unit(0.008, "npc"), type = "closed", angle = 20),
    force = 3, force_pull = 0.1, box.padding = 0.4, point.padding = 0.15,
    max.time = 10, max.iter = 2e4, max.overlaps = 20,
    min.segment.length = 0.01, segment.alpha = 0.7
  )
ggsave(paste0(figure_path, "4a.pdf"), p_4a, width = 8, height = 8, dpi = 300)

write.csv(
  df_a[order(-df_a$Max_AUC_Level2), c(
    "Factor", "Max_AUC_Level1", "annotation_Level1",
    "Max_AUC_Level2", "annotation_Level2", "residual"
  )],
  file.path(record_path, "4a_max_auc_level1_level2.csv"),
  row.names = FALSE
)

# ============================================================
# Shared with 4b, 4c and 4d: the per-lineage structure plots' GP sets and cells
# ============================================================
cluster_lineage <- cluster_lineage_map(level2_all, level1_all)
auc_level2_full <- level_2_AUC_list$auc
auc_clusters <- rownames(auc_level2_full)
prefix_lineage <- sub("[.].*$", "", auc_clusters)

# GP selection uses every cluster of the lineage, as it does in structure_plot_record.R --
# the 100-cell filter below is a display filter only.
panel_gps <- gps_above_auc_by_lineage(auc_level2_full, cluster_lineage)
gp_union <- unique(unlist(panel_gps, use.names = FALSE))
gp_union <- gp_union[order(as.integer(sub("^GP", "", gp_union)))]

message(sprintf(
  "%d GPs over %d rows (AUC > %.1f): %s",
  length(gp_union), length(panel_gps), structure_plot_auc_threshold,
  paste(sprintf("%s %d", names(panel_gps), lengths(panel_gps)), collapse = ", ")
))

# Cells and clusters drawn per lineage: the healthy non-thymocyte cells of the
# lineage, in clusters of at least structure_plot_min_cluster_cells cells.
# --- internal ---
# Unlike the structure plots these means are taken over every cell of the cluster, with
# no 2000-cell cap: the cap is there so one large cluster cannot crowd out the
# rest of a structure plot's width, and a bar of means has no width to crowd.
# The cluster *set* is the same, which is what the check at the bottom compares.
# --- end internal ---
lineage_cells_drawn <- list()
for (lineage in names(structure_plot_panels)) {
  cells <- healthy_non_thymocyte[level1_all[healthy_non_thymocyte] == lineage]
  size <- table(droplevels(factor(level2_all[cells])))
  keep <- names(size)[size >= structure_plot_min_cluster_cells]
  lineage_cells_drawn[[lineage]] <- cells[level2_all[cells] %in% keep]
}

# Every panel here drops the miniverse clusters, as Figure 1d and Extended Data
# Figure 2d do -- the exclusion comes from their EXCLUDE_LEVEL2_GROUPS rather
# than from a rule of this figure's own, so the four heatmaps cannot disagree
# about what is shown. In this cell set that is seven clusters: CD4.wM 10543
# cells, CD8.wM 6174, gdT.wM 1732, DN.wM 822, Tz.wM 706, Treg.wM 654,
# CD8aa.wM 247.
# --- internal ---
# 4c and 4d dropped them from 2026-09-09 (Ziang's request, then stated as
# temporary) but 4b did not, which left this figure disagreeing with itself and
# with Figure 1d. Made uniform on 2026-09-10, on Ziang's call, by switching to
# the shared constant.
# --- end internal ---
level2_group_all <- seurat_meta_filtered$annotation_level2_group
if (anyNA(level2_group_all[healthy_non_thymocyte])) {
  stop("annotation_level2_group is missing for some healthy non-thymocyte cells.")
}
drop_miniverse <- function(cells) {
  cells[!level2_group_all[cells] %in% EXCLUDE_LEVEL2_GROUPS]
}

# The group label and the ".wM" cluster names have to pick out the same
# clusters, or the caption and the panel would describe different sets.
excluded_clusters <- sort(unique(as.character(
  level2_all[healthy_non_thymocyte][
    level2_group_all[healthy_non_thymocyte] %in% EXCLUDE_LEVEL2_GROUPS
  ]
)))
if (!all(grepl("[.]wM$", excluded_clusters))) {
  stop(sprintf(
    "annotation_level2_group '%s' is not exactly the .wM clusters: %s",
    paste(EXCLUDE_LEVEL2_GROUPS, collapse = ", "), paste(excluded_clusters, collapse = ", ")
  ))
}
message(sprintf("dropping %d miniverse clusters: %s",
                length(excluded_clusters), paste(excluded_clusters, collapse = ", ")))

lineage_cells_drawn <- lapply(lineage_cells_drawn, drop_miniverse)

# Mean loading per cluster, over the GPs of that cluster's lineage row.
mean_loading_by_cluster <- function(cells, gps) {
  labels <- droplevels(factor(level2_all[cells]))
  sums <- rowsum(L_pm_filtered[cells, gps, drop = FALSE], group = labels, reorder = TRUE)
  counts <- as.integer(table(labels)[rownames(sums)])
  list(matrix = sweep(sums, 1L, counts, "/"), counts = counts)
}

# ============================================================
# 4b: the same GPs as a heatmap, level-2 columns grouped by level 1
# ============================================================
# --- internal ---
# The centering, the fixed [-0.2, 0.2] scale, the dominant-group row order and
# the two column annotation bars (Cell Type + Level2 Group, one legend each)
# and the lineage/level2_group/cluster column order are Extended Data Figure
# 2d's -- and through it Figure 1d's -- called from
# code/R/centered_mean_heatmap.R -- the module that figure itself uses -- rather
# than copied, so this panel cannot drift from it. It is not a subset of that
# panel's matrix: the columns here are the clusters the structure plots draw
# (>= 100 healthy non-thymocyte cells, no DP, no thymocytes), and each GP is
# centered on its mean across those columns. It keeps the .wM clusters that c
# and 4d drop.
# --- end internal ---
cells_drawn <- sort(unlist(lineage_cells_drawn, use.names = FALSE))
labels_drawn <- droplevels(factor(level2_all[cells_drawn]))
heat_means <- mean_loading_by_group(
  L_pm_filtered[cells_drawn, gp_union, drop = FALSE], labels_drawn
)
heat_raw <- heat_means$matrix          # GPs x clusters
heat_centered <- center_by_gp_mean(heat_raw)

# Figure 1's level-1 order, minus DP and thymocytes, which these rows exclude.
level1_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN")
level2_group_level1 <- level2_to_level1_map(
  meta_reference, colnames(heat_raw), level1_order
)
level2_group_group <- level2_to_group_map(meta_reference, colnames(heat_raw))

# Columns follow Figure 1d and Extended Data Figure 2d exactly: lineage, then
# annotation_level2_group as a contiguous block, then cluster alphabetically --
# so each lineage's ".P" cluster sits at the end of its lineage rather than
# mid-alphabet, and the level2_group bar reads as blocks. GP rows then follow
# the columns, in dominant-cluster blocks.
heat_order <- dominant_group_order(
  heat_raw,
  level2_group_block_order(
    colnames(heat_raw), level2_group_level1, level2_group_group, level1_order
  )
)

level1_palette <- ZemmourLib::immgent_colors$level1[level1_order]

centered_color_limit <- 0.2
render_centered_heatmap(
  heat_centered,
  NULL,
  "cluster (annotation_level2)",
  paste0(figure_path, "4b.pdf"),
  heat_order$row_order,
  heat_order$column_order,
  centered_color_limit,
  sprintf(
    paste0(
      "%d GPs with AUC > %.1f in some cluster; miniverse (.wM) clusters excluded\n",
      "level2 columns: level1 order (%s); level2_group blocks, alphabetical ",
      "within block; GP rows: dominant-cluster blocks"
    ),
    nrow(heat_centered), structure_plot_auc_threshold,
    paste(level1_order, collapse = ", ")
  ),
  group_level1 = level2_group_level1,
  level1_palette = level1_palette,
  group_annotation = level2_group_group,
  group_annotation_palette = LEVEL2_GROUP_COLORS[LEVEL2_GROUP_ORDER]
)

write.csv(
  data.frame(gp = rownames(heat_centered), heat_centered, check.names = FALSE),
  file.path(record_path, "4b_row_centered_mean_loading.csv"), row.names = FALSE
)
write.csv(
  data.frame(
    cluster = heat_means$counts$group,
    level1 = unname(level2_group_level1[heat_means$counts$group]),
    n_cells = heat_means$counts$n_cells
  ),
  file.path(record_path, "4b_column_cells.csv"), row.names = FALSE
)

# ============================================================
# 4c: one stacked bar of mean loadings per level-2 cluster, the Treg row
# ============================================================
# The Treg structure-plot row with one bar per cluster, holding that
# cluster's mean loading per GP instead of one bar per cell -- the average of
# that row's bars. Same GPs (the 11 with AUC > 0.9 in some Treg cluster), same
# per-row palette, same cells, minus Treg.wM.
#
# Bars are raw mean loadings, not rescaled to a common height: that is what
# makes them the average of those bars, since structure_plot() does not
# renormalize a subset of topics either. A tall bar means the cluster carries
# more total activity over the row's GPs, not just a different mix.
#
# --- internal ---
# Drawn for all seven lineages until 2026-09-09, when Ziang cut the panel to
# Treg to pair with d. The means are still computed for all seven, because they
# are the numbers behind 4b's columns and the check at the bottom compares all
# 1794 of them against that panel; only the Treg row is plotted. A normalized
# view (every cluster rescaled to sum to 1) was drawn while choosing and is not
# produced any more, but its numbers are the record's `proportion` column.
# --- end internal ---
#
# Unlike 4d and the structure plots this panel is not drawn at the structure
# plot's 16 x 3 in: a chart of seven bars does not need a structure plot's
# width, and Ziang asked for tall and narrow.
bar_figure_width <- 5.5  # inches
bar_figure_height <- 7

# Mean loading per cluster over the row's GPs, plus the record of every segment.
bar_means <- function(lineage) {
  gps_lineage <- panel_gps[[lineage]]
  means <- mean_loading_by_cluster(lineage_cells_drawn[[lineage]], gps_lineage)
  bar_matrix <- means$matrix
  totals <- rowSums(bar_matrix)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop(sprintf("row %s: a cluster has no loading at all.", lineage))
  }

  # This row's colours, assigned from the top of the palette without reference
  # to any other row -- the same per-row rule the structure plots use, so
  # colour means different GPs in different rows.
  colors_lineage <- structure_plot_row_colors(gps_lineage)
  if (!identical(names(colors_lineage), colnames(bar_matrix))) {
    stop(sprintf("row %s: palette order does not match its GP columns.", lineage))
  }

  long <- data.frame(
    cluster = factor(
      rep(rownames(bar_matrix), times = ncol(bar_matrix)),
      levels = rownames(bar_matrix)
    ),
    gp = factor(
      rep(colnames(bar_matrix), each = nrow(bar_matrix)),
      levels = colnames(bar_matrix)
    ),
    value = as.vector(bar_matrix),
    stringsAsFactors = FALSE
  )

  record <- data.frame(
    panel = structure_plot_panels[[lineage]],
    lineage = lineage,
    drawn_in_4c = lineage == bar_lineage,
    cluster = as.character(long$cluster),
    n_cells = means$counts[match(long$cluster, rownames(bar_matrix))],
    gp = as.character(long$gp),
    color = unname(colors_lineage[as.character(long$gp)]),
    mean_loading = long$value,
    # The same segment as a share of its cluster's total -- the normalized view's
    # number, kept so that variant can be redrawn from the record alone.
    proportion = long$value / totals[as.character(long$cluster)],
    stringsAsFactors = FALSE
  )

  list(long = long, colors = colors_lineage, record = record, gps = gps_lineage)
}

bar_lineage <- "Treg"
bar_data <- lapply(names(structure_plot_panels), bar_means)
names(bar_data) <- names(structure_plot_panels)

drawn <- bar_data[[bar_lineage]]
p_4c <- ggplot(drawn$long, aes(cluster, value, fill = gp)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = drawn$colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = "", y = "mean loading", fill = "",
    title = sprintf(
      "%s (%d GPs, AUC > %.1f)", bar_lineage, length(drawn$gps),
      structure_plot_auc_threshold
    )
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  cowplot::theme_cowplot(9) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 7)
  )
ggsave(
  paste0(figure_path, "4c.pdf"), p_4c,
  width = bar_figure_width, height = bar_figure_height, dpi = 300
)

# Every bar segment of all seven rows, raw and as a share of its cluster's
# total; drawn_in_4c marks the row this panel shows.
bar_record <- do.call(rbind, lapply(bar_data, `[[`, "record"))
write.csv(bar_record, file.path(record_path, "4c_mean_loading_by_cluster.csv"), row.names = FALSE)

# ============================================================
# 4d: the Treg structure-plot row on its own
# ============================================================
# structure_plot_record.R's loop body for one lineage, unchanged including its seeds and its
# 2000-cell cap, so this is that row and not a redrawing of it.
d_lineage <- "Treg"
gps_d <- panel_gps[[d_lineage]]
cells_d <- healthy_non_thymocyte[level1_all[healthy_non_thymocyte] == d_lineage]
cluster_size_d <- table(droplevels(factor(level2_all[cells_d])))
small_d <- names(cluster_size_d)[cluster_size_d < structure_plot_min_cluster_cells]
cells_d <- cells_d[!level2_all[cells_d] %in% small_d]
cells_d <- drop_miniverse(cells_d)   # as in 4b and 4c

set.seed(1234)
keep_d <- unlist(lapply(
  split(seq_along(cells_d), level2_all[cells_d]),
  function(idx) {
    if (length(idx) > structure_plot_max_cells_per_cluster) {
      sample(idx, structure_plot_max_cells_per_cluster)
    } else {
      idx
    }
  }
))
cells_d <- cells_d[keep_d]

fit_d <- L_pm_filtered[cells_d, gps_d, drop = FALSE]
colors_d <- structure_plot_row_colors(gps_d)
if (!identical(names(colors_d), colnames(fit_d))) {
  stop("4d: palette order does not match its GP columns.")
}

set.seed(1234)
p_4d <- structure_plot(
  fit_d, topics = gps_d, gap = 40, n = 10000, colors = colors_d,
  grouping = factor(level2_all[cells_d]),
  ggplot_call = rasterized_structure_plot_call
) +
  labs(
    y = "membership", color = "", fill = "",
    title = sprintf(
      "%s (%d GPs, AUC > %.1f)", d_lineage, length(gps_d), structure_plot_auc_threshold
    )
  ) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2)) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 5),
    legend.spacing.y = unit(0.02, "cm")
  )
ggsave(
  paste0(figure_path, "4d.pdf"), p_4d,
  width = structure_plot_width, height = structure_plot_row_height,
  dpi = 300, limitsize = FALSE
)

# ============================================================
# Do 4b, 4c and 4d show what the per-lineage structure plots show?
# ============================================================
# The point of these panels is to be those rows' GPs and cells in another form,
# so they are checked against what the rows recorded rather than against a
# second copy of the rule. The record is written by
# script/structure_plot_record.R; run it first if this read fails.
s4_gps <- read.csv("output/structure_plot_record/panel_gps.csv", stringsAsFactors = FALSE)
s4_clusters <- read.csv("output/structure_plot_record/panel_clusters.csv", stringsAsFactors = FALSE)

for (lineage in names(structure_plot_panels)) {
  recorded <- s4_gps[s4_gps$lineage == lineage, ]
  drawn_colors <- structure_plot_row_colors(panel_gps[[lineage]])
  if (!identical(recorded$gp, unname(names(drawn_colors))) ||
        !identical(recorded$color, unname(drawn_colors))) {
    stop(sprintf("row %s: GPs or colours differ from the structure-plot record.", lineage))
  }
  # Every panel here shows the clusters those rows drew, minus this lineage's
  # miniverse cluster and nothing else.
  recorded_clusters <- sort(s4_clusters$cluster[
    s4_clusters$lineage == lineage & s4_clusters$n_cells_drawn > 0
  ])
  expected_clusters <- setdiff(recorded_clusters, excluded_clusters)
  drawn_clusters <- sort(unique(as.character(level2_all[lineage_cells_drawn[[lineage]]])))
  if (!identical(drawn_clusters, expected_clusters) ||
        length(drawn_clusters) == length(recorded_clusters)) {
    stop(sprintf(
      "row %s: this figure dropped something other than the miniverse cluster.", lineage
    ))
  }
}

# 4b's rows are exactly the GPs c draws, and its uncentered entries are 4c's bar
# segments -- the two panels must not be able to show different numbers.
if (!identical(sort(rownames(heat_raw)), sort(unique(as.character(bar_record$gp))))) {
  stop("4b's GP rows are not the union of 4c's row GP sets.")
}
# 4c's clusters are a subset of 4b's columns, so every segment must be findable
# there; a missing one indexes to NA and fails this comparison.
c_vs_b <- max(abs(
  bar_record$mean_loading - heat_raw[cbind(bar_record$gp, bar_record$cluster)]
))
if (!is.finite(c_vs_b) || c_vs_b > 1e-12) {
  stop(sprintf("4c's bar segments and 4b's entries disagree (max |diff| = %g).", c_vs_b))
}

# The record's proportion column is 4c's bars rescaled per cluster, so each
# cluster's shares sum to 1.
cluster_totals <- tapply(bar_record$proportion, bar_record$cluster, sum)
if (max(abs(cluster_totals - 1)) > 1e-12) {
  stop(sprintf(
    "the record's shares do not sum to 1 per cluster (worst |sum - 1| = %g).",
    max(abs(cluster_totals - 1))
  ))
}
share_diff <- max(abs(
  bar_record$proportion -
    bar_record$mean_loading / ave(bar_record$mean_loading, bar_record$cluster, FUN = sum)
))
if (share_diff > 1e-12) {
  stop(sprintf("the record's shares are not 4c's bars rescaled (max |diff| = %g).", share_diff))
}

# 4d is that row minus Treg's miniverse cluster, so on the clusters it keeps it
# must have kept the same cells: same 100-cell filter, same 2000-cell cap, same
# seed, hence the same per-cluster counts.
d_drawn <- table(droplevels(factor(level2_all[cells_d])))
d_recorded <- s4_clusters[
  s4_clusters$lineage == d_lineage & s4_clusters$n_cells_drawn > 0 &
    !s4_clusters$cluster %in% excluded_clusters,
]
d_recorded <- d_recorded[order(d_recorded$cluster), ]
if (!identical(names(d_drawn), d_recorded$cluster) ||
      !identical(as.integer(d_drawn), as.integer(d_recorded$n_cells_drawn))) {
  stop("4d does not draw the same Treg cells as the structure-plot record's Treg row.")
}

# Did this run actually write the panels?
expected <- paste0(figure_path, c("4a.pdf", "4b.pdf", "4c.pdf", "4d.pdf"))
for (f in expected) {
  if (!file.exists(f) || file.mtime(f) < run_started_at || file.size(f) == 0) {
    stop(sprintf("%s was not written by this run (missing, empty, or older than the run)", f))
  }
}
message(sprintf(
  "wrote %s; 4b/4c/4d match output/structure_plot_record (4c vs 4b max |diff| = %g)",
  paste(basename(expected), collapse = ", "), c_vs_b
))
