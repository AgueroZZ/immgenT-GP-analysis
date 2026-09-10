# Figure S2, panel d. GP activity across T cell clusters.
#
# One panel (see analysis/FigureS2.Rmd for the caption text):
#   S2D  Row-centered mean GP activity across the 99 level-2 clusters, all
#        200 GPs, in healthy non-thymocyte cells.
#
# Panels S2A-S2C come from script/FigureS2.R; this panel is kept in its own
# script because it shares nothing with them -- a different reference cell
# set, a different rendering module -- and re-running it means loading the
# full cell-by-GP loading matrix.
#
# --- internal ---
# Became Extended Data Figure 2d on 2026-09-10, folded in from the
# single-panel Extended Data Figure 3 it had been since 2026-09-09. Before
# that it was the retired Figure S4's panel s4b, promoted to a figure of its
# own when the two-panel "GP activity across tissues and T cell clusters"
# figure was split up (the tissue half briefly shipped as Figure S5 before
# being pulled back out of Extended Data the same day -- it is now the
# internal draft in experiments/fig_n5/). The PDF had been carried over
# byte-identically through all of that; the 2026-09-10 move is the first time
# it was re-rendered, on Ziang's call to drop the miniverse clusters and to
# annotate the columns the way Figure 1d does.
# --- end internal ---
#
# For each GP, its mean loading across clusters is subtracted from every cluster
# mean, so the panel shows where a program is more or less active than its own
# average rather than how large its loading is. The centered color scale is
# fixed at [-0.2, 0.2]; values outside this range saturate at the endpoint
# colors. Level2 columns follow Figure 1's level1 order, with level2 labels
# alphabetized within each level1 block. The columns are annotated the way
# Figure 1d annotates its cells -- a level1 bar and an annotation_level2_group
# bar, each with a legend, and no per-cluster colour bar, since the column
# labels already name every cluster.
#
# Required inputs (data/) -- see code/README.md's "Data provenance" table
# for the full picture:
#   L_pm_filtered.rds                        [code/pipeline/01b_filter_cells.R]
#   igt1_96_..._ADTonly.Rds                  [primary input Seurat object]

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
source("code/R/centered_mean_heatmap.R")
source("code/R/level2_group_palette.R")

figure_path <- "figures/final-selected/Figure S2"
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Setup, centering, and ordering
# ============================================================
gp_data <- load_gp_data()
reference <- healthy_nonthymocyte_reference(gp_data)

# Drop the miniverse clusters, as Figure 1d does: they are the ".wM" clusters,
# eight of them here (Figure 1d sees seven, because it drops DP as well), so
# 107 level2 clusters become 99. The exclusion is applied to the cells before
# any mean is taken, so the centering is over the columns that remain.
keep <- !reference$meta$annotation_level2_group %in% EXCLUDE_LEVEL2_GROUPS
if (anyNA(keep) || !any(keep)) {
  stop("annotation_level2_group is missing for some healthy non-thymocyte cells.")
}
L_reference <- reference$L[keep, , drop = FALSE]
meta_reference <- reference$meta[keep, , drop = FALSE]

centered_color_limit <- 0.2
level1_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN", "DP")

level2_result <- mean_loading_by_group(L_reference, meta_reference$annotation_level2)
level2_raw <- level2_result$matrix
level2_centered <- center_by_gp_mean(level2_raw)

level2_group_level1 <- level2_to_level1_map(
  meta_reference, colnames(level2_raw), level1_order
)
level2_group_group <- level2_to_group_map(meta_reference, colnames(level2_raw))
# Columns follow Figure 1d exactly: lineage, then annotation_level2_group as a
# contiguous block, then cluster alphabetically -- so each lineage's ".P"
# cluster sits at the end of its lineage rather than mid-alphabet. GP rows then
# follow the columns, in dominant-cluster blocks.
level2_order <- dominant_group_order(
  level2_raw,
  level2_group_block_order(
    colnames(level2_raw), level2_group_level1, level2_group_group, level1_order
  )
)

stopifnot(
  nrow(level2_centered) == 200L,
  ncol(level2_centered) == 99L,
  !any(level2_group_group %in% EXCLUDE_LEVEL2_GROUPS),
  max(abs(rowMeans(level2_centered))) < 1e-12
)

level1_palette <- ZemmourLib::immgent_colors$level1[level1_order]

# --- doc:rendering ---
render_centered_heatmap(
  level2_centered,
  NULL,
  "cluster (annotation_level2)",
  file.path(figure_path, "S2D.pdf"),
  level2_order$row_order,
  level2_order$column_order,
  centered_color_limit,
  paste0(
    "all 200 GPs; level2 columns: Figure 1 level1 order ",
    "(CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP); ",
    "level2_group blocks, alphabetical within block; ",
    "miniverse (.wM) clusters excluded; ",
    "GP rows: dominant-group blocks"
  ),
  group_level1 = level2_group_level1,
  level1_palette = level1_palette,
  group_annotation = level2_group_group,
  group_annotation_palette = LEVEL2_GROUP_COLORS[LEVEL2_GROUP_ORDER]
)

summary_dir <- "output/FigureS2d"   # build intermediate (not a manuscript panel)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  data.frame(
    panel = "S2D",
    grouping = "annotation_level2",
    view = "full row-centered mean loading",
    gp_count = nrow(level2_centered),
    group_count = ncol(level2_centered),
    excluded_level2_groups = paste(EXCLUDE_LEVEL2_GROUPS, collapse = ";"),
    centered_definition = "group mean minus mean across groups for each GP",
    color_min = -centered_color_limit,
    color_mid = 0,
    color_max = centered_color_limit,
    observed_min = min(level2_centered),
    observed_max = max(level2_centered)
  ),
  file.path(summary_dir, "S2D_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

message("Wrote Figure S2 panel d to ", normalizePath(figure_path))
