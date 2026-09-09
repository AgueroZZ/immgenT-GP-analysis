# Figure S3. GP activity across T cell clusters.
#
# One panel (see analysis/FigureS3.Rmd for the caption text):
#   s3   Row-centered mean GP activity across the 107 level-2 clusters, all
#        200 GPs, in healthy non-thymocyte cells.
#
# --- internal ---
# New Extended Data Figure 3 on 2026-09-09. The panel itself is not new: it is
# the retired Figure S4's panel s4b, unchanged, promoted to a figure of its own
# when the two-panel "GP activity across tissues and T cell clusters" figure was
# split up -- the cluster half here, and a 31-GP tissue half as Figure S5. The
# PDF was carried over byte-identically rather than re-rendered; this script is
# the old FigureS4.R with its tissue half removed and its shared helpers moved
# into code/R/centered_mean_heatmap.R.
# --- end internal ---
#
# For each GP, its mean loading across clusters is subtracted from every cluster
# mean, so the panel shows where a program is more or less active than its own
# average rather than how large its loading is. The centered color scale is
# fixed at [-0.2, 0.2]; values outside this range saturate at the endpoint
# colors. Level2 columns follow Figure 1's level1 order, with level2 labels
# alphabetized within each level1 block.
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

figure_path <- "figures/final-selected/Figure S3"
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Setup, centering, and ordering
# ============================================================
gp_data <- load_gp_data()
reference <- healthy_nonthymocyte_reference(gp_data)
L_reference <- reference$L
meta_reference <- reference$meta

centered_color_limit <- 0.2
level1_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN", "DP")

level2_result <- mean_loading_by_group(L_reference, meta_reference$annotation_level2)
level2_raw <- level2_result$matrix
level2_centered <- center_by_gp_mean(level2_raw)

level2_group_level1 <- level2_to_level1_map(
  meta_reference, colnames(level2_raw), level1_order
)
level2_order <- dominant_group_order(
  level2_raw,
  level2_column_order(colnames(level2_raw), level2_group_level1, level1_order)
)

stopifnot(
  nrow(level2_centered) == 200L,
  ncol(level2_centered) == 107L,
  max(abs(rowMeans(level2_centered))) < 1e-12
)

level2_palette <- palette_for_groups(
  colnames(level2_centered),
  ZemmourLib::immgent_colors$level2,
  "annotation_level2"
)
level1_palette <- ZemmourLib::immgent_colors$level1[level1_order]

# --- doc:rendering ---
render_centered_heatmap(
  level2_centered,
  level2_palette,
  "cluster (annotation_level2)",
  file.path(figure_path, "s3.pdf"),
  level2_order$row_order,
  level2_order$column_order,
  centered_color_limit,
  paste0(
    "all 200 GPs; level2 columns: Figure 1 level1 order ",
    "(CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP); ",
    "alphabetical within level1; GP rows: dominant-group blocks"
  ),
  group_level1 = level2_group_level1,
  level1_palette = level1_palette
)

summary_dir <- "output/FigureS3"   # build intermediate (not a manuscript panel)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  data.frame(
    panel = "s3",
    grouping = "annotation_level2",
    view = "full row-centered mean loading",
    gp_count = nrow(level2_centered),
    group_count = ncol(level2_centered),
    centered_definition = "group mean minus mean across groups for each GP",
    color_min = -centered_color_limit,
    color_mid = 0,
    color_max = centered_color_limit,
    observed_min = min(level2_centered),
    observed_max = max(level2_centered)
  ),
  file.path(summary_dir, "S3_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

message("Wrote Figure S3 to ", normalizePath(figure_path))
