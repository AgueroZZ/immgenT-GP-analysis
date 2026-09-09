# Figure S5. The tissue-associated GPs, across tissues and across lineages.
#
# Panels produced (see analysis/FigureS5.Rmd for the caption text):
#   s5a  Row-centered mean GP activity across the 18 tissues, restricted to the
#        31 GPs that are active across tissues. Blue-white-red scale.
#   s5b  The same 31 GPs in the same row order, but across the 8 T cell
#        lineages. Purple-white-green scale, so the two halves of the figure
#        cannot be mistaken for each other.
#
# --- internal ---
# New Extended Data Figure 5 on 2026-09-09. Panel s5a is the retired Figure S4's
# panel s4a with its rows filtered to the 31 GPs and the dominant-group order
# recomputed on that submatrix; s5b is new. The retired figure's other panel
# (the cluster heatmap) became Figure S3.
# --- end internal ---
#
# "Active across tissues" is the selection the manuscript's "31 programs active
# across tissues" refers to: a GP is kept when its *raw* (uncentered) mean
# loading reaches 0.1 in at least one of the 18 tissues. That is a threshold on
# activity, not on differential activity -- the panels then show, per GP, how
# that activity is distributed, by subtracting the GP's mean across the groups
# shown from every group mean. Each panel's centered color scale is fixed:
# [-0.2, 0.2] for (a), matching the retired 200-GP tissue heatmap it comes from,
# and the tighter [-0.1, 0.1] for (b), because spreading a program over eight
# lineages instead of eighteen tissues gives much smaller deviations and (a)'s
# scale renders the lineage panel almost blank. Values outside each range
# saturate at the endpoint colors: 0.9% of (a)'s cells and 3.6% of (b)'s.
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

figure_path <- "figures/final-selected/Figure S5"
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)
run_started_at <- Sys.time()

# ============================================================
# Setup: the healthy non-thymocyte tissue and lineage mean matrices
# ============================================================
gp_data <- load_gp_data()
reference <- healthy_nonthymocyte_reference(gp_data)
L_reference <- reference$L
meta_reference <- reference$meta

organ_color_limit <- 0.2
level1_color_limit <- 0.1
level1_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN", "DP")

organ_raw <- mean_loading_by_group(L_reference, meta_reference$organ_simplified)$matrix
level1_raw <- mean_loading_by_group(L_reference, meta_reference$annotation_level1)$matrix

# ============================================================
# GP selection: the 31 GPs active across tissues
# ============================================================
# The same filter as experiments/healthy_nonthymus_mean_loading_heatmaps/: a GP
# is kept when at least one tissue's raw mean loading reaches the cutoff. The
# expected set is spelled out so that a change in the loadings, the cell
# filter or the cutoff fails here instead of quietly redrawing a different
# figure from the one the manuscript's "31 programs" sentence describes.
raw_mean_cutoff <- 0.1
tissue_active <- rowSums(organ_raw >= raw_mean_cutoff) > 0L
selected_gps <- rownames(organ_raw)[tissue_active]

expected_gps <- paste0("GP", c(
  1, 3, 4, 6, 8, 9, 11, 22, 23, 25, 26, 29, 30, 32, 35, 41, 43, 49, 51, 58,
  62, 63, 72, 80, 93, 100, 166, 170, 171, 174, 177
))
stopifnot(
  length(expected_gps) == 31L,
  identical(selected_gps, expected_gps)
)

organ_selected_raw <- organ_raw[tissue_active, , drop = FALSE]
level1_selected_raw <- level1_raw[tissue_active, , drop = FALSE]
organ_centered <- center_by_gp_mean(organ_selected_raw)
level1_centered <- center_by_gp_mean(level1_selected_raw)

# ============================================================
# Ordering: (a) sets the row order, (b) reuses it
# ============================================================
# Recomputing the dominant-group order on the 31-row submatrix (rather than
# inheriting the retired 200-GP order) is what makes the tissue columns reflect
# the GPs actually on display. Panel (b) then keeps (a)'s row order, so a GP sits
# on the same line in both halves and can be read across; its columns are
# Figure 1's lineage order rather than a dominant-group order.
organ_order <- dominant_group_order(organ_selected_raw)

level1_groups <- colnames(level1_centered)
if (!setequal(level1_groups, level1_order)) {
  stop("The observed lineages do not match the Figure 1 level1 order.")
}
level1_column_order <- order(match(level1_groups, level1_order))
level1_row_order <- match(
  rownames(organ_centered)[organ_order$row_order],
  rownames(level1_centered)
)

stopifnot(
  nrow(organ_centered) == 31L,
  nrow(level1_centered) == 31L,
  ncol(organ_centered) == 18L,
  ncol(level1_centered) == 8L,
  identical(rownames(organ_centered), rownames(level1_centered)),
  max(abs(rowMeans(organ_centered))) < 1e-12,
  max(abs(rowMeans(level1_centered))) < 1e-12
)

organ_palette <- palette_for_groups(
  colnames(organ_centered),
  ZemmourLib::immgent_colors$organ_simplified,
  "organ_simplified"
)
level1_palette <- palette_for_groups(
  colnames(level1_centered),
  ZemmourLib::immgent_colors$level1,
  "annotation_level1"
)

# ============================================================
# s5a: the 31 tissue-active GPs across the 18 tissues
# ============================================================
render_centered_heatmap(
  organ_centered,
  organ_palette,
  "tissue (organ_simplified)",
  file.path(figure_path, "s5a.pdf"),
  organ_order$row_order,
  organ_order$column_order,
  organ_color_limit,
  "31 tissue-active GPs; dominant-group blocks",
  palette = heatmap_palettes$blue_red
)

# ============================================================
# s5b: the same 31 GPs across the 8 lineages
# ============================================================
render_centered_heatmap(
  level1_centered,
  level1_palette,
  "lineage (annotation_level1)",
  file.path(figure_path, "s5b.pdf"),
  level1_row_order,
  level1_column_order,
  level1_color_limit,
  "the same 31 GPs, row order from (a)",
  palette = heatmap_palettes$purple_green
)

# ============================================================
# What the panels drew
# ============================================================
record_dir <- "output/FigureS5"   # build intermediate (not a manuscript panel)
dir.create(record_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  data.frame(
    GP = rownames(organ_centered),
    max_raw_tissue_mean = as.numeric(apply(organ_selected_raw, 1L, max)),
    dominant_tissue = colnames(organ_selected_raw)[max.col(organ_selected_raw, ties.method = "first")],
    dominant_lineage = colnames(level1_selected_raw)[max.col(level1_selected_raw, ties.method = "first")],
    row_in_panel = match(rownames(organ_centered), rownames(organ_centered)[organ_order$row_order])
  ),
  file.path(record_dir, "s5_selected_gps.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(
    panel = c("s5a", "s5b"),
    grouping = c("organ_simplified", "annotation_level1"),
    view = "row-centered mean loading of the 31 tissue-active GPs",
    gp_count = c(nrow(organ_centered), nrow(level1_centered)),
    group_count = c(ncol(organ_centered), ncol(level1_centered)),
    selection = paste0("raw mean loading >= ", raw_mean_cutoff, " in >= 1 tissue"),
    centered_definition = "group mean minus mean across groups for each GP",
    palette = c("blue-white-red", "purple-white-green"),
    color_min = c(-organ_color_limit, -level1_color_limit),
    color_mid = 0,
    color_max = c(organ_color_limit, level1_color_limit),
    observed_min = c(min(organ_centered), min(level1_centered)),
    observed_max = c(max(organ_centered), max(level1_centered)),
    frac_saturated = c(mean(abs(organ_centered) > organ_color_limit),
                       mean(abs(level1_centered) > level1_color_limit))
  ),
  file.path(record_dir, "S5_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

# A panel that fails to write leaves the previous PDF in place and the script
# still exits 0, so check the files rather than the exit code.
expected_panels <- file.path(figure_path, c("s5a.pdf", "s5b.pdf"))
stale <- expected_panels[!file.exists(expected_panels) |
                           file.mtime(expected_panels) < run_started_at |
                           file.size(expected_panels) == 0]
if (length(stale) > 0L) {
  stop("These panels were not written by this run: ", paste(stale, collapse = ", "))
}

message("Wrote Figure S5 to ", normalizePath(figure_path))
