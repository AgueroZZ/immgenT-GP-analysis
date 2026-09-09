# Figure S5. The tissue-associated GPs, across tissues and across lineages.
#
# Panels produced (see analysis/FigureS5.Rmd for the caption text):
#   s5a  Row-centered mean GP activity across the 18 tissues, restricted to the
#        32 tissue-associated GPs. Blue-white-red scale.
#   s5b  The same 32 GPs in the same row order, but across the 8 T cell
#        lineages. Green-white-purple scale (purple positive), so the two halves
#        of the figure cannot be mistaken for each other.
#
# --- internal ---
# New Extended Data Figure 5 on 2026-09-09. Panel s5a is the retired Figure S4's
# panel s4a with its rows filtered to the selected GPs and the dominant-group
# order recomputed on that submatrix; s5b is new. The retired figure's other
# panel (the cluster heatmap) became Figure S3. The row set started as rule A's
# 31 alone; Ziang widened it to the union with rule B the same day, so that
# every GP the manuscript calls tissue-associated appears here.
# --- end internal ---
#
# The rows are the union of the two tissue-associated GP sets the manuscript
# names -- see "GP selection" below. The panels then show, per GP, how that
# activity is distributed, by subtracting the GP's mean across the groups
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
#   organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds
#                                            [code/pipeline/02_compute_auc.R]

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
# GP selection: the union of the two tissue-associated GP sets
# ============================================================
# The manuscript names two, by different criteria, and this figure shows both:
#
#   A  "31 programs active across tissues" -- raw (uncentered) mean loading
#      reaches 0.1 in at least one of the 18 tissues. The filter from
#      experiments/healthy_nonthymus_mean_loading_heatmaps/. It is a threshold
#      on *average* activity, so a program carried by a few cells fails it.
#   B  "Seven GPs were strongly tissue-associated (AUC > 0.9)" -- Figure 5's
#      rule: one-vs-rest tissue AUC above 0.9 under Figure 5's positivity mask
#      (mean loading in the tissue above the GP's overall mean), with tissues
#      under 100 cells dropped.
#
# Six GPs satisfy both, so the union is 32. GP37 is the only member of B alone:
# its mammary-gland AUC is 0.9918, but it is active in 2.5% of mammary gland
# cells, which dilutes its mean loading there to 0.011 -- an order of magnitude
# under rule A's cutoff. Its row is therefore near-white in (a) by construction.
#
# Rule B is re-derived here rather than retyped, and checked against the literal
# `gps_of_interest` in script/Figure5.R, so the two figures cannot disagree
# about which GPs are the strongly tissue-associated ones. The union is checked
# against a literal list as well, so a change in the loadings, the AUCs, the
# cell filter or either cutoff fails here instead of quietly redrawing a
# different figure.
raw_mean_cutoff <- 0.1
auc_cutoff <- 0.9
min_tissue_cells <- 100L

# --- A: mean activity ---
tissue_active <- rowSums(organ_raw >= raw_mean_cutoff) > 0L

# --- B: tissue discrimination, on Figure 5's terms ---
organ_auc <- readRDS(paste0(
  "data/", "organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds"
))$auc
colnames(organ_auc) <- paste0("GP", seq_len(ncol(organ_auc)))
tissue_counts <- table(meta_reference$organ_simplified)
large_tissues <- setdiff(
  rownames(organ_auc), names(tissue_counts[tissue_counts < min_tissue_cells])
)
organ_auc <- organ_auc[large_tissues, , drop = FALSE]

overall_mean <- colMeans(L_reference, na.rm = TRUE)
tissue_mean <- t(vapply(
  large_tissues,
  function(tissue) {
    colMeans(L_reference[meta_reference$organ_simplified == tissue, , drop = FALSE], na.rm = TRUE)
  },
  numeric(ncol(L_reference))
))
positive_mask <- sweep(tissue_mean, 2L, overall_mean, "-") > 0
tissue_specific <- colSums((organ_auc > auc_cutoff) & positive_mask) > 0L
tissue_specific <- rownames(organ_raw) %in% colnames(organ_auc)[tissue_specific]

# Figure 5 highlights these seven throughout; rule B must reproduce exactly them.
figure5_gps_of_interest <- local({
  lines <- readLines("script/Figure5.R", warn = FALSE)
  hit <- grep("^gps_of_interest <- ", lines, value = TRUE)
  if (length(hit) != 1L) {
    stop("gps_of_interest is not assigned exactly once in script/Figure5.R")
  }
  eval(parse(text = sub("^gps_of_interest <- ", "", hit)))
})
stopifnot(setequal(rownames(organ_raw)[tissue_specific], figure5_gps_of_interest))

selected <- tissue_active | tissue_specific
selected_gps <- rownames(organ_raw)[selected]

expected_gps <- paste0("GP", c(
  1, 3, 4, 6, 8, 9, 11, 22, 23, 25, 26, 29, 30, 32, 35, 37, 41, 43, 49, 51, 58,
  62, 63, 72, 80, 93, 100, 166, 170, 171, 174, 177
))
stopifnot(
  sum(tissue_active) == 31L,
  sum(tissue_specific) == 7L,
  length(expected_gps) == 32L,
  identical(selected_gps, expected_gps)
)

organ_selected_raw <- organ_raw[selected, , drop = FALSE]
level1_selected_raw <- level1_raw[selected, , drop = FALSE]
organ_centered <- center_by_gp_mean(organ_selected_raw)
level1_centered <- center_by_gp_mean(level1_selected_raw)

# ============================================================
# Ordering: (a) sets the row order, (b) reuses it
# ============================================================
# Recomputing the dominant-group order on the 32-row submatrix (rather than
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
  nrow(organ_centered) == 32L,
  nrow(level1_centered) == 32L,
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
  "32 tissue-associated GPs; dominant-group blocks",
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
  "the same 32 GPs, row order from (a)",
  palette = heatmap_palettes$green_purple
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
    max_masked_tissue_auc = as.numeric(apply(
      ifelse(positive_mask[, selected, drop = FALSE], organ_auc[, selected, drop = FALSE], NA_real_),
      2L, max, na.rm = TRUE
    )),
    rule = ifelse(
      tissue_active[selected] & tissue_specific[selected], "A+B",
      ifelse(tissue_active[selected], "A", "B")
    ),
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
    view = "row-centered mean loading of the 32 tissue-associated GPs",
    gp_count = c(nrow(organ_centered), nrow(level1_centered)),
    group_count = c(ncol(organ_centered), ncol(level1_centered)),
    selection = paste0("raw mean loading >= ", raw_mean_cutoff,
                       " in >= 1 tissue OR masked tissue AUC > ", auc_cutoff),
    centered_definition = "group mean minus mean across groups for each GP",
    palette = c("blue-white-red", "green-white-purple"),
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
