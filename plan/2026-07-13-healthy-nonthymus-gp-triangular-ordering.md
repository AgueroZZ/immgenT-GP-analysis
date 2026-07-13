# Healthy Non-thymus GP Triangular-ordering Implementation Plan

> **For agentic workers:** Execute these steps in order and keep all outputs isolated under `experiments/healthy_nonthymus_mean_loading_heatmaps/`.

**Goal:** Add deterministic, triangular-first ordering candidates for the healthy non-thymus GP mean-loading heatmaps without replacing the existing hierarchical-clustering views.

**Architecture:** Derive one order per grouping from the within-GP normalized mean-loading matrix using the Figure 6b support-ordering method. A cell is visible when its relative mean loading is at least 0.50; columns are ordered from most to fewest visible GPs, then GP rows are ordered by the rightmost visible column, visible-support count, rarity-weighted support, and numeric GP ID. Apply the same fixed order to raw and normalized matrices within each grouping so their color-scale comparison is direct.

**Tech Stack:** R, ComplexHeatmap, circlize, PDF, ZemmourLib.

---

### Task 1: Implement deterministic support-based ordering

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`
- Test: `experiments/healthy_nonthymus_mean_loading_heatmaps/*_triangular_order.csv`

- [x] **Step 1: Define support from the normalized matrix.**

```r
support_cutoff <- 0.50
visible_mask <- normalized_matrix >= support_cutoff
stopifnot(all(rowSums(visible_mask) > 0L), all(colSums(visible_mask) > 0L))
```

- [x] **Step 2: Order columns by decreasing support and rows by the right boundary.**

```r
column_order <- order(-colSums(visible_mask), colnames(visible_mask))
ordered_mask <- visible_mask[, column_order, drop = FALSE]
rightmost <- apply(ordered_mask, 1L, function(x) max(which(x)))
rarity <- as.numeric(ordered_mask %*% seq_len(ncol(ordered_mask))^2)
gp_id <- as.integer(sub("^GP", "", rownames(ordered_mask)))
row_order <- order(-rightmost, -rowSums(visible_mask), -rarity, gp_id)
```

- [x] **Step 3: Write an order table with support diagnostics.**

```r
stopifnot(sum(diff(rightmost[row_order]) > 0L) == 0L)
stopifnot(sum(diff(colSums(visible_mask)[column_order]) > 0L) == 0L)
```

### Task 2: Render the triangular candidates

**Files:**
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_raw_mean_loading_triangular_heatmap.pdf`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_row_normalized_mean_loading_triangular_heatmap.pdf`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/annotation_level2_raw_mean_loading_triangular_heatmap.pdf`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/annotation_level2_row_normalized_mean_loading_triangular_heatmap.pdf`

- [x] **Step 1: Extend the renderer to accept explicit row and column order vectors.**

```r
Heatmap(
  display_matrix,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = triangular$row_order,
  column_order = triangular$column_order
)
```

- [x] **Step 2: Apply each grouping's fixed order to both its raw and normalized matrices.**

```r
render_heatmap(organ_raw, row_order = organ_triangular$row_order,
               column_order = organ_triangular$column_order)
render_heatmap(organ_normalized, row_order = organ_triangular$row_order,
               column_order = organ_triangular$column_order)
```

- [x] **Step 3: Preserve the existing four hierarchical-clustering PDFs.**

Expected: the experiment directory contains eight nonempty PDFs: four hierarchical and four triangular-first candidates.

### Task 3: Verify and document the candidates

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/README.md`
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/heatmap_summary.txt`
- Modify: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`

- [x] **Step 1: Run the renderer and verify all triangular PDFs and order tables.**

Run: `Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`

Expected: exit status 0, four triangular PDFs, and two `<grouping>_triangular_order.csv` files.

- [x] **Step 2: Verify ordering invariants and order identity across raw and normalized views.**

Run: `Rscript -e 'o <- read.csv("experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_triangular_order.csv"); stopifnot(all(diff(o$rightmost_visible_column[o$dimension == "GP row"]) <= 0)); cat("OK\\n")'`

Expected: `OK`; repeat against `annotation_level2_triangular_order.csv`.

- [x] **Step 3: Render the updated tissue and level-2 normalized PDFs to PNG and inspect the staircase boundary, labels, and legend.**

Run: `/usr/bin/qlmanage -t -s 1800 -o /private/tmp/gp_mean_loading_heatmaps_triangular experiments/healthy_nonthymus_mean_loading_heatmaps/annotation_level2_row_normalized_mean_loading_triangular_heatmap.pdf`

Expected: one unclipped thumbnail with a monotone right boundary in the support mask and readable labels at full PDF size.
