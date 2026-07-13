# Healthy Non-thymus GP Dominant-group Ordering Implementation Plan

> **For agentic workers:** Execute these steps in order and keep all outputs isolated under `experiments/healthy_nonthymus_mean_loading_heatmaps/`.

**Goal:** Add a Figure 4c-style dominant-group ordering for the GP-by-tissue and GP-by-level2 mean-loading heatmaps, then compare it with the existing support-based triangular order in tissue.

**Architecture:** For each GP, identify the tissue or level2 cluster with its largest raw mean loading. Order groups by how many GPs select them as dominant, and form GP blocks by dominant group. Within a block, order GPs by decreasing dominance gap, defined as the largest minus the second-largest group mean; this emphasizes group specificity without thresholding continuous values. Reuse each grouping's order for its raw and normalized heatmaps.

**Tech Stack:** R, ComplexHeatmap, circlize, ZemmourLib, PDF.

---

### Task 1: Calculate dominant-group blocks

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_dominant_group_order.csv`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/annotation_level2_dominant_group_order.csv`

- [x] **Step 1: Assign every GP to the group with its maximum raw mean loading.**

```r
dominant_index <- max.col(raw_mean_matrix, ties.method = "first")
dominant_mean <- raw_mean_matrix[cbind(seq_len(nrow(raw_mean_matrix)), dominant_index)]
second_mean <- apply(raw_mean_matrix, 1L, function(x) sort(x, decreasing = TRUE)[2])
dominance_gap <- dominant_mean - second_mean
```

- [x] **Step 2: Sort groups and GP blocks deterministically.**

```r
dominant_gp_count <- tabulate(dominant_index, nbins = ncol(raw_mean_matrix))
column_order <- order(-dominant_gp_count, -colMeans(raw_mean_matrix), colnames(raw_mean_matrix))
dominant_position <- match(dominant_index, column_order)
row_order <- order(dominant_position, -dominance_gap, -dominant_mean, gp_number)
```

- [x] **Step 3: Record the dominant group, two leading means, gap, and group-level counts.**

Expected: GP dominant-group blocks are contiguous and their block positions never decrease down the ordered rows.

### Task 2: Render dominant-group candidates

**Files:**
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_raw_mean_loading_dominant_group_heatmap.pdf`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_row_normalized_mean_loading_dominant_group_heatmap.pdf`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/annotation_level2_raw_mean_loading_dominant_group_heatmap.pdf`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/annotation_level2_row_normalized_mean_loading_dominant_group_heatmap.pdf`

- [x] **Step 1: Apply the tissue dominant-group order to both tissue color-scale variants.**

```r
render_heatmap(organ_raw, row_order = organ_dominant$row_order,
               column_order = organ_dominant$column_order)
render_heatmap(organ_normalized, row_order = organ_dominant$row_order,
               column_order = organ_dominant$column_order)
```

- [x] **Step 2: Apply the independently derived level2 order to both level2 variants.**

Expected: the same grouping has identical raw and normalized row/column order, while tissue and level2 retain independent biologically derived orders.

- [x] **Step 3: Preserve all eight existing hierarchical and support-triangular PDFs.**

Expected: four dominant-group PDFs are added; no existing PDF filename is replaced.

### Task 3: Compare tissue ordering and validate

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/README.md`
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/heatmap_summary.txt`
- Modify: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`

- [x] **Step 1: Run the renderer and require twelve nonempty heatmap PDFs.**

Run: `Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`

Expected: exit status 0 and twelve nonempty PDFs: four hierarchical, four support-triangular, and four dominant-group.

- [x] **Step 2: Check both order tables.**

Run: `Rscript -e 'x <- read.csv("experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_dominant_group_order.csv"); r <- x[x$dimension == "GP row", ]; stopifnot(all(diff(r$dominant_group_position) >= 0)); cat("OK\\n")'`

Expected: `OK`; the same check passes for level2.

- [x] **Step 3: Render and visually compare tissue normalized support-triangular and dominant-group PDFs.**

Expected: support-triangular ordering produces a monotone right boundary based on a 0.50 mask; dominant-group ordering produces contiguous tissue blocks ordered by the strongest mean-loading tissue, with within-block specificity controlled by the continuous dominance gap.
