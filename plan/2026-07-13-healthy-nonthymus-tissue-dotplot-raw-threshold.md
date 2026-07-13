# Healthy Non-thymus Tissue Dotplot Raw-threshold Implementation Plan

> **For agentic workers:** Execute these steps in order and keep all outputs isolated under `experiments/healthy_nonthymus_mean_loading_heatmaps/`.

**Goal:** Create a filtered tissue dominant-group dotplot retaining only GPs and tissues with at least one raw mean loading of 0.1 or greater.

**Architecture:** Filter the tissue raw matrix first by GP, then by tissue using the same raw-mean threshold. Recompute dominant-group ordering exclusively on the filtered submatrix so block sizes and tissue order reflect the displayed strong signals. Retain normalized means only for dot color; raw means remain the filtering and dot-area quantity.

**Tech Stack:** R, ggplot2, patchwork, PDF.

---

### Task 1: Filter and reorder the tissue matrix

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_raw_mean_ge_0.1_filter_summary.csv`

- [x] **Step 1: Retain GPs with one or more raw group means at least 0.1.**

```r
raw_mean_cutoff <- 0.1
keep_gp <- rowSums(organ_raw >= raw_mean_cutoff) > 0L
filtered_raw <- organ_raw[keep_gp, , drop = FALSE]
filtered_normalized <- organ_normalized[keep_gp, , drop = FALSE]
```

- [x] **Step 2: Retain tissues containing at least one retained GP at the same threshold.**

```r
keep_group <- colSums(filtered_raw >= raw_mean_cutoff) > 0L
filtered_raw <- filtered_raw[, keep_group, drop = FALSE]
filtered_normalized <- filtered_normalized[, keep_group, drop = FALSE]
```

- [x] **Step 3: Recompute the dominant-group order on the filtered raw matrix.**

```r
filtered_dominant <- dominant_group_order(filtered_raw)
stopifnot(nrow(filtered_raw) == 31L, ncol(filtered_raw) == 18L)
```

### Task 2: Render and verify the filtered candidate

**Files:**
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_dominant_group_normalized_color_raw_area_raw_mean_ge_0.1_dotplot.pdf`
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/README.md`
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/heatmap_summary.txt`
- Modify: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`

- [x] **Step 1: Render the filtered dotplot with normalized color and raw dot area.**

Run: `Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`

Expected: a nonempty PDF with 31 GP rows and 18 tissue columns.

- [x] **Step 2: Verify all displayed GP rows and tissue columns meet the raw threshold condition.**

Run: `Rscript -e 'x <- read.csv("experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_raw_mean_ge_0.1_filter_summary.csv"); stopifnot(all(x$max_raw_mean >= 0.1)); cat("OK\\n")'`

Expected: `OK`.

- [x] **Step 3: Render the PDF thumbnail and inspect dot size, color, label readability, and dominant-group blocks.**

Run: `/usr/bin/qlmanage -t -s 1800 -o /private/tmp/gp_mean_loading_dotplot_filtered experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_dominant_group_normalized_color_raw_area_raw_mean_ge_0.1_dotplot.pdf`

Expected: the sparse low-amplitude rows are removed while large red dots and the dominant-group structure remain readable.
