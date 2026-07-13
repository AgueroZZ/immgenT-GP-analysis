# Healthy Non-thymus Tissue Dominant Dotplot Implementation Plan

> **For agentic workers:** Execute these steps in order and keep all outputs isolated under `experiments/healthy_nonthymus_mean_loading_heatmaps/`.

**Goal:** Render one tissue dotplot candidate using the dominant-group order, with normalized mean loading encoded by dot color and raw mean loading encoded by dot area.

**Architecture:** Reuse the existing healthy non-thymocyte tissue raw and normalized matrices and the fixed dominant-group row/column order. Draw one circle per GP-tissue cell on a white background: the color is the within-GP normalized mean in `[0, 1]`, while circle area uses the unnormalized raw mean on one global tissue scale. A small canonical tissue-color strip remains above the columns.

**Tech Stack:** R, ggplot2, patchwork, ZemmourLib, PDF.

---

### Task 1: Implement the dual-scale dotplot renderer

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`

- [x] **Step 1: Convert aligned raw and normalized matrices to one dot per GP-tissue pair.**

```r
dot_df <- data.frame(
  group = factor(rep(colnames(raw_matrix), each = nrow(raw_matrix)), levels = ordered_groups),
  GP = factor(rep(rownames(raw_matrix), times = ncol(raw_matrix)), levels = rev(ordered_gps)),
  raw_mean = as.vector(raw_matrix),
  normalized_mean = as.vector(normalized_matrix)
)
```

- [x] **Step 2: Map color to normalized mean and point area to raw mean.**

```r
geom_point(aes(size = raw_mean, colour = normalized_mean), shape = 16)
scale_colour_gradient(limits = c(0, 1), low = "#FEE5D9", high = "#99000D")
scale_size_area(limits = c(0, max(raw_matrix)), max_size = 3.2)
```

- [x] **Step 3: Add the canonical tissue-color strip and retain the dominant-group order.**

Expected: no tile background; a tiny high-normalized/raw-mean value remains a tiny dot rather than a prominent red tile.

### Task 2: Render and validate the tissue candidate

**Files:**
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_dominant_group_normalized_color_raw_area_dotplot.pdf`
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/README.md`
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/heatmap_summary.txt`
- Modify: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`

- [x] **Step 1: Run the renderer from the repository root.**

Run: `Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`

Expected: exit status 0 and a nonempty dotplot PDF in the experiment directory.

- [x] **Step 2: Check the matrix/order alignment and raw-area scale.**

Run: `Rscript -e 'x <- read.csv("experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_dominant_group_order.csv"); stopifnot(sum(x$dimension == "GP row") == 200L); cat("OK\\n")'`

Expected: `OK`; the rendered plot uses the 200 rows and 18 groups in that order.

- [x] **Step 3: Render the PDF thumbnail and inspect dot visibility, legends, labels, and the tissue annotation strip.**

Run: `/usr/bin/qlmanage -t -s 1800 -o /private/tmp/gp_mean_loading_dotplot experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_dominant_group_normalized_color_raw_area_dotplot.pdf`

Expected: large red dots indicate cells with both high normalized and high raw mean; high-normalized low-raw cells remain small.
