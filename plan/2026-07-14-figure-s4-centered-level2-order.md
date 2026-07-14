# Figure S4 Centered Heatmaps and Level2 Order Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add centered mean-loading alternatives to Figure S4, order level2 groups by the Figure 1 level1 sequence and alphabetical level2 labels, and provide unfiltered 200-GP centered heatmaps in the internal experiment directory.

**Architecture:** Center each GP row by subtracting its mean across group means before applying the existing raw-mean filter for formal S4. Raw, normalized, and centered variants of each formal panel share the same retained matrix shape and row order. Level2 columns use the confirmed Figure 1 order `CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP`, with alphabetical level2 labels within every level1 block; the experiment additionally renders centered matrices without filtering any rows or columns.

**Tech Stack:** R, ComplexHeatmap, circlize, ZemmourLib, workflowr, R Markdown, Quick Look PNG rendering.

---

### Task 1: Extend Figure S4 calculations and rendering

**Files:**
- Modify: `script/FigureS4.R`
- Create: `figures/generated/Figure S4/S4a_centered_mean_loading.pdf`
- Create: `figures/generated/Figure S4/S4b_centered_mean_loading.pdf`

- [x] Add row-centering with `sweep(mean_matrix, 1L, rowMeans(mean_matrix), "-")` and use a symmetric blue-white-red color scale centered at zero.
- [x] Derive a unique `annotation_level2` to `annotation_level1` map from healthy non-thymocyte metadata; order level2 columns by `CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP`, then alphabetically within each block.
- [x] Recompute level2 row blocks using the fixed level2 column order, retain the existing tissue dominant-group order, and render raw, normalized, and centered formal alternatives from the same filtered row/column set per panel.
- [x] Add a level1 color strip above level2 heatmaps while retaining the canonical level2 strip.

### Task 2: Add Figure S4 centered tabs and static assets

**Files:**
- Modify: `analysis/FigureS4.Rmd`
- Create: `analysis/assets/FigureS4/S4a_centered_mean_loading.png`
- Create: `analysis/assets/FigureS4/S4b_centered_mean_loading.png`
- Create: `docs/assets/FigureS4/S4a_centered_mean_loading.png`
- Create: `docs/assets/FigureS4/S4b_centered_mean_loading.png`
- Modify: `docs/FigureS4.html`

- [x] Add a `Centered mean loading` tab to both panel tabsets, explaining that each GP's group means are centered by their row mean.
- [x] State the Figure 1 level1 block sequence and alphabetical within-block level2 order in the level2 panel description.
- [x] Render the two new PDFs to PNG, synchronize analysis/docs assets, rebuild the Figure S4 workflowr page, and verify the HTML contains three tabs in both panels.

### Task 3: Add full internal centered heatmaps

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/README.md`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_row_centered_mean_loading_matrix.csv`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/annotation_level2_row_centered_mean_loading_matrix.csv`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_row_centered_mean_loading_full_dominant_group_heatmap.pdf`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/annotation_level2_row_centered_mean_loading_full_level1_order_heatmap.pdf`

- [x] Add the same row-centering and symmetric centered color scale to the experiment renderer.
- [x] Render unfiltered all-200-GP centered tissue and level2 heatmaps; use the original full level2 columns in the Figure 1 level1/alphabetical order and retain all observed groups.
- [x] Document the internal outputs, centering definition, unfiltered dimensions, and level2 ordering.

### Task 4: Validate and record the update

**Files:**
- Modify: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`

- [x] Run `Rscript script/FigureS4.R` and `Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`; expect successful completion.
- [x] Verify every centered GP row has mean zero within numerical tolerance, formal pairs share rows/columns/orders, full internal centered matrices have 200 rows, and level2 columns follow the confirmed level1/alphabetical order.
- [x] Visually inspect centered tissue and level2 PDFs and the rebuilt tabbed Figure S4 page for color scale, annotation strips, labels, and legibility.
- [x] Record source files, output dimensions, validation results, and the explicit level1 sequence in the update log.
