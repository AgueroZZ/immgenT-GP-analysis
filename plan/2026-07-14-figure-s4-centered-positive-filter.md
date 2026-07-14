# Figure S4 Centered Positive-Threshold Filter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the centered Figure S4 alternatives independently retain groups and GPs using only positive row-centered mean loading >= 0.1.

**Architecture:** Calculate raw group means, calculate the full row-centered matrix by subtracting every GP's group-average mean, and apply a centered-positive filter to that matrix. The centered formal alternatives then share their own retained rows, columns, and dominant-group order, while raw and normalized alternatives retain the existing raw-positive-filtered set. The level2 centered columns retain the Figure 1 level1-first/alphabetical ordering after centered filtering.

**Tech Stack:** R, ComplexHeatmap, circlize, ZemmourLib, workflowr, Quick Look PNG rendering.

---

### Task 1: Implement and validate centered-positive filtering

**Files:**
- Modify: `script/FigureS4.R`
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`

- [x] Add a filter helper that keeps a GP if `rowSums(centered_matrix >= 0.1) > 0L`, then keeps a group if `colSums(filtered_centered >= 0.1) > 0L`.
- [x] Use the helper only for the formal centered matrices and filtered experiment centered candidates; retain the raw filtering helper for raw and normalized alternatives.
- [x] Recompute tissue and level2 centered dominant ordering from the centered-filtered matrices, with level2 columns ordered by `CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP`, then alphabetically within level1.
- [x] Add assertions that every centered-filtered row and column has at least one positive centered entry >= 0.1, every centered full row still has mean zero, and fixed orders are complete permutations.

### Task 2: Regenerate formal Figure S4 and update its explanation

**Files:**
- Modify: `analysis/FigureS4.Rmd`
- Modify: `figures/generated/Figure S4/S4a_centered_mean_loading.pdf`
- Modify: `figures/generated/Figure S4/S4b_centered_mean_loading.pdf`
- Modify: `analysis/assets/FigureS4/S4a_centered_mean_loading.png`
- Modify: `analysis/assets/FigureS4/S4b_centered_mean_loading.png`
- Modify: `docs/assets/FigureS4/S4a_centered_mean_loading.png`
- Modify: `docs/assets/FigureS4/S4b_centered_mean_loading.png`
- Modify: `docs/FigureS4.html`

- [x] State that only centered tabs use an independent positive centered threshold of 0.1 after full row-centering, and report their retained dimensions without changing raw/normalized descriptions.
- [x] Render the two centered PDFs with the independently filtered matrices, convert them to PNGs, synchronize analysis/docs assets, and rebuild `docs/FigureS4.html`.
- [x] Confirm both panel tabsets still expose raw, normalized, and centered views, and confirm the centered level2 column order follows the Figure 1 sequence among retained groups.

### Task 3: Synchronize experiment documentation and change record

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/README.md`
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/heatmap_summary.txt`
- Modify: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`

- [x] Write the centered-positive filter definition and separate retained dimensions in the experiment README and summary.
- [x] Record the raw/normalized versus centered filter distinction, exact threshold, rendered outputs, and validation results in the update log.

### Task 4: Execute and inspect

**Files:**
- Test: `script/FigureS4.R`
- Test: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`

- [x] Run `Rscript script/FigureS4.R` and `Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`; expect both to exit successfully.
- [x] Inspect centered matrix dimensions, confirm the positive threshold for every retained row/column, validate full centered-row means are within `1e-12` of zero, and inspect the centered tissue and level2 PDFs.
- [x] Run `git diff --check` on all changed source, documentation, plan, and log files.
