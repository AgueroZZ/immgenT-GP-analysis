# Filtered Raw Dominant-Group Heatmaps Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the final tissue and level2 raw dominant-group heatmaps more readable by retaining only rows and columns with at least one raw mean loading of 0.1 or greater.

**Architecture:** Keep full matrices, all normalized heatmaps, and other candidate heatmaps unchanged. Create filtered copies of each raw matrix, subset their matching group-count tables and palettes, recompute the dominant-group order on each filtered matrix, and use those filtered matrices only for the two final raw dominant-group PDF paths.

**Tech Stack:** R, ComplexHeatmap, circlize, base R matrix operations, Markdown.

---

### Task 1: Define and apply the raw-mean filter

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R:145-165,381-430`

- [x] Add `filter_raw_mean_matrix(raw_matrix, raw_mean_cutoff)` that retains a row when `rowSums(raw_matrix >= raw_mean_cutoff) > 0L`, then retains a column when `colSums(filtered_raw >= raw_mean_cutoff) > 0L`, and stops if fewer than two groups remain.
- [x] Set `raw_dominant_cutoff <- 0.1`; create filtered tissue and level2 raw matrices, matching group-count tables, matching palettes, and new dominant-group orders.
- [x] Write two filter-summary CSV files recording retained GP count, retained group count, and the cutoff.

### Task 2: Render only the final raw dominant-group PDFs from filtered matrices

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R:531-568,583-598`

- [x] Pass each filtered raw matrix, subset group-count table, subset palette, and recomputed dominant order to the existing raw dominant-group PDF path.
- [x] Use the order description `dominant-group blocks (within block: dominance gap; raw mean >= 0.1 filter)` and add retained dimensions to `heatmap_summary.txt`.
- [x] Leave the normalized dominant heatmaps and all hierarchical/triangular candidates on their full matrices and existing orders.

### Task 3: Document and validate the filtered outputs

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/README.md:55-64`
- Modify: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`

- [x] State that only the two final raw dominant-group heatmaps are filtered at raw mean loading >= 0.1 and that their orders are recomputed on the retained matrices.
- [x] Run `Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R` from the repository root; expect successful completion.
- [x] Run an R validation that verifies every retained row and column contains a value >= 0.1, the two raw dominant-group PDFs are nonempty, and the normalized dominant-group PDFs retain 200 rows.
- [x] Render the filtered raw tissue and level2 PDFs to PNGs and visually inspect their labels, legends, and dominant-group structure.
