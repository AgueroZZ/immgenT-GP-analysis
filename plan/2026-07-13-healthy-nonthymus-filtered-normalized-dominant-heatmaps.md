# Filtered Normalized Dominant-Group Heatmaps Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make each final normalized dominant-group heatmap use exactly the same retained GPs, groups, and ordering as its filtered raw counterpart.

**Architecture:** Preserve the full normalized matrices and all hierarchical/triangular candidates. Subset the existing within-GP normalized matrices using the raw-filtered matrix dimnames, reuse the raw-filtered group metadata and order, and overwrite only the two final normalized dominant-group PDF paths with these matched views.

**Tech Stack:** R, ComplexHeatmap, circlize, base R matrix subsetting, Markdown.

---

### Task 1: Construct matching normalized matrices

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R:420-470`

- [x] Subset each full normalized matrix with `rownames(raw_filtered)` and `colnames(raw_filtered)` to form tissue and level2 final normalized matrices.
- [x] Assert identical dimnames between each raw-filtered matrix and its normalized counterpart, and assert every retained GP has full-matrix raw maximum >= 0.1.
- [x] Write the two filtered normalized display matrices as CSV files with `raw_mean_ge_0.1_dominant_group` in their filenames.

### Task 2: Render the matched normalized dominant-group heatmaps

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R:615-650`

- [x] Render each normalized final PDF from the matching filtered normalized matrix, retained group metadata, retained palette, and the same filtered dominant-group order used by its raw PDF.
- [x] Label both final raw and normalized heatmaps with `focus: raw mean >= 0.1 filter` and report matched final dimensions in `heatmap_summary.txt`.

### Task 3: Document and verify matched views

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/README.md:62-73`
- Modify: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`

- [x] Document that the final raw and normalized heatmaps for each grouping share filtered rows, columns, and order; full normalized matrices and non-final candidate PDFs remain available.
- [x] Run `Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R` from the repository root; expect successful completion.
- [x] Verify equal dimnames and identical row/column order tables for each final raw/normalized pair, with every retained GP maximum raw mean >= 0.1.
- [x] Render the two final normalized PDFs to PNGs and visually inspect their labels, color bars, and matched dominant-group structure.
