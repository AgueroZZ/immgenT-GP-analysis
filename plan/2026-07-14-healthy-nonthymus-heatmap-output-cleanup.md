# Healthy Non-thymocyte Heatmap Output Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Retain only the current filtered raw, normalized, and centered heatmaps plus the all-200-GP centered internal views in the healthy non-thymocyte experiment directory.

**Architecture:** Simplify the renderer so it produces only the current selected candidate family: raw/normalized views filtered at raw mean >= 0.1, centered views independently filtered at positive centered mean >= 0.1, and unfiltered centered internal views. Preserve only CSVs required to reproduce or interpret those outputs; remove hierarchical, triangular, and obsolete full raw/normalized artifacts so future renderer runs cannot recreate clutter.

**Tech Stack:** R, ComplexHeatmap, circlize, ZemmourLib.

---

### Task 1: Restrict renderer outputs to the current heatmap family

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`

- [x] Remove `triangular_first_order()` and `write_triangular_order_csv()` and remove all hierarchical and triangular rendering calls.
- [x] Retain only the six filtered final PDFs: raw, normalized, and centered tissue/level2 heatmaps, plus the two unfiltered all-200-GP centered internal PDFs.
- [x] Write filtered raw, filtered normalized, and filtered centered matrix CSVs; retain full centered matrix CSVs only for the all-200-GP internal figures.
- [x] Retain only group-count tables, raw/centered filter summaries, final filtered order CSVs, and full-centered order CSVs.

### Task 2: Rewrite experiment documentation and summary

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/README.md`
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/heatmap_summary.txt`

- [x] List the six filtered final PDFs and two full centered internal PDFs, with their filtering rules and order definitions.
- [x] State that GP37 is not displayed because its mammary-gland raw and centered maxima are 0.01123 and 0.01058, respectively, below the current 0.1 filter thresholds.
- [x] Remove references to hierarchical clustering, triangular candidates, and superseded output files.

### Task 3: Remove obsolete artifacts and validate the clean directory

**Files:**
- Delete: obsolete `*_triangular_*`, unfiltered raw/normalized PDFs and matrices, and obsolete order tables in `experiments/healthy_nonthymus_mean_loading_heatmaps/`
- Modify: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`

- [x] Run the renderer and compare `rg --files` against the intended keep set: renderer, README, summary, eight PDFs, group counts, filter summaries, filtered matrices, and matching order CSVs.
- [x] Check all eight PDFs are one-page files, verify filter summaries match their matrix dimensions, and run `git diff --check` on source and documentation.
- [x] Record the cleanup and the GP37 threshold diagnosis in the update log.
