# Remove Tissue Dotplots Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Retain the chosen dominant-group heatmaps and remove all tissue dotplot generation and artifacts from the isolated experiment.

**Architecture:** Keep the existing matrix construction and the hierarchical, triangular, and dominant-group heatmap renderers unchanged. Delete the dotplot-specific functions, calls, metadata, README text, and generated dotplot/filter files so rerunning the experiment produces only heatmap outputs.

**Tech Stack:** R, ComplexHeatmap, base R file operations, Markdown.

---

### Task 1: Remove dotplot rendering from the experiment script

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R:286-414,687-746,788-793`

- [x] Delete `render_dominant_dotplot()` and `filter_raw_mean_matrix()` because no remaining output uses them.
- [x] Delete the three tissue dotplot render calls, the raw-mean filter computation, and the filter-summary CSV write.
- [x] Delete the dotplot-only lines from `heatmap_summary.txt`, preserving all heatmap and dominant-group summary lines.

### Task 2: Remove dotplot documentation and generated artifacts

**Files:**
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/README.md:64-82`
- Delete: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_dominant_group_normalized_color_raw_area_dotplot.pdf`
- Delete: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_dominant_group_normalized_color_raw_area_raw_mean_ge_0.1_dotplot.pdf`
- Delete: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_dominant_group_raw_color_normalized_area_raw_mean_ge_0.1_dotplot.pdf`
- Delete: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_raw_mean_ge_0.1_filter_summary.csv`

- [x] Remove the dotplot descriptions from the README while retaining the selected dominant-group heatmap documentation.
- [x] Delete the generated dotplot PDFs and their filter-only CSV.

### Task 3: Render and verify the heatmap-only experiment

**Files:**
- Verify: `experiments/healthy_nonthymus_mean_loading_heatmaps/`
- Modify: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`

- [x] Run `Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R` from the repository root; expect successful completion.
- [x] Verify that all four dominant-group heatmaps exist and are nonempty: tissue raw/normalized and level2 raw/normalized.
- [x] Verify `rg -n "dotplot" experiments/healthy_nonthymus_mean_loading_heatmaps` returns no matches, then record the cleanup and validation in the update log.
