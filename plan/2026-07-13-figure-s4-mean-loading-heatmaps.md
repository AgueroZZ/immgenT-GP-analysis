# Figure S4 Mean-Loading Heatmaps Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Figure S4 with a tissue panel and a level2 panel, each offering raw and within-GP-normalized mean-loading heatmaps as selectable workflowr tabs.

**Architecture:** Create a standalone Figure S4 renderer that uses the healthy non-thymocyte reference population, retains raw mean loading >= 0.1, and calculates one dominant-group order per grouping. Each grouping's raw and normalized PDFs use the same retained rows, columns, palette, and order; normalized values are calculated before filtering. A new workflowr page presents the two PDF-derived image alternatives within tabsets, while the index links the page.

**Tech Stack:** R, ComplexHeatmap, circlize, ZemmourLib, workflowr, R Markdown, Quick Look PNG rendering.

---

### Task 1: Create the standalone Figure S4 renderer

**Files:**
- Create: `script/FigureS4.R`
- Create: `figures/generated/Figure S4/S4a_raw_mean_loading.pdf`
- Create: `figures/generated/Figure S4/S4a_normalized_mean_loading.pdf`
- Create: `figures/generated/Figure S4/S4b_raw_mean_loading.pdf`
- Create: `figures/generated/Figure S4/S4b_normalized_mean_loading.pdf`
- Create: `figures/generated/Figure S4/S4_filter_summary.csv`

- [x] Load `code/R/setup_data.R`, restrict to `condition_broad == "healthy"` and `annotation_level1 != "thymocyte"`, and calculate raw group means for `organ_simplified` and `annotation_level2`.
- [x] For each grouping, retain a GP row only when a raw group mean is >= 0.1, then retain group columns with at least one retained GP >= 0.1; calculate normalized values on the complete grouping matrix before applying the same subset.
- [x] Assign each retained GP to its largest raw group mean; order groups by dominant-GP count and rows by dominant-group block, dominance gap, dominant mean, and GP number.
- [x] Render raw and normalized PDFs from the same retained matrix shape, palette, and explicit order; write retained dimensions and the cutoff to `S4_filter_summary.csv`.

### Task 2: Add the Figure S4 workflowr source and static image assets

**Files:**
- Create: `analysis/FigureS4.Rmd`
- Create: `analysis/assets/FigureS4/S4a_raw_mean_loading.png`
- Create: `analysis/assets/FigureS4/S4a_normalized_mean_loading.png`
- Create: `analysis/assets/FigureS4/S4b_raw_mean_loading.png`
- Create: `analysis/assets/FigureS4/S4b_normalized_mean_loading.png`
- Modify: `analysis/index.Rmd`
- Create: `docs/assets/FigureS4/S4a_raw_mean_loading.png`
- Create: `docs/assets/FigureS4/S4a_normalized_mean_loading.png`
- Create: `docs/assets/FigureS4/S4b_raw_mean_loading.png`
- Create: `docs/assets/FigureS4/S4b_normalized_mean_loading.png`

- [x] Create a Figure S4 page that links `script/FigureS4.R`, shows the code without re-executing it, and uses an R Markdown `.tabset` under panel (a) tissue and panel (b) level2.
- [x] In each tabset, label the alternatives `Raw mean loading` and `Within-GP normalized mean loading`, display the corresponding pre-rendered PNG, and state that both use the same raw-filtered GP/group set and dominant-group order.
- [x] Add a Figure S4 link to the Supplementary figures section of `analysis/index.Rmd`.
- [x] Render each generated PDF to a PNG, place matching copies under `analysis/assets/FigureS4/` and `docs/assets/FigureS4/`, and verify each source/docs asset pair is byte-identical.

### Task 3: Build and validate the formal output

**Files:**
- Create: `docs/FigureS4.html`
- Modify: `docs/index.html`
- Modify: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`

- [x] Run `Rscript script/FigureS4.R` from the repository root; expect four nonempty one-page PDFs and a two-row filter summary.
- [x] Run `Rscript -e 'workflowr::wflow_build(c("analysis/FigureS4.Rmd", "analysis/index.Rmd"))'`; expect FigureS4.html, its tab markup, and the index link in `docs/`.
- [x] Validate that raw/normalized pairs have the same retained dimensions and order, every retained GP has maximum raw mean >= 0.1, and the four PNG/PDF pairs are nonempty.
- [x] Render and visually inspect the raw and normalized alternatives for tissue and level2; check title, legend, color strip, labels, and dominant-group structure.
- [x] Record the new formal figure, tabbed alternatives, rendering/build commands, and validation result in the update log without modifying `figures/final-selected/`.
