# Figure S4 Full Centered Final Implementation Plan

Status: completed on 2026-07-15.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Replace the exploratory Figure S4 alternatives with two final, unfiltered row-centered heatmaps using a fixed symmetric color range of `[-0.2, 0.2]`.

**Architecture:** Keep `script/FigureS4.R` as the single formal renderer. Compute full tissue and level2 mean-loading matrices, center every GP row, order rows by the dominant raw-mean group, and render all 200 GPs without row or column filtering. Keep the experiment directory unchanged and remove obsolete alternatives only from the formal Figure S4 outputs and web assets.

**Tech Stack:** R, ComplexHeatmap, circlize, workflowr/rmarkdown, macOS Quick Look.

---

### Task 1: Simplify the formal renderer

**Files:**
- Modify: `script/FigureS4.R`
- Replace: `figures/generated/Figure S4/S4_summary.csv`
- Remove: `figures/generated/Figure S4/S4_filter_summary.csv`
- Remove: formal raw and normalized PDFs under `figures/generated/Figure S4/`

- [x] **Step 1: Remove formal filtering and alternative-scale code**

Keep only the raw matrices needed for ordering and the centered matrices needed for display:

```r
organ_centered <- center_by_gp_mean(organ_raw)
level2_centered <- center_by_gp_mean(level2_raw)
stopifnot(nrow(organ_centered) == 200L, nrow(level2_centered) == 200L)
```

- [x] **Step 2: Preserve the established full-matrix ordering**

```r
organ_order <- dominant_group_order(organ_raw)
level2_order <- dominant_group_order(
  level2_raw,
  level2_column_order(colnames(level2_raw), level2_group_level1, level1_order)
)
```

- [x] **Step 3: Fix the centered color range**

```r
centered_color_limit <- 0.2
color_fun <- circlize::colorRamp2(
  c(-centered_color_limit, 0, centered_color_limit),
  c("#2166AC", "#FFFFFF", "#B2182B")
)
```

Values outside the range must saturate at the endpoint colors.

- [x] **Step 4: Render only the two final PDFs**

```r
render_heatmap(organ_centered, ..., "S4a_centered_mean_loading.pdf")
render_heatmap(level2_centered, ..., "S4b_centered_mean_loading.pdf")
```

- [x] **Step 5: Write a non-filter summary and validate dimensions**

Expected summary rows:

```text
S4a,organ_simplified,200,18,-0.2,0.2
S4b,annotation_level2,200,107,-0.2,0.2
```

Run: `Rscript script/FigureS4.R`

Expected: two one-page PDFs and a summary confirming 200 rows in both panels.

### Task 2: Finalize the analysis page

**Files:**
- Modify: `analysis/FigureS4.Rmd`
- Modify: `docs/FigureS4.html`
- Keep: `analysis/assets/FigureS4/S4a_centered_mean_loading.png`
- Keep: `analysis/assets/FigureS4/S4b_centered_mean_loading.png`
- Keep: matching files under `docs/assets/FigureS4/`
- Remove: raw and normalized Figure S4 PNG assets from both asset directories

- [x] **Step 1: Remove both tabsets and all alternative-view prose**

The page must have one image and one caption per panel, with no raw, normalized, filtered, or undecided wording.

- [x] **Step 2: State the final matrix and color definitions**

```markdown
Each panel shows all 200 GPs and all observed groups. Each GP row is centered by subtracting its mean across groups. The shared color scale is fixed at -0.2 to 0.2; values outside this range are saturated.
```

- [x] **Step 3: Rebuild and synchronize the page**

Run:

```bash
Rscript -e 'rmarkdown::render("analysis/FigureS4.Rmd", output_file="FigureS4.html", quiet=TRUE)'
cp analysis/FigureS4.html docs/FigureS4.html
```

Expected: the rendered HTML references only the two centered PNGs.

### Task 3: Visual QA and documentation

**Files:**
- Modify: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`
- Modify: `plan/2026-07-15-figure-s4-full-centered-final.md`

- [x] **Step 1: Render both PDFs to PNG and inspect them**

Check that all 200 row labels and all group labels are present, vertical column labels do not overlap, the level1/level2 strips remain intact, and values near `0.2` use the deepest endpoint color.

- [x] **Step 2: Verify formal-output cleanup**

Run:

```bash
rg --files "figures/generated/Figure S4" analysis/assets/FigureS4 docs/assets/FigureS4
```

Expected: two centered PDFs, two PNGs per asset location, and `S4_summary.csv`; no raw, normalized, or filter-summary files.

- [x] **Step 3: Record the final collaborator decision**

Append the final full-centered selection, the 200-by-18 and 200-by-107 dimensions, the fixed `[-0.2, 0.2]` scale, render checks, and cleanup result to the project log.

## Validation result

- S4a is 200 GPs by 18 tissues; S4b is 200 GPs by 107 level2 clusters.
- Both PDFs are one page and use a fixed `[-0.2, 0.2]` centered scale.
- Values beyond either endpoint were programmatically confirmed to saturate at the endpoint color.
- Quick Look inspection confirmed deeper colors, complete row labels, vertical group labels, and intact annotations.
- Formal outputs contain only two centered PDFs, one summary CSV, and the two matching PNGs in each web asset directory.
- The experiment directory has no new diff from this update.
