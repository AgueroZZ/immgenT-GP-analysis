# Healthy Non-thymus GP Mean-loading Heatmaps Implementation Plan

> **For agentic workers:** Execute these steps in order and keep all analysis outputs isolated under `experiments/healthy_nonthymus_mean_loading_heatmaps/`.

**Goal:** Create reproducible, large GP-by-group heatmaps of mean loading for all 200 GPs across healthy non-thymocyte `organ_simplified` tissues and `annotation_level2` clusters.

**Architecture:** The experiment script loads the shared filtered cell-loading matrix and current Seurat metadata through `load_gp_data()`, restricts cells using the same healthy non-thymocyte definition as the AUC pipeline, aggregates each GP to a group mean, and renders raw and within-GP-max-normalized matrices. Each of the four matrices is independently clustered on rows and columns with complete-linkage Euclidean hierarchical clustering; the script also writes the numeric matrices, group sizes, and resulting orders for inspection.

**Tech Stack:** R, ComplexHeatmap, circlize, ZemmourLib, Seurat metadata, PDF.

---

### Task 1: Define the isolated experiment and healthy reference population

**Files:**
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/README.md`

- [x] **Step 1: Load the shared data and select the exact cell universe.**

```r
source("code/R/setup_data.R")
gp_data <- load_gp_data()

keep_cells <-
  gp_data$seurat_meta_filtered$condition_broad == "healthy" &
  gp_data$seurat_meta_filtered$annotation_level1 != "thymocyte"
L_reference <- gp_data$L_pm_filtered[keep_cells, , drop = FALSE]
meta_reference <- gp_data$seurat_meta_filtered[keep_cells, , drop = FALSE]
stopifnot(ncol(L_reference) == 200L, nrow(L_reference) == nrow(meta_reference))
```

- [x] **Step 2: Aggregate all 200 GP columns to one mean per observed group.**

```r
mean_loading_by_group <- function(L_mat, labels) {
  labels <- droplevels(factor(as.character(labels)))
  group_sums <- rowsum(L_mat, group = labels, reorder = TRUE)
  group_counts <- as.integer(table(labels)[rownames(group_sums)])
  list(
    matrix = t(sweep(group_sums, 1L, group_counts, "/")),
    counts = data.frame(group = rownames(group_sums), n_cells = group_counts)
  )
}
```

- [x] **Step 3: Implement row-max normalization without changing the raw values.**

```r
normalize_by_gp_max <- function(mean_matrix) {
  row_max <- apply(mean_matrix, 1L, max)
  stopifnot(all(is.finite(row_max)), all(row_max > 0))
  sweep(mean_matrix, 1L, row_max, "/")
}
```

### Task 2: Render and record the four clustered heatmaps

**Files:**
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_raw_mean_loading_heatmap.pdf`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_row_normalized_mean_loading_heatmap.pdf`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/annotation_level2_raw_mean_loading_heatmap.pdf`
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/annotation_level2_row_normalized_mean_loading_heatmap.pdf`

- [x] **Step 1: Use the canonical annotation palettes and stop if a displayed label lacks a color.**

```r
palette_for_groups <- function(groups, palette) {
  missing <- setdiff(groups, names(palette))
  if (length(missing) > 0L) {
    stop("Missing palette entries: ", paste(missing, collapse = ", "))
  }
  palette[groups]
}

organ_colors <- palette_for_groups(colnames(organ_matrix), ZemmourLib::immgent_colors$organ_simplified)
level2_colors <- palette_for_groups(colnames(level2_matrix), ZemmourLib::immgent_colors$level2)
```

- [x] **Step 2: Render each matrix with complete-linkage Euclidean row and column clustering.**

```r
ht <- ComplexHeatmap::Heatmap(
  display_matrix,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete"
)
drawn_ht <- ComplexHeatmap::draw(ht)
```

- [x] **Step 3: Write the raw matrix, normalized matrix, group-size table, and the row/column order returned by each draw.**

```r
write.csv(cbind(GP = rownames(display_matrix), display_matrix), matrix_csv, row.names = FALSE)
write.csv(group_counts, count_csv, row.names = FALSE)
write.csv(data.frame(GP = rownames(display_matrix)[ComplexHeatmap::row_order(drawn_ht)]), row_order_csv, row.names = FALSE)
write.csv(data.frame(group = colnames(display_matrix)[ComplexHeatmap::column_order(drawn_ht)]), column_order_csv, row.names = FALSE)
```

### Task 3: Validate visual and numeric output

**Files:**
- Create: `experiments/healthy_nonthymus_mean_loading_heatmaps/heatmap_summary.txt`
- Modify: `experiments/healthy_nonthymus_mean_loading_heatmaps/README.md`

- [x] **Step 1: Run the renderer from the repository root.**

Run: `Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R`

Expected: exit status 0; four nonempty PDFs and all CSV/summary files exist.

- [x] **Step 2: Verify dimensions and normalization.**

Run: `Rscript -e 'x <- read.csv("experiments/healthy_nonthymus_mean_loading_heatmaps/organ_simplified_row_normalized_mean_loading_matrix.csv", check.names = FALSE); stopifnot(nrow(x) == 200L, all(abs(apply(x[-1L], 1L, max) - 1) < 1e-10)); cat("OK\\n")'`

Expected: `OK`; equivalent checks pass for the level-2 normalized matrix, and both raw matrices have 200 GP rows.

- [x] **Step 3: Render every PDF to PNG and inspect for clipped labels or unreadable cells.**

Run: `mkdir -p /private/tmp/gp_mean_loading_heatmaps && pdftoppm -png -r 120 experiments/healthy_nonthymus_mean_loading_heatmaps/annotation_level2_row_normalized_mean_loading_heatmap.pdf /private/tmp/gp_mean_loading_heatmaps/level2_normalized`

Expected: one rendered page with all 200 GP labels and all observed level-2 labels visible.

### Task 4: Record the project handoff

**Files:**
- Create: `log/2026-07-13-healthy-nonthymus-gp-mean-loading-heatmaps.md`

- [x] **Step 1: Record the reference population, clustering method, generated files, and completed validation in the update log.**

```markdown
## Validation

- Renderer: `Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R` exited successfully.
- Numeric checks: all raw matrices have 200 GP rows; each normalized GP row has maximum 1.
- Visual checks: all four PDF pages were rendered to PNG and inspected.
```
