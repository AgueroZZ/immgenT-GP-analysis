# CITE-seq Normalization and Fixed-loading EBMF Implementation Plan

> **For agentic workers:** Implement the checked steps in order. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restore the reproducible 20260206 CITE-seq protein-normalization and fixed-loading EBMF workflow, and document it as a dedicated Methods page.

**Architecture:** One cluster-oriented script will extract raw, CLR-normalized, and LogNormalize ADT matrices plus a versioned metadata snapshot from the 20260206 ADT-only Seurat object. A second script will align measured CITE-seq cells to the scRNA GP loading matrix, initialize protein scores by OLS, fix the cell loadings, and reproduce the 20/40/80/120/160/200 flashier checkpoints. The workflowr Methods page and provenance table will distinguish CLR from the LogNormalize matrix actually used for EBMF.

**Tech Stack:** R, Seurat, Matrix, flashier, workflowr

---

### Task 1: Restore 20260206 protein-matrix preparation

**Files:**
- Create: `code/other/prepare_citeseq_protein_matrices_20260206.R`

- [x] Read `igt1_96_withtotalvi20260206_clean_ADTonly.Rds` and validate the ADT assay plus `nCount_ADT`, `cellID`, and `cite_seq` metadata.
- [x] Save `seurat_meta_20260206.rds` and the transposed raw ADT matrix as `protein_mat.rds`.
- [x] Run ADT CLR normalization with `margin = 2` and save the transposed result as the legacy-named `protein_mat_normalized.rds`.
- [x] Compute `adt_scale_factor <- round(mean(seurat_meta$nCount_ADT))`, run ADT LogNormalize, and save `protein_mat_normalized_lognorm.rds`.
- [x] Validate matching dimensions and dimnames across all three matrices and report the observed scale factor (3472 for the audited data).

### Task 2: Restore the fixed-loading protein EBMF checkpoints

**Files:**
- Create: `code/other/fit_citeseq_fixed_loading_ebmf_20260206.R`

- [x] Load `seurat_meta_20260206.rds`, `protein_mat_normalized_lognorm.rds`, and `flashier_snmf_summary.rds`.
- [x] Restrict to `cite_seq == TRUE` cells present in both the protein matrix and scRNA loading matrix, preserving identical row order.
- [x] Compute the OLS protein-score initialization and save `protein_projection_OLS_lognorm.rds`.
- [x] Initialize flashier with point-Laplace protein scores and fix every scRNA-derived loading column.
- [x] Run and save checkpoints at cumulative iterations 20, 40, 80, 120, 160, and 200, using extrapolation only from iteration 80 onward.

### Task 3: Add the CITE-seq Methods page and provenance links

**Files:**
- Create: `analysis/Methods_FlashierFit_Citeseq.Rmd`
- Modify: `analysis/index.Rmd`
- Modify: `code/README.md`
- Modify: provenance comments in scripts that currently label the LogNormalize/backfit200 files as gaps.

- [x] Explain the 20260206 ADT-only input, the distinction between CLR and LogNormalize, and the exact 3472 scale factor.
- [x] Document the fixed-loading model, OLS initialization, flashier prior/variance choices, checkpoint schedule, and outputs.
- [x] Link the new Methods page from the site index and link the two source scripts from the Methods page.
- [x] Replace the resolved `_lognorm` and `backfit200` provenance-gap statements with the restored producers.

### Task 4: Validate and render

- [x] Parse both R scripts without executing the cluster-scale analysis.
- [x] Re-run a sampled numerical reconstruction check against the cached LogNormalize matrix and require floating-point agreement.
- [x] Render `Methods_FlashierFit_Citeseq.Rmd` and `index.Rmd` into `docs/`.
- [x] Inspect the rendered page for headings, code blocks, formulas, filenames, and links.
- [x] Record files, evidence, and validation results in `log/2026-07-14-citeseq-normalization-ebmf-methods.md`.

### Task 5: Simplify the CITE-seq Methods narrative

**Files:**
- Modify: `analysis/Methods_FlashierFit_Citeseq.Rmd`
- Modify: `docs/Methods_FlashierFit_Citeseq.html`

- [x] Reframe the page around the executable workflow: read the provided Seurat object, normalize ADT measurements, align cells, initialize by OLS, and fit fixed-loading EBMF.
- [x] Remove object-version and Seurat-object background exposition from the rendered Methods narrative.
- [x] Render and inspect the revised HTML structure, then update the work log.
