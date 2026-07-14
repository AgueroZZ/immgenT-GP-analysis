# 2026-07-14: Restored CITE-seq normalization and fixed-loading EBMF methods

## Scope

Restored the missing 20260206 CITE-seq protein-matrix preparation and the exact
fixed-loading flashier checkpoint workflow. Added a dedicated workflowr Methods
page and resolved the corresponding `_lognorm` and `backfit200` provenance gaps.
The cluster-scale scripts were documented and syntax-validated but were not run
end-to-end locally.

## Added source files

- `code/other/prepare_citeseq_protein_matrices_20260206.R`
- `code/other/fit_citeseq_fixed_loading_ebmf_20260206.R`
- `analysis/Methods_FlashierFit_Citeseq.Rmd`
- `docs/Methods_FlashierFit_Citeseq.html`
- `plan/2026-07-14-citeseq-normalization-ebmf-methods.md`

## Restored normalization

The authoritative input is
`igt1_96_withtotalvi20260206_clean_ADTonly.Rds`. The preparation script saves:

- `seurat_meta_20260206.rds`;
- raw cells-by-proteins ADT counts as `protein_mat.rds`;
- CLR-normalized ADT (`margin = 2`) under the historical filename
  `protein_mat_normalized.rds`;
- Seurat LogNormalize ADT under
  `protein_mat_normalized_lognorm.rds`.

The LogNormalize scale factor is `round(mean(nCount_ADT))`, which is 3472 for
the 20260206 object. A validation run reconstructed 10,000 sampled cells x 180
proteins directly from the cached raw counts; the maximum absolute difference
from the cached LogNormalize matrix was `8.881784e-16`.

## Restored protein EBMF

The fixed-loading script uses only `cite_seq == TRUE` cells shared by the
20260206 metadata, LogNormalize protein matrix, and scRNA GP loading matrix. It
initializes protein scores by QR-based OLS, fixes all scRNA-derived cell loading
columns, uses point-Laplace protein scores with `var_type = 2`, and saves full
plus compact checkpoints at 20, 40, 80, 120, 160, and 200 cumulative backfit
iterations. Extrapolation is disabled at 20/40 and enabled at 80/120/160/200.

## Version audit

The full and ADT-only 20260206 objects have identical 682,935 cell IDs,
682,935 x 66 metadata, 180 ADT feature identities/dimensions, and reduction
names; the ADT-only object omits the RNA assay. A 20,000-cell comparison also
confirmed that the cached raw protein matrix matches the current ADT counts
exactly for shared cells. The older cached matrix contains 16 additional cells,
all removed safely by the explicit current-metadata intersection.

## Documentation and validation

- Both new R scripts passed `parse()`.
- `workflowr::wflow_build()` successfully rendered the new Methods page and
  rebuilt the site index.
- A real-browser inspection confirmed the long title, floating TOC, MathJax
  formulas, source links, code blocks, and checkpoint table render correctly.
- `code/README.md`, Figure 1, Figure 6, Figure S2, Figure S6, and the shorter
  historical protein-projection pipeline now point to the restored producers.

## Methods-page narrative revision

The CITE-seq Methods page was revised to follow the structure of the upstream
`Methods_FlashierFit` page. It now begins with reading a provided Seurat object
and presents five executable steps: ADT normalization, cell alignment, OLS
initialization, fixed-loading flashier fitting, and the backfitting schedule.
The rendered page no longer describes the identity or reduced-copy status of a
specific Seurat object, nor the historical CLR matrix that is outside the EBMF
fit path. `workflowr::wflow_build()` completed successfully, and the rendered
HTML contains the expected five headings, source-script links, equations, and
checkpoint table. Direct in-app-browser inspection of the local `file://` page
was blocked by browser policy after rendering.
