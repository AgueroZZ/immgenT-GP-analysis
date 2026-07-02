# code/

Reusable R code backing `script/`. Nothing here writes to
`figures/` or the original `code/`/`script/` directories — see
`../script/README.md` for how the figure scripts use this.

## R/ — sourced helpers

| File | Used by | Contents |
|---|---|---|
| `setup_data.R` | (generic) | `load_gp_data()`: consolidated L/F/metadata loading, replacing the ~15-line block duplicated at the top of every original `Figure_*.R`. Not all figure scripts use it directly (several load slightly different subsets inline), but it's the reference pattern. |
| `roc_auc.R` | Figure2, Figure4, pipeline/02 | `compute_auc_by_group()`, `compute_auc_matrix()`, `compute_auc_threshold_by_group()`, `compute_auc_threshold_matrix()` — one-vs-rest AUC/threshold per GP. From `ROC.R`, unchanged. |
| `plot_utils.R` | Figure1, Figure2, FigureS2 | `scale_cols()`, `filter_cells_by_total_membership()`, `tukey_outliers()`, `lineage_colors()`. |
| `lineage_plots.R` | Figure2 | `plot_gp_swarm()`, `plot_loadings_on_mde()`. |
| `volcano_helpers.R` | Figure1, Figure2 | `plot_gp_signature_volcano()` ("signature volcano" per-gene view). **Not from `script/`/`code/`** — ported from the separate `immgen-signature` Shiny app (`.../Immgen/webapps/immgen-signature/R/mod_signature.R`), a sibling project outside this repo. See the file header for the full provenance note and the caveat about its `n_label` parameter. |
| `cross_gp_helpers.R` | Figure2 | `plot_cross_gp_heatmap()` (multi-GP top-feature comparison heatmap). **Not from `script/`/`code/`** — ported from the same `immgen-signature` app (`R/mod_cross_gp.R`'s `.build_heatmap()`), used for Figure 2 panel 2M. Has one addition beyond the app: an optional `pin_top` argument to manually force specific rows to the top (used to pin Fcer1g/Ccl5/Cd7 in 2M) — this is our own extension, not part of the app. |
| `activation_shared_setup.R` | Figure3, FigureS3 | CD4/CD8 resting-vs-activated cell groups, the curated activation-GP set + semantic color grouping, `F_pm_filtered_norm`, and the GP-group ordering used by both figures' heatmaps. |
| `tf_network.R` | Figure3, FigureS3 | `optimize_bipartite_order()`, `plot_tf_gp_network_v2()` (bipartite TF-GP network layout/plot). |
| `gated_protein_helpers.R` | Figure6, FigureS6 | `MyDimPlotHighlightDensity_df()`, `plot_gated_gp_vs_protein()` (protein-gate vs. GP-loading comparison). |
| `citeseq_shared_setup.R` | Figure6, FigureS6 | CITE-seq cell/protein filtering, curated marker-override table (`df_markers2`), cell exclusions, `well_aligned_gps`. |

## pipeline/ — data preparation (upstream of the figure scripts)

Numbered in dependency order. **These are provenance/documentation, not a
one-command rebuild** — several steps are cluster-scale (need the full
~683k-cell raw Seurat object and packages like `flashier`) and are not
meant to be run against the trimmed `data/` directory shipped in this repo,
which already contains their cached outputs.

1. `01_extract_data.R` — raw Seurat object -> gene-filtered RNA counts,
   CITE-seq protein matrix, condensed factorization summary.
2. `02_compute_auc.R` — per-GP AUC for predicting lineage/organ (feeds
   Figure2, Table S1).
3. `03_protein_thresholds.R` — GMM + manual protein positivity thresholds
   (feeds Figure6, FigureS6).
4. `04_protein_projection.R` — projects CITE-seq protein data onto the
   fixed scRNA cell loadings (Figure 6's caption panel (a) schematic).
5. `05_igt_validation.R` — per-batch (IGT) reproducibility validation
   (feeds FigureS1).

### Known gaps (flagged inline in each script, repeated here for visibility)

- **`flashier_snmf.rds`** (the raw flashier semi-NMF fit) and
  **`flashier_snmf_matrix.qs`** (its input matrix) are read by
  `01_extract_data.R` / `05_igt_validation.R` but produced by no script in
  this repository — the model-fitting step was run interactively/on a
  cluster and that script wasn't preserved here. You (the repo owner)
  mentioned you have this script and will add it separately.
- **`protein_flash_selected_summary_lognorm_backfit200.rds`** (consumed by
  Figure6/FigureS6) is not reproduced exactly by `04_protein_projection.R`
  — that script produces the same object without the `_backfit200` suffix
  (presumably a longer-`maxiter` rerun that wasn't saved back under a
  script).
- **`level_1_AUC_list_figure_no_thymocytes_healthy.rds`** and its
  `organ_simplified` counterpart (consumed by Figure4) are a differently-named
  variant of `02_compute_auc.R`'s output, despite that script already
  restricting to non-thymocyte/healthy cells.

None of these gaps block reproducing the 9 figures from the `data/` files
already in this repo — they only matter if regenerating those cached
files from more raw inputs.

### Explicitly out of scope

The pre-existing `code/` validation suite (`halves_validation*.R`,
`quarter_validation.R`, `eighth_validation.R`, `reproducibility.R`,
`00_/01_bi_cross_validation.R`, `compare_K200_K300.R`,
`replicate_RQVI_*.R`, and the exploratory pre-production fitting scripts
`fit_nmf*.R`/`run_irlba.R`) was confirmed (via grep across all figure
scripts) to feed **no panel in any of the 9 figures reproduced here** — it
supports methods-text robustness claims, not a figure. Left untouched in
`code/`, per your scoping decision.
