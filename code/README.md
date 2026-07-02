# code/

Reusable R code backing `script/`. Nothing here writes to `figures/` — see
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

## Data provenance

Every `data/*` file read anywhere in `script/` or `code/pipeline/`, and
which step (if any) in this repo produces it. Figure scripts also flag
their own inputs inline via a `Required inputs` header comment that points
back to this table.

| `data/` file | Produced by | Notes |
|---|---|---|
| `igt1_96_withtotalvi20260206_clean_ADTonly.Rds` | *(primary input)* | Processed Seurat object (RNA + CITE-seq ADT); not produced by anything in this repo, it's the starting point. |
| `flashier_snmf.rds` | **gap** | Raw flashier semi-NMF fit. Fit interactively/on a cluster; that fitting script is not yet in this repo (to be added). |
| `flashier_snmf_matrix.qs` | **gap** | Input matrix to the fit above; same status. |
| `flashier_snmf_fitted_prior.rda` | **gap** | Saved alongside `flashier_snmf.rds` from the same interactive fit; not reproduced here. Used by Figure1.R's PVE-type panel. |
| `flashier_snmf_summary.rds` | `01_extract_data.R` | Condenses `flashier_snmf.rds` into `L_pm`/`F_pm`/`elbo`/`pve`. |
| `L_pm_filtered.rds`, `F_pm_filtered.rds` | **gap** | These are `flashier_snmf_summary.rds`'s `L_pm`/`F_pm` after `filter_cells_by_total_membership()` (`code/R/plot_utils.R`) — the script that ran that filter and saved these exact files is not preserved here. |
| `shifted_log_counts.qs`, `counts.qs` | `01_extract_data.R` | Gene-filtered RNA counts (raw and shifted-log). |
| `shifted_log_counts_subset.rds` | **gap** | A subset of `shifted_log_counts.qs` (used by Figure4's panel c); the subsetting script isn't preserved. |
| `mean_shifted_log_expr.rds` | **gap** | Per-gene mean shifted-log expression, most likely `colMeans()` of `shifted_log_counts.qs`; not scripted here. |
| `protein_mat.rds`, `protein_mat_normalized.rds` | `01_extract_data.R` | CITE-seq ADT matrix (raw and CLR-normalized). |
| `protein_mat_normalized_lognorm.rds` | **gap** | A differently-named/derived variant of `protein_mat_normalized.rds`; the exact transform isn't preserved. |
| `Thresholds_Selected_Proteins.csv`, `GMM_Thresholds_Summary.csv` | `03_protein_thresholds.R` | Per-protein positivity thresholds. |
| `TableS4_citeseq_qc_20250513.csv` | *(external)* | The manuscript's own Supplementary Table S4 (manually reviewed protein QC classifications) — not computationally derived. |
| `CITEseq_markers_full.rds` | **gap** | Curated per-GP marker-protein table (positive/negative signature); presumably manually reviewed, not scripted here. |
| `protein_flash_selected_summary_lognorm.rds` | `04_protein_projection.R` | Re-estimated protein factor matrix U (Figure 6 panel a schematic). |
| `protein_flash_selected_summary_lognorm_backfit200.rds` | **gap** | The variant actually consumed by Figure6/FigureS6 — likely the same script re-run with a larger `maxiter` and saved under a different name; not reproduced exactly. |
| `level_1_AUC_list_figure.rds`, `level_2_AUC_list_figure.rds`, `organ_simplified_AUC_list_figure.rds` | `02_compute_auc.R` | Per-GP AUC/threshold for predicting level-1/level-2/organ, restricted to non-thymocyte healthy cells. |
| `level_1_AUC_list_figure_no_thymocytes_healthy.rds`, `level_2_AUC_list_figure_no_thymocytes_healthy.rds`, `organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds` | **gap** | Consumed by Figure4; a differently-named variant of `02_compute_auc.R`'s output despite that script already restricting to non-thymocyte/healthy cells. Both variants exist side-by-side in `data/`, suggesting a later rerun changed the naming without an updated script being preserved. |
| `condition_detailed_AUC_list_figure.rds` | **gap** | Consumed by TableS1; the same AUC pattern as `02_compute_auc.R` but grouped by `condition_detailed` instead of level_1/level_2/organ — not scripted here, though it would follow the same pattern as the other three. |
| `umap_result.rds` | **gap** | The MDE/UMAP embedding coordinates used throughout. The dimensionality-reduction step itself is not scripted anywhere in this repo. |
| `igt_specific_cosine_scores.csv`, `igt_specific_validated_matrix.csv` | `05_igt_validation.R` (Stage B) | Per-IGT reproducibility score matrix (feeds FigureS1). |
| `igt_specific/*.qs` (per-IGT flashier fits) | **gap** (Stage A) | Already present in `data/`; the per-IGT re-fitting step that produced them is not scripted here (see `05_igt_validation.R`'s header). |
| `GSEA_signatures_select_toplot.csv` | *(external)* | A curated/downloaded gene-set collection (e.g. from MSigDB), not generated by any script here. |

None of these gaps block reproducing the 9 figures from the `data/` files
already in this repo — they only matter if regenerating those cached files
from more raw inputs.

### Explicitly out of scope

A validation suite (split-sample reproducibility, bi-cross-validation, rank
comparison, and pre-production model-fitting scripts) was confirmed via
grep across all figure scripts to feed **no panel in any of the 9 figures
reproduced here** — it supports methods-text robustness claims, not a
figure, so it isn't included in this repository.
