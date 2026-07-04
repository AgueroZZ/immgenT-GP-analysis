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
2. `01b_filter_cells.R` — filters cells by total GP membership, producing
   `L_pm_filtered.rds`/`F_pm_filtered.rds` (unlike `01_extract_data.R`,
   this one *is* runnable against the `data/` in this repo).
3. `02_compute_auc.R` — per-GP AUC for predicting lineage/organ, computed
   under two different cell restrictions (healthy-only for Figure2/Figure4;
   non-thymocyte-only, no healthy filter, for Table S1 -- see the "Data
   provenance" table below for why these differ).
4. `03_protein_thresholds.R` — GMM + manual protein positivity thresholds
   (feeds Figure6, FigureS6).
5. `04_protein_projection.R` — projects CITE-seq protein data onto the
   fixed scRNA cell loadings (Figure 6's caption panel (a) schematic), and
   derives `CITEseq_markers_full.rds`.
6. `05_igt_validation.R` — per-batch (IGT) reproducibility validation
   (feeds FigureS1).

## Data provenance

Every `data/*` file read anywhere in `script/` or `code/pipeline/`, and
which step (if any) in this repo produces it. Figure scripts also flag
their own inputs inline via a `Required inputs` header comment that points
back to this table.

`data/` is a symlink into the shared `immgen-t-factors` project; most files
there are raw/primary inputs, but `code/pipeline/` scripts also write their
outputs there directly (e.g. the AUC files below) -- none of `data/` is
pushed to GitHub either way (it's local-only, not part of what this repo
shows readers), so there's no need to separate "input" from "generated"
across two directories.

| `data/` file | Produced by | Notes |
|---|---|---|
| `igt1_96_withtotalvi20260206_clean_ADTonly.Rds` | *(primary input)* | Processed Seurat object (RNA + CITE-seq ADT); not produced by anything in this repo, it's the starting point. |
| `flashier_snmf.rds` | *(external, to be provided separately -- e.g. Zenodo)* | Raw flashier semi-NMF fit (EBMF). Fitting this takes a long time; the fitting script isn't in this repo yet, and the plan is to distribute this file directly rather than the code that produces it. |
| `flashier_snmf_matrix.qs` | *(derivable from `flashier_snmf.rds`)* | `flashier_snmf$flash_fit$Y` -- the input matrix to the fit above, saved alongside it (see `code/extract_data.R`, the actual extraction script run on the source data, not yet ported into `code/pipeline/`). |
| `flashier_snmf_fitted_prior.rda` | *(external, to be provided separately -- same as `flashier_snmf.rds`)* | Saved alongside `flashier_snmf.rds` from the same fit. Used by Figure1.R's PVE-type panel. |
| `flashier_snmf_summary.rds` | `01_extract_data.R` | Condenses `flashier_snmf.rds` into `L_pm`/`F_pm`/`elbo`/`pve`. |
| `L_pm_filtered.rds`, `F_pm_filtered.rds` | `01b_filter_cells.R` | `flashier_snmf_summary.rds`'s `L_pm`/`F_pm`, restricted to cells present in the current Seurat object (18 of `flashier_snmf_summary.rds`'s 682,953 cells aren't -- the same kind of data-version drift already noted for `seurat_meta.rds`), then filtered by `filter_cells_by_total_membership()` (`code/R/plot_utils.R`). Recovered from a commented-out block in the original `Figure_Overview.R`; verified **byte-identical** to the cached files (max abs diff 0, same 681,423 cells in the same order). |
| `shifted_log_counts.qs`, `counts.qs` | `01_extract_data.R` | Gene-filtered RNA counts (raw and shifted-log). |
| `shifted_log_counts_subset.rds` | **gap** | A subset of `shifted_log_counts.qs` (used by Figure4's panel c); the subsetting script isn't preserved. |
| `mean_shifted_log_expr.rds` | **gap** | Per-gene mean shifted-log expression, most likely `colMeans()` of `shifted_log_counts.qs`; not scripted here. |
| `protein_mat.rds`, `protein_mat_normalized.rds` | `01_extract_data.R` | CITE-seq ADT matrix (raw and CLR-normalized). |
| `protein_mat_normalized_lognorm.rds` | **gap** | A differently-named/derived variant of `protein_mat_normalized.rds`; the exact transform isn't preserved. |
| `Thresholds_Selected_Proteins.csv`, `GMM_Thresholds_Summary.csv` | `03_protein_thresholds.R` | Per-protein positivity thresholds. |
| `TableS4_citeseq_qc_20250513.csv` | *(external)* | The manuscript's own Supplementary Table S4 (manually reviewed protein QC classifications) — not computationally derived. |
| `CITEseq_markers_full.rds` | `04_protein_projection.R` | Per-GP positive/negative protein markers (\|score\| >= 0.5 on the same filtered/scaled protein factor matrix as Figure 6 panel b). Recovered from live (not commented-out) code in the original `Figure_CITEseq.R`; the marker-selection logic itself is verified byte-identical to the cached file when fed the same upstream input -- the only divergence is that this step uses the non-`backfit200` protein summary (see that row's caveat above). |
| `protein_flash_selected_summary_lognorm.rds` | `04_protein_projection.R` | Re-estimated protein factor matrix U (Figure 6 panel a schematic). |
| `protein_flash_selected_summary_lognorm_backfit200.rds` | `citeseq_run.R` (not yet ported into `code/pipeline/`) | Six-stage `flash_backfit()` sequence (20+20+40+40+40+40 iterations, checkpointed at each stage) starting from `04_protein_projection.R`'s OLS-initialized fit; the final checkpoint is this file. Doesn't depend on cell-level annotation (protein values/embeddings only), so not affected by the stale-`seurat_meta.rds` issue found for the AUC files below. |
| `level_1_AUC_list_figure_no_thymocytes_healthy.rds`, `level_2_AUC_list_figure_no_thymocytes_healthy.rds`, `organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds` | `02_compute_auc.R` | Per-GP AUC/threshold for predicting level-1/level-2/organ, restricted to non-thymocyte **healthy-only** cells. Used by Figure2.R (2A) and Figure4.R. |
| `level_1_AUC_list_figure_no_thymocytes.rds`, `level_2_AUC_list_figure_no_thymocytes.rds`, `organ_simplified_AUC_list_figure_no_thymocytes.rds` | `02_compute_auc.R` (copied from the old no-suffix files below, not recomputed -- verified equivalent) | Per-GP AUC/threshold for predicting level-1/level-2/organ, restricted to non-thymocyte cells **without** an additional healthy-only restriction (healthy and diseased together). Used by TableS1.R -- confirmed against the published `Table S1.xlsx`: its Organ column includes disease-specific sites (SLO, prostate, pancreas, synovial fluid) that only exist in this broader population, and a full column-by-column diff against the published table matches on Signature genes 200/200, Level1 199/200, Organ 189/200, Level2 169/200 (residual few-row differences are AUC values sitting right at the 0.8 cutoff, consistent with the ~16-18 cell version drift already noted elsewhere in this table). |
| `level_1_AUC_list_figure.rds`, `level_2_AUC_list_figure.rds`, `organ_simplified_AUC_list_figure.rds` (no suffix) | *(superseded by the `_no_thymocytes.rds` files above, same computation, kept only for the ambiguous old name)* | Confirmed identical to the `_no_thymocytes` files above (`organ_simplified`: byte-identical; `level_1`/`level_2`: diff ~1e-12, floating-point noise) -- simply copied under the clearer name rather than recomputed. Earlier documentation here wrongly guessed these were superseded/wrong; they are in fact what TableS1.R should use (see above). Left in `data/` untouched, no longer read directly (scripts read the `_no_thymocytes.rds` copies instead). |
| `condition_detailed_AUC_list_figure.rds` | **gap** | TableS1's published Condition column (per-GP AUC against `condition_detailed`, non-thymocyte cells, no healthy restriction -- same restriction as the `_no_thymocytes` family above, not degenerate for this variable since it's not additionally healthy-filtered) isn't reproduced by anything in this repo. This cached file predates the current Seurat object (dated Dec 2 2025) and was never reverified against current annotation; a fresh run is ~200 GPs x 127 `condition_detailed` categories x 635k cells, timed at ~4.5 min for just 5 GPs (~3 hours for all 200), so it's deferred rather than rerun. `TableS1.R` currently has no Condition column as a result. |
| `umap_result.rds` | *(derivable, not currently scripted)* | Byte-identical to `Seurat::Embeddings(seurat_obj, "mde2_totalvi_20241006")` for the cells present in both. The cached file has 682,951 cells vs. the current Seurat object's 682,935 (16 extra) -- same class of drift as elsewhere, but the values for the overlapping cells match exactly, confirmed. |
| `igt_specific_cosine_scores.csv`, `igt_specific_validated_matrix.csv` | `05_igt_validation.R` (Stage B) | Per-IGT reproducibility score matrix (feeds FigureS1). |
| `igt_specific/*.qs` (per-IGT flashier fits) | `validate_by_experiments.R` (not yet ported into `code/pipeline/`; Stage A) | Already present in `data/`. Per-IGT flashier fits, run separately for each of the ~80 IGT batches via `mclapply()`; depends on `flashier_snmf_matrix.qs` (see above) and is cluster-scale, same as `01_extract_data.R`. |
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
