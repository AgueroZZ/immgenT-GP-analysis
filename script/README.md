# script/

One R script per published figure, plus one per Extended Data table. Each script:

- reads only from `data/` (never writes there),
- writes its panels into `figures/final-selected/Figure N/` (or `Figure SN/`)
  under the *current* figure numbering, which for Figures 1-5 is NOT the
  numbering `figures/Previous/bits/` still uses -- see "Which published panel
  is the ground truth" below before diffing anything (the
  `ExtendedDataTable*.R` scripts instead write a single CSV into
  `figures/final-selected/` -- except Table 8, which writes an .xlsx, since its
  per-GP meta-column grouping needs a merged header row a flat CSV can't
  express),
- sources shared plotting/data-loading helpers from `code/R/`
  rather than redefining them.

Run any script from the repo root, e.g. `Rscript script/Figure4.R`.

## Figure -> script map

| Figure | Script | Primary original source(s) |
|---|---|---|
| 1 | `Figure1.R` | `Figure_Overview.R` (+ panel 1A from `Figure_Lineage.R`) |
| 2 | `Figure2.R` | `Figure_Overview.R` |
| 3 | `Figure3.R` | `Figure_Lineage.R` |
| 4 | `Figure4.R` | `Figure_Activation.R` |
| 5 | `Figure5.R` | `Figure_Organ.R` |
| 6 | `Figure6.R` | `Figure_CITEseq.R` + `gated_protein_loading_plot.R` |
| S1 | `FigureS1.R` | `Figure_Saturation.R` + `Figure_batch.R` |
| S2 | `FigureS2.R` | `Figure_Lineage.R` |
| S3 | `FigureS3.R` | `Figure_Activation.R` |
| S4 | `FigureS4.R` | new (no published counterpart) |
| S5 | `FigureS5.R` + `FigureS5_plot.py` | new (no published counterpart) |
| S6 | `FigureS6.R` | `gated_protein_loading_plot.R` |

Figures 7 and S7 are **out of scope** (see below).

## Which published panel is the ground truth

`figures/Previous/bits/` is the published panel set and still uses the OLD
numbering. Figure 1 was split into Figures 1+2 and everything from the old
Figure 2 onward shifted up by one, so a same-letter filename comparison is
wrong for Figures 1-5:

| ours | published (`figures/Previous/bits/`) |
|---|---|
| 1A, 1B | 1A, 1B |
| 1C | *(new panel -- the GP network, no published counterpart)* |
| 1D | **1C** (the giant loading heatmap; ours adds the highlight bar/boxes) |
| 2A-2F | **1D-1I** |
| 3A-3M | **2A-2M** (letter for letter) |
| 4a-4e | **3c-3g** (the old 3a/3b moved to Figure S3's reserved a/b slots) |
| 5a-5e | **4a-4e** |
| 6*, 7*, S1A-S1D, S2*, S3c-S3h, S6* | same letters |
| S1E, S4*, S5* | *(new panels)* |

There is no `Previous/bits/Figure 5/` at all.

## Extended Data table -> script map

Each table is written as one CSV into `figures/final-selected/` by its script and
previewed on a matching `analysis/ExtendedDataTable*.Rmd` page. Together these
replace the retired Table S1/S2/S3, all now regenerated from code. Table 8 is
the exception: it ships in two layouts from two scripts, a tidy long table
(the default, `ExtendedDataTable8.html`) and a per-GP-column wide table
(`ExtendedDataTable8_wide.html`) that needs a merged-header .xlsx instead of a
flat CSV -- see the row below.

| Extended Data table | Script | Primary original source(s) |
|---|---|---|
| 1 Summary of GP characteristics and annotations | `ExtendedDataTable1_GP_summary.R` | `Supplement_Table1.R` (reworked: healthy non-thymocyte AUC, protein signatures added) |
| 2 GP AUC lineage | `ExtendedDataTable2_GP_AUC_lineage.R` | `runAUC.R` (healthy AUC, via `02_compute_auc.R`) |
| 3 GP AUC tissue | `ExtendedDataTable3_GP_AUC_tissue.R` | `runAUC.R` (healthy AUC, via `02_compute_auc.R`) |
| 4 GP AUC cluster | `ExtendedDataTable4_GP_AUC_cluster.R` | `runAUC.R` (healthy AUC, via `02_compute_auc.R`) |
| 5 GP during activation | `ExtendedDataTable5_GP_during_activation.R` | `Figure_Activation.R` (`GP_activation_summary`) |
| 6 Protein factor matrix | `ExtendedDataTable6_protein_factor_matrix.R` | `Figure_CITEseq.R` (`Protein_F_pm`) |
| 7 Protein gating | `ExtendedDataTable7_protein_gating.R` | `Figure_CITEseq.R` (`compute_alignment_scores`) |
| 8 Comprehensive gene signature matrix (long, default) | `ExtendedDataTable8_gene_signature_matrix_long.R` | new (full list behind Table 1's Top Genes +/-, capped at 100 per direction) |
| 8 Comprehensive gene signature matrix (wide) | `ExtendedDataTable8_gene_signature_matrix_wide.R` | same data as the long version, per-GP-column .xlsx layout |

## Verification status

**How to verify.** Comparing PDF *file sizes* is not a verification -- it once
passed a panel that plotted the wrong variable (5d, 0.96% off in size) and a
gallery whose panels were permuted (6d-6f). Nor is md5: every regenerated PDF
embeds a fresh `/CreationDate`, so only the unmodified copies (1B, 6a, all of
Figure 7) ever match byte-for-byte. Compare **pixels**, against the panel the
mapping above says is the counterpart:

```
magick -density 110 "figures/final-selected/Figure 5/5d.pdf[0]" \
  -background white -alpha remove -colorspace Gray a.png
magick -density 110 "figures/Previous/bits/Figure 4/4d.pdf[0]" \
  -background white -alpha remove -colorspace Gray b.png
magick compare -metric RMSE a.png b.png null:
```

Note if you script this in zsh: the page selector must be written
`"${pdf}[0]"`, since `"$pdf[0]"` is parsed as an array subscript and fails with
a confusing "no decode delegate" error.

**Current state.** Every script was re-run and every panel pixel-compared
against its published counterpart. All panels reproduce the published figure
except where listed below. Nothing is left unexplained.

### Deliberate deviations from the published panels

- **Figure 1, panel 1C** is a new panel (the GP-gene network); the published
  1C is our 1D, which additionally carries the top color bar, colored column
  labels and boxes marking the GPs highlighted in 1C.
- **Figure 2, panels 2B/2C**: the "active cell / active gene" definitions were
  replaced with hard thresholds (normalized loading > 0.1; |normalized score| >
  0.25) and the axes are now log-scaled, so these no longer resemble the
  published 1E/1F linear base-R histograms. The captions state the current
  definitions. `FigureS1.R`'s S1E uses the same definitions.
- **Figure 3, panel 3M** is a genuine "top 30 up-regulated genes" panel: it
  ranks candidates by max(score) across GP3/GP29/GP22 (`rank_by = "pos"`),
  where the published 2M ranked by max|score| and so admitted 9 genes on a
  large *negative* score, each let past the "positive somewhere" gate by a
  token +0.01..+0.12 elsewhere (CT010467.1, Cmss1, Tmsb4x, Ms4a4b, Cd52,
  Mir6236, Ppia, Malat1, Ly6e -- the blue block at the bottom of 2M). Their
  replacements are the next 9 genuinely up-regulated genes: Anxa2, Prkch,
  Itm2b, Klra1, Cd3e, Nr4a2, Klre1, Chn2, Klrk1, which brings GP22's
  NK-receptor set (Klra1/Klre1/Klrk1, next to the Klra7/Klrd1 already there)
  into the figure. 21 of the 30 are unchanged. Rows are additionally pinned
  (Fcer1g/Ccl5/Cd7/Ctsw on top) and the columns ordered GP3/GP29/GP22.
- **Figure 5, panel 5c** drops the across-organ expression dotplot that formed
  the left half of the published 4c and ships the gene-score heatmap alone, at
  half the width. The caption describes the heatmap only.
- **Figure 6, panel 6b** was reworked into the triangular-first protein-row
  ordering described in its caption; the published 6b is a plain clustered
  heatmap with a different aspect ratio.
- **Figure S3, panel s3c** now reports activated CD4/CD8 cells at threshold
  0.1, not all cells at 0.2 as the published panel did. The caption states it.
- **Figure S1, panel S1C** is computed over the 18 standard-spleen IGTs that
  contribute >= 500 cells, not all 35, so that it shares one per-IGT
  mean-loading matrix with S1D and its ten labelled GPs are exactly S1D's ten
  rows. `Figure_batch.R` built that matrix twice on different IGT sets, and the
  two rankings disagree (GP2: 7th on the subset, 11th over all 35; GP25: 6th
  over all 35, 23rd on the subset), so the published S1C labels GP25 where the
  published S1D shows GP2 -- "the top ten from (C)" could not be followed
  across the panels. Restricting both to the 18 also drops GP25 out of the top
  ten, which is the point: its rank over all 35 came from IGTs with very few
  cells. Note the consequence -- **S1D is pixel-identical to the published
  panel**, and S1C is the panel that moves (its variance axis roughly doubles,
  since it is no longer diluted by small IGTs).
- **Figures 3C/3E/3G** and **s3c** had no exact source in the original scripts
  and were reconstructed (3C/3E/3G with `plot_loadings_on_mde()`, the styling
  used for 3J/3K/3L). Content matches; canvas size and the UMAP->MDE / K->GP
  relabeling differ.
- **Figure 3, panels 3J/3K/3L** are lettered GP3 = J, GP29 = K, GP22 = L, so
  the figure keeps one GP order throughout (it matches 3M's GP3/GP29/GP22
  columns). The published 2J/2K/2L were GP22/GP29/GP3, so **J and L are
  swapped**: our 3J is the published 2L and our 3L is the published 2J. The
  panel *contents* are unchanged apart from the UMAP->MDE relabeling below, and
  `verify_panels.sh` pairs them by content (3J<->2L, 3L<->2J) so the swap does
  not read as a regression.
- **Figure 3, panel 3I** and the J-L panels are otherwise identical to 2I-2L
  apart from the UMAP->MDE axis/legend relabeling, and **4e**'s title was
  renumbered.

### Panels not code-generated here

- **1B** and **6a** are hand-finished schematics (confirmed via PDF `/Creator`
  metadata: the only two non-R-generated panels among the 9 figures). They are
  copied through untouched -- no script writes them.

### Panels from the `immgen-signature` Shiny app, not the analysis scripts

**Figure 2, panel 2A** and **Figure 3, panels 3D/3F/3H** (GP1/GP68/GP30/GP58
signature volcanoes) and **Figure 3, panel 3M** come from a sibling project,
`.../Immgen/webapps/immgen-signature/`, not from the pre-refactor analysis
code. The plotting functions were ported into `code/R/volcano_helpers.R`
(`plot_gp_signature_volcano()`, matching `mod_signature.R::.make_static_plot()`)
and `code/R/cross_gp_helpers.R` (`plot_cross_gp_heatmap()`, matching
`mod_cross_gp.R`'s `.build_heatmap()`).

Both are driven by sidebar inputs whose exported values are not recorded in any
file, so they had to be recovered by matching the published panels:

- the volcanoes' `n_label` (number of labeled genes) was bisected per panel;
  GP68 and GP58 match exactly, GP1 and GP30 are within ~30 bytes.
- 3M's `direction` is **"pos"** ("Positive only"), not the app's `"both"`
  default. This is what the published 2M was exported with: under `"both"` the
  top-30-by-max|score| pulls in 11 all-negative housekeeping genes
  (Tpt1/Actb/Eef1a1/Tmsb10/...) and pushes out
  Ikzf2/Klrd1/Il2rb/Itgae/Dapk2/Junb/Cd3g/Ly6e/Ppia/Malat1/Mir6236. With
  `direction = "pos"` and the app's ranking, our port reproduced the published
  30-gene set exactly -- that is how the setting was recovered.

  3M then deliberately moves off it via `rank_by = "pos"` (our own extension to
  `plot_cross_gp_heatmap()`, see the deviations list above). The app has no such
  control: it always ranks on max|score|, so "Positive only" gates which genes
  are eligible without affecting the order they are ranked in, and a strongly
  down-regulated gene can be selected as long as it is positive somewhere.
  `experiments/fig3m_updown_ranking/` has the three-way comparison
  (max|score| gated / max(score) / max|score| ungated) that settled it.

### Reproducibility

Every RNG consumer is seeded, so re-running a script reproduces its panels
bit-for-bit (modulo the PDF creation date):

- cell subsampling / MDE sampling: `set.seed()` at the top of each panel's
  section, carried over from the original scripts;
- swarm jitter: `position_jitter(seed = )` via `plot_gp_swarm(seed = 42)`;
- graph layouts: `set.seed()` before each `ggraph()` call; the TF-GP layout in
  `code/R/tf_network.R` is a deterministic barycenter ordering, no RNG;
- **label repulsion: every `geom_text_repel()` call passes `seed = 42`.**
  Without it `ggrepel` picks a different layout on every run, which is what
  made panels drift between the PDFs and the site PNGs.

`figures/Previous/bits/` was produced before the repel seeds existed, so label
*placement* (never label content or data) still differs slightly from the
published panels in the label-heavy scatterplots (3A, 3D, 3F, 3H, 4a, 5a, 5d,
6g, 6h, S1C). `4b`'s and `s3h`'s node placement likewise differs from the
published layout because of `graphlayouts`/`igraph` version drift, not RNG.

### Site PNGs

`analysis/assets/**.png` are **converted from the panel PDFs** by
`script/render_site_assets.sh`, not re-plotted. Re-plotting was how the site
and the published panels silently diverged. Run it after any script re-run and
before `workflowr::wflow_build()`.

### Panels rasterized by `scattermore`

**6c-6f** (protein-gate vs. GP-loading galleries) draw a rasterized density
layer, so their PDFs differ from the published ones by 0.5-5% in size and
~0.04 RMSE with no visible difference. Their GP-to-letter assignment is
hard-coded (`fig6_letter`) to the published order c=GP171, d=GP12, e=GP80,
f=GP23 -- do not let it fall back to the order `GPs_fig6` happens to list.

## Note on re-running FigureS6.R

`cairo_pdf()` was found during verification to not reliably truncate an
existing output file of a different size in place (a stale/partial file
from an earlier interrupted run could persist and get silently reused).
`FigureS6.R` now explicitly `unlink()`s each page's output file before
opening the `cairo_pdf()` device, so re-running it is safe. If you ever see
a script's output that doesn't match expectations, deleting its output
files before re-running is good practice in general.

## Why Figure 7 and Figure S7 are excluded

Confirmed via PDF `/Creator` metadata during this refactor:

- Figure 7, panel A is a **Figma** schematic (not code-generated).
- Figure 7, panels C/D/F and **all 11** Figure S7 panels are **Matplotlib**
  (Python) output — but no `.py` file exists anywhere in this repository.
- The remaining Figure 7 R-generated panels (7B, 7E, 7G) trace to
  `replicate_RQVI_cells.R` and related RQVI validation code, but that
  entire validation track was scoped out of this refactor (see
  `code/README.md`).

These are not reproducible from this repository as it stands. If the
missing Python source turns up, it can be added as a `code/python/`
sibling; until then, Figure 7/S7 have no script here.
