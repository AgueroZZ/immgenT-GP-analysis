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
| 7, panel b only | `Figure7b.R` + `Figure7b_rematch.py` + `Figure7b_plot.py` | new (was Extended Data Figure 5 until 2026-07-28) |
| S1 | `FigureS1.R` | `Figure_Saturation.R` + `Figure_batch.R` |
| S2 | `FigureS2.R` | `Figure_Lineage.R` |
| S3 | `FigureS3.R` | `Figure_Activation.R` |
| S4 | `FigureS4.R` | new (no published counterpart) |
| S5 | `FigureS5.R` | `Figure_CITEseq.R` (its panel a) |
| S6 | `FigureS6.R` | `Figure_CITEseq.R` + `gated_protein_loading_plot.R` |

The rest of Figure 7, and all of S7, are **out of scope** (see below). The S5
slot was emptied on 2026-07-28, when what was published there became main
Figure 7b, and refilled on 2026-07-30 by the protein-program heatmap that had
been Figure S6's first panel.

## The 2026-07-28 / 07-29 re-lettering

Four figures were reordered at once. Each script's header carries its own
old-to-new table; in brief:

- **Figure 6** was reordered to a-j: the KLRG1 pair and the CD69 gene heatmap
  moved to the front (published g, h, i -> b, c, d), the four gating panels moved
  back (published c-f -> e-h), and two gating panels were added, (i) GP77 and
  (j) GP8. The published 6b, 6j and 6k moved out to Figure S6.
- **Figure S6** stopped being a two-page gallery of 23 GPs. It became seven
  panels: s6a-s6c were the published 6b/6j/6k, and s6d-s6g four of the
  gallery's GPs (GP29, GP58, GP22, GP68) as standalone panels. Two more gallery
  GPs were promoted into Figure 6; the remaining 17 are no longer published.
  **Re-lettered again on 2026-07-30** (below): it is now six panels, s6a-s6f.
- **Figure S3**'s s3a and s3b were merged into one s3a (they were always one
  panel), so every panel after them dropped a letter: published s3c-s3h are now
  s3b-s3g.
- **Figure 4** (2026-07-29) dropped its bipartite TF-GP network, the published
  3e, so the two heatmaps after it moved up a letter: published 3f, 3g are now
  4c, 4d and the figure is a-d. That panel was the only caller anywhere of
  `code/R/tf_network.R`, which is therefore now unused -- Figure S3's s3g draws
  its own TF-GP network inline. The helper is kept for provenance, and
  `FigureS3.R`'s `source()` of it (which never called into it) was removed.
- **Extended Data Figure 5** became main **Figure 7b**, assembled from its two
  half-panels into one PDF, replacing a different collaborator panel.

## The 2026-07-30 S5/S6 split

Figure S6's first panel -- the protein-program heatmap, published as Figure 6b --
became a standalone figure, `FigureS5.R` -> `figures/final-selected/Figure S5/s5.pdf`,
taking the Extended Data Figure 5 slot that Figure 7b had vacated. Figure S6's
remaining panels each dropped one letter: its s6b-s6g are now **s6a-s6f** (the
two CD69 mean-activity heatmaps are a, b; the four gating panels c-f). Only
letters and file paths moved; no panel was re-rendered, so every PDF and site PNG
is the same file under a new name.

## Which published panel is the ground truth

`figures/Previous/bits/` is the published panel set and still uses the OLD
numbering. Figure 1 was split into Figures 1+2 and everything from the old
Figure 2 onward shifted up by one, and Figures 6 / S3 / S6 were re-lettered on
2026-07-28 and Figure 4 on 2026-07-29, so a same-letter filename comparison is
wrong almost everywhere:

| ours | published (`figures/Previous/bits/`) |
|---|---|
| 1A, 1B | 1A, 1B |
| 1C | *(new panel -- the GP network, no published counterpart)* |
| 1D | **1C** (the giant loading heatmap; ours adds the highlight bar/boxes) |
| 2A-2F | **1D-1I** |
| 3A-3M | **2A-2M** (letter for letter) |
| 4a, 4b | Figure 3/**3c, 3d** (the old 3a/3b moved to Figure S3's a/b slot) |
| 4c, 4d | Figure 3/**3f, 3g** (the published 3e, a TF-GP network, was dropped from the figure on 2026-07-29 and has no counterpart here) |
| 5a-5e | **4a-4e** |
| 6a | Figure 6/**6a** |
| 6b, 6c | Figure 6/**6g, 6h** (KLRG1) |
| 6d | Figure 6/**6i** (CD69 gene heatmap) |
| 6e-6h | Figure 6/**6c-6f** (gating: GP171, GP12, GP80, GP23) |
| 6i, 6j | *(new panels -- gating for GP77 and GP8)* |
| S5/s5 | Figure 6/**6b** (protein-program heatmap; Figure S6's s6a until 2026-07-30) |
| S6/s6a, s6b | Figure 6/**6j, 6k** (CD69 GPs by tissue / lineage) |
| S6/s6c-s6f | *(new standalone panels; their GPs were only ever drawn inside the retired Figure S6/s6-1 and s6-2 gallery pages)* |
| S3/s3b-s3g | Figure S3/**s3c-s3h** |
| 7A, 7C-7G | Figure 7/**7A, 7C-7G** (straight copies, out of scope) |
| 7B | *(ours since 2026-07-28 -- the former Extended Data Figure 5, assembled; it replaced a different published 7B, which is therefore NOT its counterpart)* |
| S1A-S1D, S2* | same letters |
| S1E, S4* | *(new panels)* |

There is no `Previous/bits/Figure 5/` at all, and no Figure S3 s3a/s3b (that
panel is not produced in this repository).

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
| 7 Manual protein positivity thresholds | `ExtendedDataTable7_protein_thresholds.R` | `protein_thresholding_manual.R` (`Threshold_manual`, via `03_protein_thresholds.R`) |
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

**A re-run can silently not write.** `FigureS6.R` was seen to exit 0, with a
complete log, having regenerated nothing: `unlink()` removed the old
`s6-1.pdf`/`s6-2.pdf` (confirmed gone by an external `stat`), the cairo device
opened and closed without error, and the previous files then reappeared
*byte-identical with their original mtime*. Narrowed to the combination of a
memory-heavy session, `cairo_pdf()`, and a multi-megabyte write to the published
path -- the same script rendering to a scratch directory overwrites happily,
`ggsave()` over a tracked panel PDF is unaffected, and a small `cairo_pdf()`
write to the same directory in the same session is unaffected. Cause is outside
R; treat it as an environment hazard, not a script bug.

The two gallery pages that triggered it were retired on 2026-07-28, so no script
here writes multi-megabyte `cairo_pdf()` output any more. `FigureS6.R` keeps a
cheaper form of the mitigation instead: it stamps `run_started_at` at the top and
`stop()`s at the end if any of its six panels is missing, empty, or older than
that stamp (`FigureS5.R`, split out of it, carries the same guard for its one
panel). If a panel script ever appears to run cleanly without updating its
output, check the output's **mtime** before trusting the run, and reach for
either pattern (the stronger one -- render to `tempfile()`, `file.copy()` into
place, compare sizes -- is in this file's git history at commit `5248bd6`).

**Verifying published *numbers*, not just panels.** Where a table publishes a
value that a figure also computes with, the two have to be checked against each
other -- nothing in the build enforces it. `verify_thresholds.R` does this for
the protein positivity thresholds: it sources the real
`code/R/citeseq_shared_setup.R` (not a copy of its logic), takes the
`threshold_df` that `Figure6.R:170` and `FigureS6.R:55` actually gate on, and
diffs it against the values Extended Data Table 7 publishes -- protein set,
every value, and whether each published threshold is even reachable through
`select_proteins` and the ADT matrix. It also reports the reverse gap: markers
that are gated on but have no threshold and so are silently skipped. Exits
non-zero on any mismatch; both the pass and fail paths are tested.

```
Rscript script/verify_thresholds.R
```

`verify_cd69_gp_ranking.R` does the same for a *claim* rather than a value: the
CD69 GP subset (`cd69_top_gps_subset`, in `code/R/citeseq_shared_setup.R`) is
hand-picked, so the two captions describing it (Fig. 6d and Fig. S6a, b) and the
correlation ranks quoted beside it are all unenforced. The script sources the real
setup, recomputes the CD69 Spearman correlation over all 200 GPs on the figures'
own cells, and checks that the curated 10 really are among the most strongly
correlated (all within the top 20 by |rho|), that they are *not* a true top-10
under either ranking, that both signs are still represented, that the shared
setup's comment quotes the current ranks, that neither figure script keeps its own
copy of the list, and that neither caption has dropped its hedge or drifted back
to the old "the ten GPs most associated with CD69" wording.

```
Rscript script/verify_cd69_gp_ranking.R
```

`verify_gating_gps.R` covers the ten protein-gating panels, which are split
across two figures (6e-6j and s6c-s6f) and so across two scripts. It reads each
script's GP-to-letter map rather than retyping it, and checks that every gated GP
comes from the curated `well_aligned_gps` pool, that the two figures show disjoint
GPs, that every declared letter has a PDF on disk with no orphans left from an
earlier lettering, and -- reproducing the gate's own logic -- that no marker in
those ten signatures is silently skipped for want of an ADT column or a manual
threshold. That last check is what lets Fig. S6c-f's caption and Extended Data
Table 7's caption both state that no panel subtitle omits a marker.

```
Rscript script/verify_gating_gps.R
```

**Current state.** Every script was re-run and every panel pixel-compared
against its published counterpart. All panels reproduce the published figure
except where listed below. Nothing is left unexplained. `verify_thresholds.R`
passes: all 42 published thresholds are identical to the gated ones.
`verify_cd69_gp_ranking.R` and `verify_gating_gps.R` pass.

After the 2026-07-28 re-lettering, every panel that only changed letter was
additionally pixel-compared against its own pre-move file and came out
**identical** (Figure 6's b-h, Figure S6's a-c as they were then, Figure S3's
b-g, and Figure 7b
against the previously assembled Extended Data Figure 5), so the re-lettering
moved letters and nothing else. The one panel whose content changed on purpose is
6i (GP77) -- see the gating-helper note under "Deliberate deviations".

### Deliberate deviations from the published panels

- **Figure 1, panel 1C** is a new panel (the GP-gene network); the published
  1C is our 1D, which additionally carries the top color bar, colored column
  labels and boxes marking the GPs highlighted in 1C.
- **Figure 2, panels 2B/2C**: the "active cell / active gene" definitions were
  replaced with hard thresholds (normalized loading > 0.1; |normalized score| >
  0.25), so the bar heights differ from the published 1E/1F base-R histograms;
  the x axes are linear as they are there (they were log-scaled here until
  2026-07-29, when the collaborators asked for ordinary axes back -- which moved
  both panels *closer* to the published pair, RMSE 0.239 -> 0.154 for 2B and
  0.209 -> 0.160 for 2C). Both use 60 bins, raised from 40 in the same round:
  worth it for 2C (~7 genes wide over a 1..421 range; much past 60 and single
  GPs out of 200 start showing as their own bar), near-invisible for 2B, which
  is dominated by its first bar at any binning because 152 of the 200 GPs are
  active in under 2.5% of cells. The captions state the current definitions.
  `FigureS1.R`'s S1E uses the same definitions -- only the definitions: it keeps
  its log-log axes, which its caption states, since it spans several orders of
  magnitude on both axes.
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
- **Figure S5** (the published Figure 6b, and Figure S6's s6a until 2026-07-30)
  was reworked into the triangular-first protein-row ordering described in its
  caption; the published panel is a plain clustered heatmap with a different
  aspect ratio.
- **Figure 6, panel 6d** (the published 6i, the CD69 up/down gene heatmap) plots
  the same genes, GPs and values as the published panel; only the two legends'
  vertical sizing differs, from ggplot2 version drift. ~0.047 RMSE, no data
  difference.
- **Figure 6, panel 6i** (GP77) is new, and it is also the panel that exposed a
  bug in `plot_gated_gp_vs_protein()`. The loading-matched set used to be taken as
  `quantile(loadings, 1 - n_prot/N)` with the quantile clamped to 0.9999, so any
  gate below 0.01% of the eligible cells got the top ~0.01% instead of `n_prot`:
  GP77's 18-cell gate was matched against 47 loading-gated cells while the panel
  title said "Matched n=18". The helper now takes exactly the `n_prot`
  highest-loading cells. This is pixel-neutral for the other nine gating panels
  (verified), because their quantile was already exact.
- **Figure S3, panel s3b** (the published s3c) now reports activated CD4/CD8
  cells at threshold 0.1, not all cells at 0.2 as the published panel did. The
  caption states it.
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
- **Figures 3C/3E/3G** and **s3b** had no exact source in the original scripts
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
  apart from the UMAP->MDE axis/legend relabeling, and **4d**'s title was
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
- graph layouts: `set.seed()` before each `ggraph()` call (the barycenter
  TF-GP layout in `code/R/tf_network.R` was deterministic and RNG-free anyway,
  and is unused since Figure 4's TF-GP panel was dropped);
- **label repulsion: every `geom_text_repel()` call passes `seed = 42`.**
  Without it `ggrepel` picks a different layout on every run, which is what
  made panels drift between the PDFs and the site PNGs.

`figures/Previous/bits/` was produced before the repel seeds existed, so label
*placement* (never label content or data) still differs slightly from the
published panels in the label-heavy scatterplots (3A, 3D, 3F, 3H, 4a, 5a, 5d,
6b, 6c, S1C). `4b`'s and `s3g`'s node placement likewise differs from the
published layout because of `graphlayouts`/`igraph` version drift, not RNG.

### Site PNGs

`analysis/assets/**.png` are **converted from the panel PDFs** by
`script/render_site_assets.sh`, not re-plotted. Re-plotting was how the site
and the published panels silently diverged. Run it after any script re-run and
before `workflowr::wflow_build()`.

The single exception is **Figure 7/7B**: `Figure7b_plot.py` writes the published
PDF and `analysis/assets/Figure7/7B.png` from the same matplotlib figure in one
call, so there is no second plotting path that could drift, and re-deriving the
PNG from the PDF at 72 dpi would only make the preview coarser.
`render_site_assets.sh` therefore skips all of `Figure 7/`.

### Panels rasterized by `scattermore`

The ten protein-gate vs. GP-loading panels (**6e-6j** and **s6c-s6f**) draw a
rasterized density layer, so their PDFs differ from the published ones by 0.5-5%
in size and ~0.04 RMSE with no visible difference. Their GP-to-letter assignment
lives in one named vector per script (`fig6_gating` in `Figure6.R`,
`figs6_gating` in `FigureS6.R`) and the loops iterate over its *names*, so a GP
cannot be drawn under another GP's letter. Earlier versions kept the GP list and
the letters in two separate vectors and assigned the letters positionally, which
silently permuted three published panels -- do not reintroduce that shape.

## Note on re-running FigureS6.R (and FigureS5.R)

`cairo_pdf()` was found during verification to not reliably truncate an
existing output file of a different size in place (a stale/partial file
from an earlier interrupted run could persist and get silently reused). That
applied to the two gallery pages, which were retired on 2026-07-28; the current
six panels are written by `ggsave()`, and the script's closing
mtime/size assertion (see "A re-run can silently not write" above) turns any
non-write into a hard error. `FigureS5.R`, split out of it on 2026-07-30, writes
its one panel with `pdf()` and carries the same assertion. If you ever see a script's output that doesn't match
expectations, deleting its output files before re-running is good practice in
general.

## Why most of Figure 7, and all of Figure S7, are excluded

Confirmed via PDF `/Creator` metadata during this refactor:

- Figure 7, panel A is a **Figma** schematic (not code-generated).
- Figure 7, panels C/D/F and **all 11** Figure S7 panels are **Matplotlib**
  (Python) output — but no `.py` file exists anywhere in this repository.
- Figure 7's R-generated panels 7E and 7G trace to `replicate_RQVI_cells.R` and
  related RQVI validation code, but that entire validation track was scoped out
  of this refactor (see `code/README.md`).

These are not reproducible from this repository as it stands. If the
missing Python source turns up, it can be added as a `code/python/`
sibling; until then, those panels have no script here.

**Panel 7B is the exception, and is ours as of 2026-07-28.** The published 7B (a
cell-level EBMF-vs-RQVI correlation heatmap from the collaborator's pipeline) was
dropped and replaced with what had been Extended Data Figure 5, assembled from its
two half-panels into a single PDF. It is built by `Figure7b.R` ->
`Figure7b_rematch.py` -> `Figure7b_plot.py`, documented on
`analysis/Figure7.Rmd`, and has no published counterpart to diff against.
