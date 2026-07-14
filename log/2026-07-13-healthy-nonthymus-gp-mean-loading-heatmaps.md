# Healthy Non-thymus GP Mean-loading Heatmaps

## Change

Added the isolated experiment directory
`experiments/healthy_nonthymus_mean_loading_heatmaps/`. It renders all 200 GP
cell-loading columns against two healthy non-thymocyte groupings:

- tissue: `organ_simplified`;
- cluster: `annotation_level2`.

The reference population is exactly the one used by the healthy AUC workflow:
`condition_broad == "healthy"` and `annotation_level1 != "thymocyte"`.

Each grouping has two independently clustered PDF heatmaps:

- direct raw mean loading;
- within-GP normalized mean loading, calculated as each row divided by its
  largest group mean.

Both dimensions use Euclidean-distance, complete-linkage hierarchical
clustering. The group-color strip uses the canonical `ZemmourLib` palette.

## Outputs

The renderer writes four PDFs, four displayed matrices, two group-size tables,
and the row/column clustering order for every heatmap. The raw tissue and
level-2 heatmaps share a single upper color limit, enabling direct comparison
of absolute mean-loading intensity.

## Validation

- `Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R` completed successfully.
- The reference contains 337,060 cells, 200 GPs, 18 observed
  `organ_simplified` groups, and 107 observed `annotation_level2` groups.
- All raw matrices have 200 GP rows and finite values. Every normalized GP row
  has maximum value 1 within numerical tolerance.
- All four PDFs are nonempty. The large tissue PDFs are 12.2 x 28.3 inches;
  the large level-2 PDFs are 22.8 x 28.3 inches.
- Current PDF layouts were rendered with Quick Look and visually checked for
  dendrogram, label, color-strip, legend, and cell-grid visibility.

## Triangular-first follow-up

Added four triangular-first candidate PDFs without replacing the original
hierarchical-clustering PDFs. The fixed order is derived separately for tissue
and level-2 matrices from the normalized means, using the Figure 6b-style
support mask `normalized mean >= 0.50`.

- Columns are ordered from most to fewest supported GPs.
- GP rows are ordered by their rightmost supported column, then by support
  count and a right-weighted support score.
- Each grouping's raw and normalized candidate shares exactly one order.

The two `*_triangular_order.csv` files record the order and diagnostics. Both
the tissue and level-2 row boundaries have zero reverse increases, confirming
the intended monotone step-wise envelope. Updated candidate PDFs were rendered
with Quick Look and visually inspected after shortening the tissue title to
avoid clipping.

## Figure 4c-style dominant-group follow-up

Added four threshold-free dominant-group candidates. Each GP is assigned to
the tissue or level-2 cluster with its largest raw mean loading. Groups are
ordered by their number of dominant GPs; GP rows are then arranged into
contiguous dominant-group blocks and, within each block, by decreasing
dominance gap (largest minus second-largest group mean). The raw and normalized
versions of one grouping share this fixed order.

For tissue, 17 of 18 groups are dominant for at least one GP and the largest
block contains 40 GPs. For level2, 69 of 107 groups are represented and the
largest block contains 13 GPs. Visual comparison of the tissue normalized
figures showed that support-triangular ordering emphasizes the cutoff-defined
envelope, whereas dominant-group ordering produces the Figure 4c-like diagonal
block structure without thresholding any loading. Both candidate families are
retained for comparison.

## Tissue normalized-color/raw-area dotplot candidate

Added `organ_simplified_dominant_group_normalized_color_raw_area_dotplot.pdf`.
It uses the tissue dominant-group order and no tile background. Dot color is
the within-GP normalized mean loading, while dot area is the unnormalized raw
mean loading on a single tissue-wide scale with maximum 0.44278271. Raw means
are deliberately not rescaled per GP.

Visual inspection confirmed the intended behavior: an entry can remain red
when it is the relative maximum for a GP, but it stays visually negligible when
its raw mean is trivial. Large red dots require both relative enrichment and a
substantial absolute average loading. The resulting plot is intentionally
sparser than the tile heatmap because it exposes the large dynamic range of
raw mean loading.

## Raw mean >= 0.1 filtered tissue dotplot

Added `organ_simplified_dominant_group_normalized_color_raw_area_raw_mean_ge_0.1_dotplot.pdf`.
The plot retains a GP only when at least one tissue has raw mean loading at
least 0.1, then retains tissues only when at least one retained GP meets that
same threshold. The dominant-group order is recomputed on this filtered matrix,
not inherited from the full 200-GP matrix.

The filter retained 31 of 200 GPs. All 18 tissues remained because each has at
least one retained GP with raw mean loading at least 0.1. The filtered PDF uses
an adaptive 10-inch height rather than the 28.3-inch height needed for the
full 200-GP plot. Visual inspection confirmed readable labels, compact rows,
and the intended dual encoding: normalized color plus unscaled raw dot area.

## Filtered tissue dotplot with reversed encoding

Added `organ_simplified_dominant_group_raw_color_normalized_area_raw_mean_ge_0.1_dotplot.pdf`.
It uses exactly the same raw-mean cutoff (0.1), retained 31 GPs and 18 tissues,
and recomputed dominant-group order as the existing filtered dotplot. The only
change is the visual mapping: dot color is the raw mean loading, and dot area
is the within-GP normalized mean loading.

The renderer now exposes the color and area metrics explicitly, while retaining
the previous normalized-color/raw-area mapping as its default. The new PDF is
nonempty (38,908 bytes), and a Quick Look rendering was visually checked for
the color gradient, normalized-area legend, ordering, and label legibility.

## Final selection: dominant-group heatmaps only

Removed all tissue dotplot code, documentation, PDFs, and the dotplot-only
filter summary after selecting the heatmaps as the final display. The selected
dominant-group ordering remains unchanged, and the experiment continues to
render the four requested dominant-group heatmaps:

- tissue (`organ_simplified`): raw and within-GP normalized mean loading;
- level2 (`annotation_level2`): raw and within-GP normalized mean loading.

The renderer completed successfully after cleanup. All four final PDFs are
nonempty, and fresh Quick Look checks of the normalized tissue and level2
heatmaps confirmed the dominant-block structure, labels, legends, and color
strips render correctly. `rg -n -i 'dotplot' experiments/healthy_nonthymus_mean_loading_heatmaps`
returns no matches.

## Filtered final raw dominant-group heatmaps

The final raw dominant-group heatmaps now retain only a GP row with at least
one raw mean loading >= 0.1 and only a group column with at least one retained
GP meeting that cutoff. Dominant-group ordering is recomputed after filtering.
The normalized dominant-group heatmaps, full matrices, and hierarchical and
triangular candidates remain unchanged.

The tissue raw heatmap retained 31 of 200 GPs and all 18 tissues. The level2
raw heatmap retained 64 of 200 GPs and all 107 level2 clusters; every retained
column has at least one retained GP at or above the cutoff. Both PDF files are
nonempty, and fresh Quick Look renderings confirmed legible titles, labels,
legends, group-color strips, and the recomputed dominant-group structure.

## Matched final normalized dominant-group heatmaps

Updated the final normalized dominant-group heatmaps to use the exact same
raw-filtered GP rows, group columns, and dominant-group order as the matching
raw heatmaps. Within-GP normalization is still calculated from the full raw
matrix before filtering; therefore every GP shown has a full-matrix raw maximum
of at least 0.1, avoiding relative patterns based on near-zero maxima.

The final tissue pair is 31 GPs by 18 tissues, and the level2 pair is 64 GPs
by 107 clusters. The two filtered normalized matrices are saved as CSV files,
and both normalized PDFs were regenerated and visually checked with Quick Look
for matching dominant-block structure, labels, legends, and color strips.

## Formal Figure S4

Added a standalone formal renderer at `script/FigureS4.R` and its workflowr
page at `analysis/FigureS4.Rmd`. Figure S4a is the tissue
(`organ_simplified`) heatmap and Figure S4b is the level2
(`annotation_level2`) heatmap. Each panel exposes two clickable alternatives:
raw mean loading and within-GP normalized mean loading.

For a given panel, the two alternatives share exactly the same raw-filtered
GPs, groups, canonical color strip, and dominant-group order. The filter is
raw mean loading >= 0.1, applied before order calculation; normalization is
calculated on the full matrix before the matching subset is taken. S4a retains
31 GPs and 18 tissues, while S4b retains 64 GPs and 107 clusters.

`Rscript script/FigureS4.R` generated four one-page PDFs and
`S4_filter_summary.csv` under `figures/generated/Figure S4/`. Their PNG
renderings were synchronized to both `analysis/assets/FigureS4/` and
`docs/assets/FigureS4/`. `workflowr::wflow_build()` then built
`docs/FigureS4.html` and the updated index. The HTML contains both tabsets,
all source/docs image pairs are byte-identical, and raw/normalized tissue and
level2 alternatives were visually checked. No `figures/final-selected/` files
were changed because the raw-versus-normalized choice remains open.

## Centered alternatives and Figure 1 level2 ordering

Added a third, row-centered alternative to each Figure S4 panel. Centered
values are calculated from the full mean-loading matrix as each GP's group mean
minus that GP's mean across all groups, before taking the same raw-filtered
subset used by the raw and normalized views. The centered color scale is
symmetric blue-white-red around zero. Figure S4 now contains clickable raw,
normalized, and centered tabs for both tissue and level2 panels.

The level2 columns in all three formal S4b alternatives now follow the Figure
1 level1 sequence `CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP`, with
`annotation_level2` labels alphabetized within each level1 block. GP rows stay
in dominant-level2 blocks relative to this fixed column order. Two annotation
strips identify the level1 block and individual level2 group.

The experiment renderer also writes unfiltered, internal centered PDFs with
all 200 GPs and all observed groups:

- `organ_simplified_row_centered_mean_loading_full_dominant_group_heatmap.pdf`
  (200 GPs by 18 tissues);
- `annotation_level2_row_centered_mean_loading_full_level1_order_heatmap.pdf`
  (200 GPs by 107 level2 groups).

The corresponding centered matrix CSVs and row/column order CSVs are retained
in the experiment directory. The full level2 order was checked against the
healthy non-thymocyte metadata: its distinct level1 blocks match the confirmed
sequence exactly and each block is alphabetically ordered. The maximum absolute
centered-row mean was `2.523942e-16`. Both formal and internal renderers
completed successfully; all six formal PDFs and the four centered experiment
PDFs are one-page files. Quick Look inspection confirmed the centered
diverging scale, labels, and the level1/level2 annotation strips.

## Centered positive-threshold filtering

Changed the centered alternatives so they no longer inherit the raw-mean
filter. The full mean-loading matrix is row-centered first. A GP is then kept
only if at least one positive centered entry is >= 0.1; a group is kept only if
at least one retained GP has a positive centered entry >= 0.1. Centered rows
and columns are ordered independently after this filter, while raw and
normalized alternatives continue to share the original raw-mean-filtered set.

The resulting formal centered dimensions are 22 GPs by 16 tissues for S4a and
61 GPs by 100 level2 clusters for S4b. The centered level2 columns retain the
Figure 1 level1-first, alphabetical-within-level1 sequence among the 100 kept
clusters. The all-200-GP internal centered PDFs remain unfiltered.

Both renderers completed successfully. Filter summary CSVs, filtered centered
matrix CSVs, and centered order CSVs were regenerated in the experiment
directory. Quick Look inspection confirmed the new centered tissue and level2
layouts, diverging color scale, and level1/level2 strips.

## Heatmap typography update

Increased heatmap typography across the formal Figure S4 and the experiment
renderer. Column labels are now vertical, and row/column label font sizes are
calculated from the available heatmap cell dimensions, with a readability-first
range of 9--14 pt. Titles, annotation labels, row titles, and legend text were
also enlarged. This keeps labels as large as each matrix layout can support
without column-label overlap, including the wide level2 heatmaps.

All formal PDFs and the experiment outputs were regenerated. Quick Look
inspection of the tissue and level2 centered panels confirmed clearly larger
labels, vertical non-overlapping column names, and readable legends.

## Current-output cleanup and GP37 diagnosis

Removed obsolete hierarchical, triangular, and unfiltered raw/normalized
artifacts from `experiments/healthy_nonthymus_mean_loading_heatmaps/`. The
renderer now creates only the current six filtered raw/normalized/centered
heatmaps, the two unfiltered all-200-GP centered internal heatmaps, and the
matrices, group counts, filter summaries, and order CSVs needed to interpret
those eight PDFs. A future run no longer recreates the retired candidates.

GP37 was checked because it is mammary-gland-specific but absent from the
filtered tissue panels. Its largest raw tissue mean is `0.01122608` and its
largest centered tissue mean is `0.01057881`, both in mammary gland. These are
below the current raw and positive-centered display thresholds of `0.1`, so its
absence is caused by filtering rather than annotation or ordering.

## Per-GP SD centered filtering and compact figure sizing

Replaced the fixed centered-mean `0.1` filter with a GP-specific variability
rule. For each GP, the renderer calculates the SD of that GP's mean loading
across all groups. A centered entry is supported when it is at least `2 x` that
GP-specific SD. Centered heatmaps retain GP rows and group columns containing
at least one supported entry, then recompute the dominant-group row order.
Raw and normalized views remain unchanged and continue to share the raw-mean
`>= 0.1` subset.

The centered tissue view now retains 182 GPs and 17 tissues; only lung is
removed. The centered level2 view retains all 200 GPs and 103 of 107 clusters;
the removed clusters are `CD4.E`, `DP.wA`, `gdT.K`, and `Tz.C`. The unfiltered
internal centered views remain 200 GPs by 18 tissues and 200 GPs by 107 level2
clusters.

GP37 now passes the centered tissue filter. Its maximum centered tissue mean is
`0.01057881`, while its GP-specific SD is `0.00264075` and its `2 x SD` cutoff
is `0.00528150`. In level2, its maximum centered mean is `0.01434731` and its
`2 x SD` cutoff is `0.00515319`.

Heatmap height is now proportional to the number of displayed GPs, with a
160 mm minimum and 3.5 mm per row, instead of the previous fixed 480 mm
minimum. This reduces the filtered tissue raw/normalized PDF height from the
oversized layout to 9.8 inches while keeping row labels at the largest
non-overlapping size. Centered panels expand vertically as needed for their
larger retained GP sets.

## Fixed 0.01 centered filtering

Replaced the per-GP `2 x SD` centered filter after review showed that it retained
too many GP rows. The centered display now uses a fixed positive cutoff: a GP
row is retained when at least one centered group mean is `>= 0.01`, and a group
column is retained when at least one retained GP meets the same cutoff.

This produces a substantially smaller centered display: 70 GPs by all 18
tissues and 112 GPs by all 107 level2 clusters. GP37 remains in both centered
views. Its maximum centered tissue mean is `0.01057881`, just above the fixed
cutoff, and its maximum centered level2 mean is `0.01434731`.

Raw and normalized views remain unchanged at 31 GPs by 18 tissues and 64 GPs
by 107 level2 clusters using the raw mean `>= 0.1` filter. The internal full
centered views also remain unfiltered at 200 GPs by 18 tissues and 200 GPs by
107 level2 clusters. Experiment artifact names now use `centered_mean_ge_0.01`;
the retired `centered_mean_ge_2sd` matrices, summaries, and order files were
removed.
