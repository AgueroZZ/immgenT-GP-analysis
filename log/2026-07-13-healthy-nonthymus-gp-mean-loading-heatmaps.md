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
