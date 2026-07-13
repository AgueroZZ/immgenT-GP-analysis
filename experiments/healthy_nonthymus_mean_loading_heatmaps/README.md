# Healthy non-thymus GP mean-loading heatmaps

This isolated experiment draws large heatmaps for all 200 gene programs (GPs)
using the same healthy non-thymocyte reference population as the GP AUC
pipeline: `condition_broad == "healthy"` and
`annotation_level1 != "thymocyte"`.

It does not modify formal figures, analysis pages, or the workflowr site.

## Matrices

For each group, every cell's fitted GP loading is averaged. The script renders
two views of each grouping:

- **Raw mean loading:** the direct per-group mean of each GP's fitted loading.
  The raw tissue and level-2 figures use one shared color maximum, so their
  color intensity is comparable.
- **Within-GP normalized mean loading:** each GP row is divided by its largest
  group mean. Every row therefore ranges from zero to one and emphasizes the
  relative tissue or cluster preference of that GP; it does not preserve
  differences in absolute loading magnitude between GPs.

Both views include an Euclidean-distance, complete-linkage hierarchical
clustering version. They also include a triangular-first candidate, with a
single fixed order shared between the raw and normalized version of each
grouping.

The four PDF outputs are:

- `organ_simplified_raw_mean_loading_heatmap.pdf`
- `organ_simplified_row_normalized_mean_loading_heatmap.pdf`
- `annotation_level2_raw_mean_loading_heatmap.pdf`
- `annotation_level2_row_normalized_mean_loading_heatmap.pdf`

The four triangular-first candidates are:

- `organ_simplified_raw_mean_loading_triangular_heatmap.pdf`
- `organ_simplified_row_normalized_mean_loading_triangular_heatmap.pdf`
- `annotation_level2_raw_mean_loading_triangular_heatmap.pdf`
- `annotation_level2_row_normalized_mean_loading_triangular_heatmap.pdf`

For the triangular order, a normalized mean loading of at least `0.50` is
considered visible. Groups are ordered from the most to fewest visible GPs.
GPs are then ordered by the position of their rightmost visible group, then by
their visible-group count and a right-weighted support score. This prioritizes
a step-wise right boundary over strict monotonicity of the total support count.
`*_triangular_order.csv` records the exact order and support diagnostics.

The four dominant-group candidates are:

- `organ_simplified_raw_mean_loading_dominant_group_heatmap.pdf`
- `organ_simplified_row_normalized_mean_loading_dominant_group_heatmap.pdf`
- `annotation_level2_raw_mean_loading_dominant_group_heatmap.pdf`
- `annotation_level2_row_normalized_mean_loading_dominant_group_heatmap.pdf`

For this Figure 4c-style order, every GP is assigned to the group with its
largest raw mean loading. Groups are ordered by their number of dominant GPs.
GPs form contiguous dominant-group blocks and are sorted within each block by
the continuous dominance gap, largest group mean minus second-largest group
mean. No loading threshold is used. The raw and normalized views of a grouping
share one dominant-group order before the final raw-heatmap filtering.
`*_dominant_group_order.csv` records the full assignment and group counts.

The two final raw dominant-group PDFs use a readability filter: retain a GP
only if at least one raw mean loading is at least `0.1`, then retain a group
only if at least one retained GP reaches that cutoff. The matching final
normalized dominant-group PDFs use exactly the same filtered GPs, groups, and
dominant-group orders. Normalized values are calculated from the full matrix
before subsetting, so every retained GP has a raw mean maximum of at least
`0.1`. All full matrices and the hierarchical and triangular candidates remain
unfiltered. The two `*_raw_mean_ge_0.1_dominant_group_filter_summary.csv`
files record the retained dimensions, the matching order CSVs record the
filtered orders, and the two filtered normalized matrix CSVs record the values
shown in the final normalized heatmaps.

The colored strip above each matrix uses the canonical
`ZemmourLib::immgent_colors` palette for the corresponding biological label.

## Run

From the `immgenT-GP-analysis` repository root:

```bash
Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R
```

The renderer also writes the four displayed matrices as CSV files, the number
of reference cells per group, the hierarchical-clustering orders, the two
triangular-order tables, and the two dominant-group order tables.
`heatmap_summary.txt` records the reference-population size, observed group
counts, color scale, ordering methods, and PDF dimensions.
