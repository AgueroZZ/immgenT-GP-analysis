# Healthy non-thymus GP mean-loading heatmaps

This experiment uses healthy non-thymocyte cells only:
`condition_broad == "healthy"` and `annotation_level1 != "thymocyte"`.
For each tissue (`organ_simplified`) or level2 cluster (`annotation_level2`),
it averages fitted GP loadings across cells.

## Current outputs

The final filtered heatmaps are retained in raw, normalized, and centered
versions for both groupings:

- `organ_simplified_raw_mean_loading_dominant_group_heatmap.pdf`
- `organ_simplified_row_normalized_mean_loading_dominant_group_heatmap.pdf`
- `organ_simplified_row_centered_mean_loading_dominant_group_heatmap.pdf`
- `annotation_level2_raw_mean_loading_dominant_group_heatmap.pdf`
- `annotation_level2_row_normalized_mean_loading_dominant_group_heatmap.pdf`
- `annotation_level2_row_centered_mean_loading_dominant_group_heatmap.pdf`

Raw and normalized views keep a GP when its raw group mean is at least `0.1`
in one group, then retain groups with at least one retained GP meeting that
threshold. Centered views first subtract each GP's mean across groups, then
independently retain rows and columns containing at least one centered mean of
`0.01` or greater. Raw/normalized and centered views therefore can retain
different matrices and use separate dominant-group orders.

This centered rule retains 70 GPs and all 18 tissues, and 112 GPs with all 107
level2 clusters. The raw/normalized views remain 31 GPs by 18 tissues and 64
GPs by 107 level2 clusters.

For level2 heatmaps, columns follow the Figure 1 level1 order `CD8`, `CD4`,
`Treg`, `gdT`, `CD8aa`, `Tz`, `DN`, and `DP`; labels are alphabetized within
each level1 block. GP rows are grouped by their dominant group and ordered by
decreasing dominance gap within group.

Two unfiltered internal centered PDFs retain all 200 GPs and all observed
groups for collaborator review:

- `organ_simplified_row_centered_mean_loading_full_dominant_group_heatmap.pdf`
- `annotation_level2_row_centered_mean_loading_full_level1_order_heatmap.pdf`

The matching filtered matrices, full centered matrices, group-count tables,
filter summaries, and row/column order CSVs are retained alongside the PDFs.

## GP37 note

GP37 is a mammary-gland-specific program in the tissue matrix, but it does not
meet the raw-display threshold: its maximum raw tissue mean is `0.01123` in
mammary gland. It is retained in the centered tissue view because its maximum
centered tissue mean is `0.01058`, which is just above the centered cutoff of
`0.01`. The fixed centered cutoff therefore captures this specific
low-amplitude program without changing the raw/normalized filter.

## Run

From the `immgenT-GP-analysis` repository root:

```bash
Rscript experiments/healthy_nonthymus_mean_loading_heatmaps/render_mean_loading_heatmaps.R
```
