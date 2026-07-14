# Figure S4 centered per-GP SD filtering plan

Status: completed on 2026-07-14.

## Goal

Update the centered Figure S4 views so that specificity is evaluated relative to each GP's own across-group variability, while retaining the existing raw/normalized filtering and the unfiltered 200-GP internal centered views.

## Tasks

1. Replace the fixed centered-mean cutoff with a per-GP rule.
   - Compute `sd(mean_loading[GP, all groups])` separately for each GP.
   - Mark an entry as supported when `centered_mean >= 2 * GP_SD`.
   - Retain rows and columns containing at least one supported entry.
   - Validate the rule with row-specific cutoffs and confirm that GP37 is retained.

2. Preserve all other Figure S4 definitions.
   - Keep raw and normalized tabs on the shared raw-mean `>= 0.1` row/column set.
   - Keep tissue dominant-group ordering.
   - Keep level2 columns in Figure 1 level1 order and alphabetical level2 order within level1.
   - Keep the internal full centered views unfiltered at all 200 GPs.

3. Improve figure proportions.
   - Size heatmap height from the number of GP rows rather than imposing the previous 480 mm minimum.
   - Preserve readable, vertical group labels and the largest non-overlapping row-label font.

4. Regenerate and document outputs.
   - Rebuild formal PDFs, PNG previews, and `docs/FigureS4.html`.
   - Rebuild the cleaned experiment directory and rename centered-filter artifacts to describe the `2 SD` rule.
   - Update the analysis narrative, experiment README, output summaries, and project log.
   - Visually inspect representative tissue and level2 PDFs and verify dimensions, retained sets, and GP37 values.

## Validation result

- Centered tissue: 182 GPs by 17 tissues; GP37 retained.
- Centered level2: 200 GPs by 103 clusters; GP37 retained.
- Internal full centered views remain 200 GPs by 18 tissues and 200 GPs by 107 level2 clusters.
- All six formal PDFs are one page, source/docs PNG pairs are byte-identical, and the rebuilt HTML reports the per-GP `2 x SD` rule.
