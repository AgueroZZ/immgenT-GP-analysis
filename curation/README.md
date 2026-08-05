# Curated inputs

Hand-maintained files that feed the reproducible pipeline. Unlike everything in
`data/` (a symlink to the shared `immgen-t-factors/data`), these are small text
files tracked in this repo so that every edit is reviewable in the git history.

| File | Consumed by | What it is |
| --- | --- | --- |
| `GP_manual_annotations.csv` | `script/ExtendedDataTable1_GP_summary.R` | One row per GP (`GP1`..`GP200`) with a free-text `Comments`. Blank means "no comment". |

## `GP_manual_annotations.csv`

The manual GP comments, taken from the finalized `Comments` column of
`tables/ExtendedDataTable1_target.xlsx` on 2026-08-05 (82 of the 200 GPs carry
one). That spreadsheet is a local hand edit -- `tables/` is gitignored -- so
this CSV, not it, is what the build reads. They started life as the `annotation` column of the retired
hand-maintained `Table S1.xlsx` (the collaborator's spreadsheet, keyed
`K1`..`K200`; `K`*n* is the same factor as `GP`*n*), but every one of those 22
entries has since been rewritten or expanded. **This CSV is the single source of
truth**, and Extended Data Table 1 is where the comments are published.

To add or revise a comment, edit the `Comments` cell for that GP here
(the file is plain CSV -- Excel-friendly, all 200 GPs pre-listed so no GP name
has to be typed), then re-run
`Rscript script/ExtendedDataTable1_GP_summary.R` and rebuild the
`ExtendedDataTable1` page. The script errors out if the `GP` column is not
exactly `GP1`..`GP200`, so a row can't be silently dropped or renamed.
