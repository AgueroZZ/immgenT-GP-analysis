# immgenT-GP-analysis

Gene program (GP) factorization analysis of the ImmGen-T single-cell
RNA-seq/CITE-seq data: ~200 GPs learned by empirical Bayes matrix
factorization (EBMF), and how they relate to cell lineage, tissue,
activation state, transcription factors, and surface protein expression.

## Layout

- `data/` -- input data (symlinked from `../immgen-t-factors/data`; not
  tracked in git, see `.gitignore`).
- `data-local/` -- small corrected/recomputed data files that supersede a
  stale cached file in `data/`. Kept separate (and git-tracked) because
  `data/` is a symlink into the shared `immgen-t-factors` project, so this
  repo never writes there.
- `code/` -- shared R helper functions and the data-prep pipeline,
  sourced by the figure scripts. See `code/README.md`.
- `script/` -- one R script per published figure
  (`Rscript script/FigureN.R`), writing panels into `figures/generated/`.
  See `script/README.md` for the figure -> script map and per-panel
  provenance notes.
- `figures/generated/` -- regenerated figure panels, one PDF per panel,
  named to match `figures/final-selected/`.
- `figures/final-selected/` -- the published ground-truth panels, kept for
  comparison.
- `analysis/` -- the workflowr site source (one page per figure: code,
  image, caption). Build with `workflowr::wflow_build()`; output goes to
  `docs/`.
