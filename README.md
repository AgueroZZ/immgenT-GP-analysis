# immgenT-GP-analysis

Gene program (GP) factorization analysis of the ImmGen-T single-cell
RNA-seq/CITE-seq data: ~200 GPs learned by empirical Bayes matrix
factorization (EBMF), and how they relate to cell lineage, tissue,
activation state, transcription factors, and surface protein expression.

## Layout

- `data/` -- input data (symlinked from `../immgen-t-factors/data`; not
  tracked in git, see `.gitignore`).
- `code-refactor/` -- shared R helper functions and the data-prep pipeline,
  sourced by the figure scripts. See `code-refactor/README.md`.
- `script-refactor/` -- one R script per published figure
  (`Rscript script-refactor/FigureN.R`), writing panels into
  `figure-refactor/`. See `script-refactor/README.md` for the figure ->
  script map and per-panel provenance notes.
- `figure-refactor/` -- regenerated figure panels, one PDF per panel, named
  to match `figures/final-selected/bits/`.
- `figures/final-selected/` -- the published ground-truth panels, kept for
  comparison.
- `analysis/` -- the workflowr site source (one page per figure: code,
  image, caption). Build with `workflowr::wflow_build()`; output goes to
  `docs/`.
- `code/`, `script/` -- the original, pre-refactor analysis code, kept for
  provenance/history.

