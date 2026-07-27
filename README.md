# immgenT-GP-analysis

Gene program (GP) factorization analysis of the ImmGen-T single-cell
RNA-seq/CITE-seq data: ~200 GPs learned by empirical Bayes matrix
factorization (EBMF), and how they relate to cell lineage, tissue,
activation state, transcription factors, and surface protein expression.

## Layout

- `data/` -- input data (symlinked from `../immgen-t-factors/data`; not
  tracked in git, see `.gitignore`). Also holds outputs of `code/pipeline/`
  scripts (e.g. the AUC files below) -- not pushed to GitHub either way,
  since it's all local-only data, not part of what this repo shows readers.
- `code/` -- shared R helper functions and the data-prep pipeline,
  sourced by the figure scripts. See `code/README.md`.
- `script/` -- one R script per published figure
  (`Rscript script/FigureN.R`), writing panels into `figures/final-selected/`.
  See `script/README.md` for the figure -> script map and per-panel
  provenance notes.
- `figures/final-selected/` -- the current selected panel set, one PDF per
  panel: everything `script/` regenerates, plus the few panels not produced
  by R in this repo (Figure 7, 1B, 6a -- see `script/README.md`).
- `figures/Previous/` -- the earlier published panels, kept as the
  byte-comparison ground truth.
- `analysis/` -- the workflowr site source (one page per figure: code,
  image, caption). Build with `workflowr::wflow_build()`; output goes to
  `docs/`.
