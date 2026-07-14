# Figure S4 centered fixed 0.01 filtering plan

Status: completed on 2026-07-14.

## Goal

Replace the centered per-GP SD filter with a fixed positive centered-mean cutoff that produces a more selective display while retaining GP37.

## Tasks

1. Update both renderers to retain a centered entry when `centered mean >= 0.01`.
2. Retain GP rows and group columns containing at least one qualifying centered entry, then recompute the dominant-group order.
3. Keep raw and normalized filtering unchanged at raw mean `>= 0.1`.
4. Keep the internal full centered views unfiltered at all 200 GPs and all observed groups.
5. Rename centered experiment artifacts from the retired `2sd` rule to the fixed `0.01` rule and remove stale artifacts.
6. Regenerate formal PDFs and web previews, rebuild `docs/FigureS4.html`, visually inspect the layouts, and document the result in the analysis, README, and log.

## Validation result

- Centered tissue: 70 GPs by all 18 tissues; GP37 retained.
- Centered level2: 112 GPs by all 107 clusters; GP37 retained.
- Raw/normalized and full internal centered dimensions are unchanged.
- All formal PDFs are one page, the centered layouts were visually inspected, and source/docs PNG pairs are byte-identical.
