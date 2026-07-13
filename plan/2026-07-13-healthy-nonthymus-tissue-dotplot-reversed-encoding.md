# Reversed encoding candidate for the filtered tissue dotplot

## Goal

Create a directly comparable tissue dotplot in which dot area shows the within-GP normalized mean loading and color shows the raw mean loading.

## Plan

- [x] Generalize the tissue-dotplot renderer so the color and area metrics can be selected explicitly.
- [x] Render the reversed-encoding candidate using the existing raw-mean threshold and dominant-group order.
- [x] Verify the output dimensions, filtering, and visual legibility.
- [x] Record the new candidate in the experiment README and update log.
