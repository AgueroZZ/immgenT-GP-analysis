# Figure S1. GP reproducibility across IGTs.

**(a)** GPs were tested for reproducibility in individual datasets (IGTs): for each IGT a dataset-specific EBMF factorization was computed and its factors were matched to the global GPs by Hungarian assignment on the cosine similarity of gene-score vectors (restricted to shared genes, with columns scaled), and a GP is counted as "validated" in an IGT when this cosine similarity reaches a given threshold. Adding IGTs one at a time in index order, the curves show the cumulative number of GPs (of 200) validated in at least one IGT included so far, for cosine thresholds of 0.2–0.8.

**(b)** Using the same per-IGT cosine matching, the number of GPs validated by at least X IGTs as a function of X (both axes log-scaled), for thresholds of 0.2–0.8.

**(c)** Between-IGT variability of GP loadings, computed on a standard spleen subset used for this purpose. Each GP's mean loading is computed within every IGT; each point is a GP, plotting the mean across IGTs of these per-IGT mean loadings (x-axis) against their variance across IGTs (y-axis). The ten GPs with the highest between-IGT variance are labeled.

**(d)** Heatmap of the per-IGT mean loading (same standard spleen subset, restricted to IGTs with ≥ 500 cells) for the ten GPs with the highest between-IGT variance from (c). Rows (GPs) are hierarchically clustered; color runs from white (low) to red (high mean loading).
