# Extended Data Table 7: protein gating.
#
# One row per GP: its curated positive/negative protein-marker signature, and
# how well a protein gate built from that signature recovers cells with high GP
# loading. For each GP we gate cells by its markers (positive markers above
# their threshold, negative markers at/below), then measure the proportion of
# gated cells that are "positive" for the GP (loading > 0.1; the True Discovery
# Proportion) versus that proportion among all cells, and flag GPs whose gate is
# well aligned with the GP loading.
#
# Reproduces the published data/CITEseq_alignment_scores_manual.csv (Table S3):
# compute_alignment_scores() + format_scores_table() (code/R/gated_protein_helpers.R)
# run on the manually-curated marker set df_markers2 from
# code/R/citeseq_shared_setup.R. Here we keep the presentation subset of columns
# used in Table S3.

library(dplyr)
library(Matrix) # protein matrices are dgCMatrix; must be attached for `[` to dispatch

data_path <- "data/"
output_path <- "figures/generated/"

source("code/R/gated_protein_helpers.R") # compute_alignment_scores(), format_scores_table()
source("code/R/citeseq_shared_setup.R")  # df_markers2, L_pm_filtered, protein_mat_normalized_lognorm,
                                         # threshold_results_subset_manual, select_proteins, *_cells

scores_manual <- compute_alignment_scores(
  df_m = df_markers2,
  protein_mat = protein_mat_normalized_lognorm,
  loading_mat = L_pm_filtered,
  threshold_df = threshold_results_subset_manual,
  selected_proteins = select_proteins,
  exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells),
  missing_threshold_action = "skip"
)
scores_manual_table <- format_scores_table(scores_manual, df_markers2)

# "Well aligned" in the published Table S3 is the formula-positive set
# (prop_pos_gated >= 0.25 & > prop_pos_all; 35 GPs here, matching
# data/CITEseq_alignment_scores_manual.csv) hand-trimmed to those whose protein
# gate and GP-loading gate also agree visually on the MDE embedding -- the 25
# GPs below. We reproduce that curated flag (intersected with the formula set so
# it stays internally consistent if the underlying numbers ever shift). Note
# this 25-GP table set is deliberately distinct from Figure 6b's separately
# curated 27-GP well_aligned_gps in code/R/citeseq_shared_setup.R.
well_aligned_final <- c(
  "GP3", "GP8", "GP10", "GP12", "GP22", "GP25", "GP26", "GP27", "GP29", "GP30",
  "GP32", "GP35", "GP41", "GP58", "GP63", "GP68", "GP77", "GP80", "GP107",
  "GP117", "GP163", "GP166", "GP170", "GP171", "GP181"
)
scores_manual_table[["Well aligned"]] <-
  scores_manual_table[["Well aligned"]] & scores_manual_table$GP %in% well_aligned_final

# Presentation subset used in Table S3, with the two proportion columns given
# their published, more descriptive headers.
protein_gating <- scores_manual_table[, c(
  "GP",
  "Positive markers",
  "Negative markers",
  "Prop. positive (gated)",
  "Prop. positive (all)",
  "Well aligned"
)]
colnames(protein_gating) <- c(
  "GP",
  "Positive markers",
  "Negative markers",
  "Prop. positive among gated (True Discovery Proportion)",
  "Prop. positive among all cells",
  "Well aligned"
)

write.csv(
  protein_gating,
  file = paste0(output_path, "ExtendedDataTable7_protein_gating.csv"),
  row.names = FALSE
)
