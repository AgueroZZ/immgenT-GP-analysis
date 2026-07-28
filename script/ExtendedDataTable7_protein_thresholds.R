# Extended Data Table 7: manual protein positivity thresholds.
#
# One row per surface protein that carries a manually reviewed positivity
# threshold, and the hand-set value used to call a cell positive for that
# protein. A cell counts as positive when its log-normalized ADT value is
# strictly above the threshold (and negative when it is at or below), which is
# exactly how the protein gates in Figure 6, Figure S6, and the CITE-seq
# alignment scores are built.
#
# Provenance of the values: pipeline step 3 (code/pipeline/03_protein_thresholds.R)
# first fits a 2-component Gaussian mixture per protein and takes the upper edge
# of the negative component as an automatic cutoff. Those automatic cutoffs were
# then reviewed by eye against each protein's ADT histogram and replaced by the
# rounded, hand-set values collected here (column Threshold_manual of
# data/Thresholds_Selected_Proteins.csv) -- the values that
# code/R/citeseq_shared_setup.R actually loads and gates on.
#
# Note: this table is the published threshold set as loaded at run time. The
# hard-coded manual vector inside code/pipeline/03_protein_thresholds.R is a
# different curation round (only 12 of its 46 entries match the CSV), so it is
# deliberately not used here.

data_path <- "data/"
output_path <- "figures/final-selected/"

thresholds <- read.csv(
  paste0(data_path, "Thresholds_Selected_Proteins.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)
thresholds <- thresholds[, c("Protein", "Threshold_manual")]

# CD62L is thresholded at 3 in code/R/citeseq_shared_setup.R rather than in the
# CSV (pipeline step 3 leaves it out of the manually-reviewed subset), so it is
# appended here to match the threshold set the gating code actually sees.
thresholds <- rbind(
  thresholds,
  data.frame(Protein = "CD62L", Threshold_manual = 3)
)

# radix ordering so the row order is locale-independent
thresholds <- thresholds[order(thresholds$Protein, method = "radix"), ]
rownames(thresholds) <- NULL
colnames(thresholds) <- c("Protein", "Manual threshold")

write.csv(
  thresholds,
  file = paste0(output_path, "ExtendedDataTable7_protein_thresholds.csv"),
  row.names = FALSE
)
