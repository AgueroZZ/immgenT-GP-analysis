# Pipeline step 3: per-protein positivity thresholds for CITE-seq gating.
#
# 1. GMM-based threshold per protein (2-component Gaussian mixture on values
#    > 0.5; threshold = upper bound of the "negative" component). Writes
#    GMM_Thresholds_Summary.csv.
# 2. Consistency check against the hand-curated override thresholds in
#    Thresholds_Selected_Proteins.csv. That file is a curated INPUT, not an
#    output of this step -- see below. Its Threshold_manual column is what
#    R/citeseq_shared_setup.R gates on, and therefore what Figure 6,
#    Figure S6, and Extended Data Table 7 use.
#
# The per-protein diagnostic histogram/scatter PNG galleries from the
# original scripts (one file per protein, purely for visual QC) are dropped
# here -- they aren't required to reproduce any figure panel.
#
# Source: merged from script/Protein_Binarization.R (GMM thresholding) and
# script/protein_thresholding_manual.R (manual overrides).
#
# WHY STEP 2 NO LONGER WRITES Thresholds_Selected_Proteins.csv: the original
# protein_thresholding_manual.R generated that file from a hard-coded 46-entry
# vector, and the file was then hand-revised -- a second reviewer's pass, kept
# in the CSV's Threshold_david column with free-text notes, 5 proteins dropped
# (BTLA.CD272, GR1-LY6G-LY6C1-LY6C2, LY108, SLAM.CD150, THY1.2), and
# Threshold_manual updated to match. Only 12 of the 46 original entries survive
# that revision, so re-running the old write would silently roll every protein
# gate back to the superseded draft. The revised CSV is authoritative; this
# step now only verifies it is still in sync with the GMM run above.

library(mclust)

data_path <- "data/"

L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
protein_mat_normalized_lognorm <- readRDS(paste0(data_path, "protein_mat_normalized_lognorm.rds"))
proteins_quality <- read.csv(paste0(data_path, "TableS4_citeseq_qc_20250513.csv"), header = TRUE, stringsAsFactors = FALSE, skip = 1)
cells_citeseq <- seurat_meta_filtered$cellID[seurat_meta_filtered$cite_seq]

protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[intersect(rownames(protein_mat_normalized_lognorm), cells_citeseq), ]
poor_proteins <- proteins_quality$protein[proteins_quality$classification == "poor"]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[, !colnames(protein_mat_normalized_lognorm) %in% poor_proteins]

thymocyte_cells <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 == "thymocyte"]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[!rownames(protein_mat_normalized_lognorm) %in% thymocyte_cells, ]

all_proteins <- colnames(protein_mat_normalized_lognorm)

# ============================================================
# 1. GMM-based threshold per protein
# ============================================================
threshold_results <- data.frame(Protein = character(), Threshold = numeric(), stringsAsFactors = FALSE)
message("Starting GMM-based thresholding for counts > 0.5...")
for (prot in all_proteins) {
  prot_vals <- protein_mat_normalized_lognorm[, prot]
  non_zero_vals <- prot_vals[prot_vals > 0.5]
  if (length(non_zero_vals) < 50) {
    message(sprintf("Skipping %s: too few entries > 0.5.", prot))
    next
  }
  set.seed(42)
  gmm_fit <- Mclust(non_zero_vals, G = 2, verbose = FALSE)
  means <- gmm_fit$parameters$mean
  neg_cluster <- which.min(means)
  threshold <- max(non_zero_vals[gmm_fit$classification == neg_cluster], na.rm = TRUE)
  threshold_results <- rbind(threshold_results, data.frame(Protein = prot, Threshold = threshold))
}
write.csv(threshold_results, paste0(data_path, "GMM_Thresholds_Summary.csv"), row.names = FALSE)
message("GMM Thresholding complete.")

# ============================================================
# 2. Consistency check against the curated manual thresholds
#
# Read-only. Thresholds_Selected_Proteins.csv is hand-curated and must not be
# regenerated here (see the header note). What we can check is that its
# Threshold column -- which was populated from an earlier run of the GMM above
# -- still matches the GMM thresholds just computed. A mismatch means the
# curated file predates a change in the protein matrix or the GMM step, and
# that the manual values were reviewed against different histograms than the
# ones this run would produce.
# ============================================================
curated_path <- paste0(data_path, "Thresholds_Selected_Proteins.csv")
curated <- read.csv(curated_path, header = TRUE, stringsAsFactors = FALSE)

message(sprintf(
  "Curated manual thresholds: %d proteins in %s (column Threshold_manual).",
  sum(!is.na(curated$Threshold_manual)), basename(curated_path)
))

gmm_now <- threshold_results$Threshold[match(curated$Protein, threshold_results$Protein)]
drift <- abs(gmm_now - curated$Threshold) > 1e-6 | is.na(gmm_now)
if (any(drift)) {
  warning(sprintf(
    paste0(
      "%s is out of sync with this GMM run for %d protein(s): %s.\n",
      "  The curated manual thresholds were set against the older GMM values. ",
      "Re-review them before trusting the gates; do NOT overwrite the file."
    ),
    basename(curated_path), sum(drift), paste(curated$Protein[drift], collapse = ", ")
  ))
} else {
  message("Curated file is in sync with this GMM run; manual thresholds unchanged.")
}
