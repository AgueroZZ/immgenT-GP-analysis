# Pipeline step 3: per-protein positivity thresholds for CITE-seq gating.
#
# 1. GMM-based threshold per protein (2-component Gaussian mixture on values
#    > 0.5; threshold = upper bound of the "negative" component).
# 2. Manually-reviewed override thresholds for a curated subset of proteins
#    (used preferentially over the GMM threshold in the figure scripts).
# Produces GMM_Thresholds_Summary.csv and Thresholds_Selected_Proteins.csv,
# consumed by Figure6.R, FigureS6.R, and Protein_to_GP.R-derived analyses.
#
# The per-protein diagnostic histogram/scatter PNG galleries from the
# original scripts (one file per protein, purely for visual QC) are dropped
# here -- they aren't required to reproduce any figure panel.
#
# Source: merged from script/Protein_Binarization.R (GMM thresholding) and
# script/protein_thresholding_manual.R (manual overrides). The manual-override
# script referenced `target_proteins` before it was ever assigned in that
# script (it only worked because Protein_Binarization.R happened to leave
# the same variable name in the global environment from a prior run in the
# same session) -- fixed here by computing it once, before first use.

library(dplyr)
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
target_proteins <- setdiff(all_proteins, "CD62L")

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
# 2. Manual overrides for a curated subset of proteins
# ============================================================
threshold_results_subset <- threshold_results[threshold_results$Protein %in% target_proteins, ]
rownames(threshold_results_subset) <- seq_len(nrow(threshold_results_subset))

threshold_results_subset_manual <- c(
  "B220" = 2.75, "BTLA.CD272" = 3.25, "CD2" = 5.75, "CD4" = 3.85, "CD5" = 3.3,
  "CD11A" = 5, "CD24" = 4, "CD27" = 4, "CD29" = 3.5, "CD31" = 4.5, "CD38" = 4,
  "CD39" = 3, "CD44" = 5, "CD49B" = 3, "CD103" = 3.35, "CD155.PVR" = 3.25,
  "CD160" = 3, "CD55.DAF" = 4.25, "CD73.5NTD" = 5, "CD80" = 3, "CD86" = 4,
  "CD8A" = 3.95, "CD8B" = 4, "FR4" = 5, "GITR.CD357" = 3,
  "GR1-LY6G-LY6C1-LY6C2" = 2.5, "ICAM1" = 3, "ICOS.CD278" = 3,
  "IL7RA.CD127" = 4, "ITA4.CD49D" = 3.95, "ITAM.CD11B" = 2.25,
  "ITAX.CD11C" = 2, "KLRG1" = 2.85, "LY49A" = 3, "LY108" = 3, "CD45RB" = 7,
  "NEUROPILIN1.CD304" = 3, "SCA1" = 5, "SLAM.CD150" = 3.75, "TCRGD" = 3.1,
  "TCRVG2" = 2.3, "TCRVG3" = 3.25, "THY1.2" = 7, "IL2RA.CD25" = 2.15,
  "ITB7" = 5, "CD69" = 4
)
threshold_results_subset$Threshold_manual <- threshold_results_subset_manual[match(threshold_results_subset$Protein, names(threshold_results_subset_manual))]

write.csv(threshold_results_subset, paste0(data_path, "Thresholds_Selected_Proteins.csv"), row.names = FALSE)
