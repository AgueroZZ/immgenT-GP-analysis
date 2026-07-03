# Pipeline step 2: per-GP predictive AUC.
#
# Computes one-vs-rest AUC (+ optimal threshold) of each GP's loading for
# predicting major lineage (level_1), sub-lineage (level_2), and organ,
# restricted to non-thymocyte, healthy cells. Feeds Figure2.R (2A), Figure4.R,
# and TableS1.R.
#
# RESOLVED (previously logged here as a naming-mismatch GAP): two families of
# cached AUC files exist in data/ -- level_1/2/organ_simplified_AUC_list_figure.rds
# (no suffix) and the same names with a _no_thymocytes_healthy suffix. These
# are NOT the same computation: the no-suffix files turned out to be a stale
# artifact of an older, pre-refactor script (a chunk in
# analysis/old/Figures_Manuscript_v1.rmd) that read the now-known-stale cached
# seurat_meta.rds and did not restrict to healthy cells. The
# _no_thymocytes_healthy files match this script's logic (and this script's
# output) essentially exactly (max abs AUC diff ~1e-13). Figure2.R and
# TableS1.R were updated to read the _no_thymocytes_healthy files, matching
# Figure4.R and making all three figures/tables consistent. The no-suffix
# files are superseded and no longer read by anything in script/.
#
# Source: ported from runAUC.R, unchanged apart from path variables and
# output filenames (see above).

library(dplyr)
source("code/R/roc_auc.R")

data_path <- "data/"

seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]

# Restrict to non-thymocyte, healthy cells
cells_thymocyte <- which(seurat_meta_filtered$annotation_level1 == "thymocyte")
L_pm_no_thymocytes <- L_pm_filtered[-cells_thymocyte, ]
seurat_meta_no_thymocytes <- seurat_meta_filtered[-cells_thymocyte, ]

cells_healthy <- which(seurat_meta_no_thymocytes$condition_broad == "healthy")
L_pm_no_thymocytes_healthy <- L_pm_no_thymocytes[cells_healthy, ]
seurat_meta_no_thymocytes_healthy <- seurat_meta_no_thymocytes[cells_healthy, ]

level_1_AUC_list <- compute_auc_threshold_matrix(loading_mat = L_pm_no_thymocytes_healthy, group_info = seurat_meta_no_thymocytes_healthy$annotation_level1)
rownames(level_1_AUC_list$auc) <- sort(unique(seurat_meta_no_thymocytes_healthy$annotation_level1))
rownames(level_1_AUC_list$threshold) <- sort(unique(seurat_meta_no_thymocytes_healthy$annotation_level1))
saveRDS(level_1_AUC_list, file = paste0(data_path, "level_1_AUC_list_figure_no_thymocytes_healthy.rds"))

level_2_AUC_list <- compute_auc_threshold_matrix(loading_mat = L_pm_no_thymocytes_healthy, group_info = seurat_meta_no_thymocytes_healthy$annotation_level2)
rownames(level_2_AUC_list$auc) <- sort(unique(seurat_meta_no_thymocytes_healthy$annotation_level2))
rownames(level_2_AUC_list$threshold) <- sort(unique(seurat_meta_no_thymocytes_healthy$annotation_level2))
saveRDS(level_2_AUC_list, file = paste0(data_path, "level_2_AUC_list_figure_no_thymocytes_healthy.rds"))

organ_AUC_list <- compute_auc_threshold_matrix(loading_mat = L_pm_no_thymocytes_healthy, group_info = seurat_meta_no_thymocytes_healthy$organ_simplified)
rownames(organ_AUC_list$auc) <- sort(unique(seurat_meta_no_thymocytes_healthy$organ_simplified))
rownames(organ_AUC_list$threshold) <- sort(unique(seurat_meta_no_thymocytes_healthy$organ_simplified))
saveRDS(organ_AUC_list, file = paste0(data_path, "organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds"))

# NOTE: an earlier version of this script also computed a
# condition_detailed_AUC_list_figure.rds (grouped by condition_detailed,
# restricted to non-thymocyte cells only -- not additionally healthy-only,
# since condition_detailed is the condition variable itself). That column
# was dropped from TableS1.R: restricting a condition-predicting AUC to
# healthy-only cells is degenerate (condition_broad would have just one
# value left), and the cached data/condition_detailed_AUC_list_figure.rds
# (which did NOT restrict to healthy) turned out to have the same
# stale-seurat_meta.rds problem as level_1/2/organ above -- recovered
# source was an `eval=FALSE` chunk in analysis/old/Figures_Manuscript_v1.rmd.
# Since the column is gone, this file is no longer produced or read by
# anything in script/.
