# Pipeline step 2: per-GP predictive AUC.
#
# Computes one-vs-rest AUC (+ optimal threshold) of each GP's loading for
# predicting major lineage (level_1), sub-lineage (level_2), and organ,
# restricted to non-thymocyte, healthy cells. Feeds Figure2.R (2A) and
# TableS1.R.
#
# GAP: Figure4.R (script-refactor/) consumes
# level_1_AUC_list_figure_no_thymocytes_healthy.rds /
# organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds -- differently
# named from this script's output despite this script already restricting
# to non-thymocyte/healthy cells. Both variants exist side-by-side in
# data/, suggesting a later rerun changed the output naming without an
# updated script being preserved. This step reproduces the *_figure.rds
# names actually consumed by Figure2.R/TableS1.R; the *_no_thymocytes_healthy
# variants are treated as an existing input for Figure4.R.
#
# Source: ported from code/runAUC.R, unchanged apart from path variables.

library(dplyr)
source("code-refactor/R/roc_auc.R")

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
saveRDS(level_1_AUC_list, file = paste0(data_path, "level_1_AUC_list_figure.rds"))

level_2_AUC_list <- compute_auc_threshold_matrix(loading_mat = L_pm_no_thymocytes_healthy, group_info = seurat_meta_no_thymocytes_healthy$annotation_level2)
rownames(level_2_AUC_list$auc) <- sort(unique(seurat_meta_no_thymocytes_healthy$annotation_level2))
rownames(level_2_AUC_list$threshold) <- sort(unique(seurat_meta_no_thymocytes_healthy$annotation_level2))
saveRDS(level_2_AUC_list, file = paste0(data_path, "level_2_AUC_list_figure.rds"))

organ_AUC_list <- compute_auc_threshold_matrix(loading_mat = L_pm_no_thymocytes_healthy, group_info = seurat_meta_no_thymocytes_healthy$organ_simplified)
rownames(organ_AUC_list$auc) <- sort(unique(seurat_meta_no_thymocytes_healthy$organ_simplified))
rownames(organ_AUC_list$threshold) <- sort(unique(seurat_meta_no_thymocytes_healthy$organ_simplified))
saveRDS(organ_AUC_list, file = paste0(data_path, "organ_simplified_AUC_list_figure.rds"))
