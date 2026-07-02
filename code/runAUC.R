### Re-run AUC analysis
library(ggplot2)
library(fastTopics)
library(dplyr)
library(ggrepel)
library(tidyr)
library(tibble)
library(cowplot)
library(patchwork)

# data_path <- "/project2/mstephens/immgent/"
data_path <- "data/"
code_path <- "code/"
source(paste0(code_path, "ROC.R"))
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered),]

# For the AUC analysis, let's do not consider thymocytes
cells_thymocyte <- which(seurat_meta_filtered$annotation_level1 == "thymocyte")
L_pm_no_thymocytes <- L_pm_filtered[-cells_thymocyte, ]
seurat_meta_no_thymocytes <- seurat_meta_filtered[-cells_thymocyte, ]

# For the AUC analysis, let's focus on healthy samples only
cells_healthy <- which(seurat_meta_no_thymocytes$condition_broad == "healthy")
L_pm_no_thymocytes_healthy <- L_pm_no_thymocytes[cells_healthy, ]
seurat_meta_no_thymocytes_healthy <- seurat_meta_no_thymocytes[cells_healthy, ]





# Level 1:
level_1_AUC_list <- compute_auc_threshold_matrix(loading_mat = L_pm_no_thymocytes_healthy, group_info = seurat_meta_no_thymocytes_healthy$annotation_level1)
rownames(level_1_AUC_list$auc) <- sort(unique(seurat_meta_no_thymocytes_healthy$annotation_level1))
rownames(level_1_AUC_list$threshold) <- sort(unique(seurat_meta_no_thymocytes_healthy$annotation_level1))
saveRDS(level_1_AUC_list, file = paste0(data_path, "level_1_AUC_list_figure.rds"))

# Level 2:
level_2_AUC_list <- compute_auc_threshold_matrix(loading_mat = L_pm_no_thymocytes_healthy, group_info = seurat_meta_no_thymocytes_healthy$annotation_level2)
rownames(level_2_AUC_list$auc) <- sort(unique(seurat_meta_no_thymocytes_healthy$annotation_level2))
rownames(level_2_AUC_list$threshold) <- sort(unique(seurat_meta_no_thymocytes_healthy$annotation_level2))
saveRDS(level_2_AUC_list, file = paste0(data_path, "level_2_AUC_list_figure.rds"))

# organ:
organ_AUC_list <- compute_auc_threshold_matrix(loading_mat = L_pm_no_thymocytes_healthy, group_info = seurat_meta_no_thymocytes_healthy$organ_simplified)
rownames(organ_AUC_list$auc) <- sort(unique(seurat_meta_no_thymocytes_healthy$organ_simplified))
rownames(organ_AUC_list$threshold) <- sort(unique(seurat_meta_no_thymocytes_healthy$organ_simplified))
saveRDS(organ_AUC_list, file = paste0(data_path, "organ_simplified_AUC_list_figure.rds"))








