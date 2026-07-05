# Pipeline step 2: per-GP predictive AUC.
#
# Computes one-vs-rest AUC (+ optimal threshold) of each GP's loading for
# predicting major lineage (level_1), sub-lineage (level_2), and organ.
# Two different cell restrictions are used by different downstream
# consumers, and both live in data/ under unambiguous names:
#
#   level_{1,2}_AUC_list_figure_no_thymocytes_healthy.rds
#   organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds
#     -- non-thymocyte, HEALTHY-ONLY cells (condition_broad == "healthy").
#     Feeds Figure2.R (2A), Figure4.R (organ mapping in healthy tissue), and
#     the Extended Data GP-AUC tables (lineage / tissue / cluster) plus the
#     Extended Data GP-summary annotations (script/ExtendedDataTable1-4_*.R).
#
#   level_{1,2}_AUC_list_figure_no_thymocytes.rds
#   organ_simplified_AUC_list_figure_no_thymocytes.rds
#     -- non-thymocyte cells, healthy AND diseased together (no healthy
#     restriction). Fed the retired TableS1.R (now removed): confirmed against
#     the published Table S1.xlsx, whose Organ column includes disease-specific
#     sites (SLO, prostate, pancreas, synovial fluid) that only exist in this
#     broader population (e.g. GP4/GP35 -> SLO, GP61 -> pancreas, several GPs ->
#     synovial fluid). The Extended Data tables that replaced Table S1 use the
#     healthy-only variant above instead, so this broader family is currently
#     unused by any table but kept in data/ for reference.
#
# These _no_thymocytes.rds files are identical to the old, ambiguously-named
# level_1/2/organ_simplified_AUC_list_figure.rds (no suffix) -- confirmed via
# direct recomputation earlier (organ_simplified: byte-identical; level_1/
# level_2: diff ~1e-12, floating-point noise) -- just copied under a clearer
# name rather than recomputed. The old no-suffix files are left in data/
# untouched (redundant, but harmless) and are no longer read by anything in
# script/ or code/.
#
# NOTE: the old Table S1.xlsx also had a Condition column (per-GP AUC against
# condition_detailed, non-thymocyte cells, no healthy restriction). The Extended
# Data GP-summary that replaced Table S1 intentionally has no Condition column,
# so it is not reproduced below. (For reference: the cached
# data/condition_detailed_AUC_list_figure.rds predates the current Seurat object
# and was never reverified; a fresh run of 127 condition_detailed categories x
# 200 GPs takes ~3 hours, timed at ~4.5 min for 5 GPs.)
#
# Source: ported from runAUC.R, unchanged apart from path variables and
# output filenames (see above).

library(dplyr)
source("code/R/roc_auc.R")

data_path <- "data/"

seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]

# Restrict to non-thymocyte cells (shared by both variants below)
cells_thymocyte <- which(seurat_meta_filtered$annotation_level1 == "thymocyte")
L_pm_no_thymocytes <- L_pm_filtered[-cells_thymocyte, ]
seurat_meta_no_thymocytes <- seurat_meta_filtered[-cells_thymocyte, ]

compute_and_save <- function(loading_mat, group_info, out_file) {
  auc_list <- compute_auc_threshold_matrix(loading_mat = loading_mat, group_info = group_info)
  rownames(auc_list$auc) <- sort(unique(group_info))
  rownames(auc_list$threshold) <- sort(unique(group_info))
  saveRDS(auc_list, file = paste0(data_path, out_file))
}

## -- Variant 1: non-thymocyte, HEALTHY-ONLY (Figure2 2A, Figure4) --
cells_healthy <- which(seurat_meta_no_thymocytes$condition_broad == "healthy")
L_pm_no_thymocytes_healthy <- L_pm_no_thymocytes[cells_healthy, ]
seurat_meta_no_thymocytes_healthy <- seurat_meta_no_thymocytes[cells_healthy, ]

compute_and_save(L_pm_no_thymocytes_healthy, seurat_meta_no_thymocytes_healthy$annotation_level1, "level_1_AUC_list_figure_no_thymocytes_healthy.rds")
compute_and_save(L_pm_no_thymocytes_healthy, seurat_meta_no_thymocytes_healthy$annotation_level2, "level_2_AUC_list_figure_no_thymocytes_healthy.rds")
compute_and_save(L_pm_no_thymocytes_healthy, seurat_meta_no_thymocytes_healthy$organ_simplified, "organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds")

## -- Variant 2: non-thymocyte, ALL cells healthy + diseased (retired TableS1) --
## Already computed and cached under the old ambiguous no-suffix names;
## copied here rather than rerun (verified equivalent, see header note). No
## current table consumes these -- kept for reference only.
# compute_and_save(L_pm_no_thymocytes, seurat_meta_no_thymocytes$annotation_level1, "level_1_AUC_list_figure_no_thymocytes.rds")
# compute_and_save(L_pm_no_thymocytes, seurat_meta_no_thymocytes$annotation_level2, "level_2_AUC_list_figure_no_thymocytes.rds")
# compute_and_save(L_pm_no_thymocytes, seurat_meta_no_thymocytes$organ_simplified, "organ_simplified_AUC_list_figure_no_thymocytes.rds")

## -- condition_detailed (old TableS1 Condition column, not reproduced) --
# compute_and_save(L_pm_no_thymocytes, seurat_meta_no_thymocytes$condition_detailed, "condition_detailed_AUC_list_figure_no_thymocytes.rds")
