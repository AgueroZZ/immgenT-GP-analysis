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
# are NOT the same computation, despite using the same (current) annotation:
# the no-suffix files turn out to be missing the healthy-only restriction
# below (i.e. computed on ALL non-thymocyte cells, healthy and diseased
# together) -- confirmed exactly by rerunning with that one restriction
# dropped (organ_simplified: byte-identical, max diff 0; level_1/level_2:
# max diff ~1e-12, floating-point noise). This also explains why the
# no-suffix organ_simplified file has 4 extra categories (SLO, prostate,
# pancreas, synovial fluid) not present in the healthy-only version -- these
# are disease-model-specific sample sites. (An earlier theory blamed a stale
# cached seurat_meta.rds from an older script -- ruled out: the category
# labels in the no-suffix files are the same current-format labels as
# everywhere else, not the old numeric-cluster style that stale file used.)
# Figure2.R and TableS1.R were updated to read the _no_thymocytes_healthy
# files, matching Figure4.R and making all three figures/tables consistent.
# The no-suffix files are superseded and no longer read by anything in
# script/.
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
# value left), so the "missing healthy filter" explanation above doesn't
# apply here -- whatever's wrong with the cached
# data/condition_detailed_AUC_list_figure.rds (dated Dec 2 2025, predating
# even the current Seurat object) was never independently root-caused, since
# dropping the column made it moot. Recovered source logic was an
# `eval=FALSE` chunk in analysis/old/Figures_Manuscript_v1.rmd. This file is
# no longer produced or read by anything in script/.
