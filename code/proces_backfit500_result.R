library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(fastTopics)
library(dplyr)
library(pheatmap)
library(scales)
library(scattermore)
library(qs)

data_path <- "../data/"
code_path <- "../code/"
data_path <- "data/"
code_path <- "code/"
source(paste0(code_path, "ROC.R"))
source(paste0(code_path, "filtering_membership.R"))

seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
# flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
flashier_snmf_summary <- qs::qread(paste0(data_path, "fit_K200_longer.qs"))
cells_flashier <- rownames(flashier_snmf_summary$L_pm)
cells_seurat <- seurat_meta$cellID
L_pm <- flashier_snmf_summary$L_pm[cells_flashier %in% cells_seurat, ]
F_pm <- flashier_snmf_summary$F_pm
colnames(L_pm) <- paste0("K", seq_len(ncol(L_pm)))
colnames(F_pm) <- paste0("F", seq_len(ncol(F_pm)))
D <- diag(1 / apply(L_pm, 2, function(x) max(x)))
L <- L_pm %*% D
cells <- filter_cells_by_total_membership(L,numiter = 12)
seurat_meta_filtered <- seurat_meta[cells,]
L_pm_filtered <- L_pm[cells,]
d <- apply(L_pm_filtered,2,max)
L_pm_filtered <- scale_cols(L_pm_filtered,1/d)
F_pm_filtered <- scale_cols(F_pm,d)
L_pm_filtered2 <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
saveRDS(L_pm_filtered, paste0(data_path, "L_pm_filtered_500.rds"))
saveRDS(F_pm_filtered, paste0(data_path, "F_pm_filtered_500.rds"))

