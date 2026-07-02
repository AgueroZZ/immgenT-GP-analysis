library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(fastTopics)
library(dplyr)
library(pheatmap)
library(scales)
library(scattermore)

#####################################################
#####################################################
####################################################
##### Defining directory and loading functions
####################################################
#####################################################
# data_path <- "../data/"
# code_path <- "../code/"
data_path <- "data/"
code_path <- "code/"
output_path <- "outputs/"

#####################################################
#####################################################
#####################################################
### Loading Data
#####################################################
#####################################################
load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
# flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered),]
level_1_AUC_list <- readRDS(file = paste0(data_path, "level_1_AUC_list_figure.rds"))
level_2_AUC_list <- readRDS(file = paste0(data_path, "level_2_AUC_list_figure.rds"))
condition_detailed_AUC_list_figure <- readRDS(file = paste0(data_path, "condition_detailed_AUC_list_figure.rds"))
organ_simplified_AUC_list_figure <- readRDS(file = paste0(data_path, "organ_simplified_AUC_list_figure.rds"))

# Normalize F_pm_filtered such that each column has max abs value of 1
F_pm_filtered <- apply(F_pm_filtered, 2, function(x) x / max(abs(x)))
colnames(F_pm_filtered) <- colnames(L_pm_filtered)

#####################################################
### Produce a table that shows the top 10 genes (by absolute value) for each GP in F_pm_filtered
### (only consider genes with absolute value > 0.25)
### Also, in each column, produce the conditions that pass auc > 0.8 for each GP
#####################################################
top_genes_list <- list()
for (gp in colnames(F_pm_filtered)) {
  top_genes <- rownames(F_pm_filtered)[which(abs(F_pm_filtered[, gp]) > 0.25)]
  top_genes <- top_genes[order(abs(F_pm_filtered[top_genes, gp]), decreasing = TRUE)]
  top_genes_list[[gp]] <- head(top_genes, 10)
}

#####################################################
### Build supplementary table: one row per GP
### Columns: GP, Level1, Level2, Condition, Organ, Signature_Genes
#####################################################

gps <- colnames(F_pm_filtered)

get_passing_categories <- function(auc_list, gps, L_mat, threshold = 0.8) {
  auc_mat <- auc_list$auc
  thr_mat <- auc_list$threshold
  # Median loading per GP across all cells
  gp_median_loading <- apply(L_mat, 2, median)
  # rows = categories, cols = GPs
  lapply(gps, function(gp) {
    if (!gp %in% colnames(auc_mat)) return("")
    auc_vals <- auc_mat[, gp]
    thr_vals <- thr_mat[, gp]
    med_load <- gp_median_loading[gp]
    # Keep only categories predicted by HIGH loading (threshold >= median loading)
    cats <- rownames(auc_mat)[auc_vals > threshold & thr_vals >= med_load]
    paste(cats, collapse = "; ")
  })
}

level1_cats    <- get_passing_categories(level_1_AUC_list,                  gps, L_pm_filtered)
level2_cats    <- get_passing_categories(level_2_AUC_list,                  gps, L_pm_filtered)
condition_cats <- get_passing_categories(condition_detailed_AUC_list_figure, gps, L_pm_filtered)
organ_cats     <- get_passing_categories(organ_simplified_AUC_list_figure,   gps, L_pm_filtered)

sig_genes_pos <- lapply(gps, function(gp) {
  vals <- F_pm_filtered[, gp]
  candidates <- names(vals)[vals > 0.25]
  top <- head(candidates[order(vals[candidates], decreasing = TRUE)], 5)
  paste(top, collapse = "; ")
})

sig_genes_neg <- lapply(gps, function(gp) {
  vals <- F_pm_filtered[, gp]
  candidates <- names(vals)[vals < -0.25]
  top <- head(candidates[order(abs(vals[candidates]), decreasing = TRUE)], 5)
  paste(top, collapse = "; ")
})

supp_table <- data.frame(
  GP                   = gps,
  Level1               = unlist(level1_cats),
  Level2               = unlist(level2_cats),
  Condition            = unlist(condition_cats),
  Organ                = unlist(organ_cats),
  Signature_Genes_Pos  = unlist(sig_genes_pos),
  Signature_Genes_Neg  = unlist(sig_genes_neg),
  stringsAsFactors = FALSE
)

write.csv(supp_table,
          file = paste0(output_path, "Supplementary_Table1_GP_summary.csv"),
          row.names = FALSE)

