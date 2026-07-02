library(ggplot2)
library(fastTopics)
library(dplyr)
library(ggrepel)
# data_path <- "../data/"
# code_path <- "../code/"
# data_path <- "data/"
# code_path <- "code/"
data_path <- "/project2/mstephens/immgent/"
code_path <- ""

source(paste0(code_path, "ROC.R"))

num_cells_lower <- 1000

###############
seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
cells_flashier <- rownames(flashier_snmf_summary$L_pm)
cells_seurat <- seurat_meta$cellID
which(!cells_flashier %in% cells_seurat)
L_pm <- flashier_snmf_summary$L_pm[cells_flashier %in% cells_seurat, ]
all(rownames(L_pm) == cells_seurat)
selected_condition <- table(seurat_meta$condition_detailed)
selected_condition <- names(selected_condition[selected_condition >= num_cells_lower])
selected_level2 <- table(seurat_meta$annotation_level2)
selected_level2 <- names(selected_level2[selected_level2 >= num_cells_lower])
selected_cells <- which(seurat_meta$condition_detailed %in% selected_condition & seurat_meta$annotation_level2 %in% selected_level2)
seurat_meta_selected <- seurat_meta[selected_cells, ]
# redefine the levels
seurat_meta_selected$condition_detailed <- factor(seurat_meta_selected$condition_detailed)
seurat_meta_selected$annotation_level2 <- factor(seurat_meta_selected$annotation_level2)
compute_coef_per_factor <- function(factor_id){
  Y <- L_pm[selected_cells, factor_id]
  condition <- seurat_meta_selected$condition_detailed
  level2 <- seurat_meta_selected$annotation_level2
  model_full <- lm(Y ~ condition + level2)
  coefs <- coef(model_full)
  # Extract coefficients for condition
  condition_coefs <- coefs[grep("^condition", names(coefs))]
  # Extract coefficients for annotation_level2
  level2_coefs <- coefs[grep("^level2", names(coefs))]

  return(
    list(
      factor_id = factor_id,
      condition_coefs = condition_coefs,
      level2_coefs = level2_coefs
    )
  )
}
factor_coef_list <- lapply(1:ncol(L_pm), compute_coef_per_factor)
saveRDS(factor_coef_list, file = paste0(data_path, "conditional_analysis_factor_condition_detailed.rds"))



################

seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))

cells_flashier <- rownames(flashier_snmf_summary$L_pm)
cells_seurat <- seurat_meta$cellID
which(!cells_flashier %in% cells_seurat)
L_pm <- flashier_snmf_summary$L_pm[cells_flashier %in% cells_seurat, ]
all(rownames(L_pm) == cells_seurat)

healthy_cells <- which(seurat_meta$condition_broad == "healthy")
L_pm <- L_pm[healthy_cells, ]
seurat_meta <- seurat_meta[healthy_cells, ]

selected_organ <- table(seurat_meta$organ_simplified)
selected_organ <- names(selected_organ[selected_organ >= num_cells_lower])
selected_level2 <- table(seurat_meta$annotation_level2)
selected_level2 <- names(selected_level2[selected_level2 >= num_cells_lower])
selected_cells <- which(seurat_meta$organ_simplified %in% selected_organ & seurat_meta$annotation_level2 %in% selected_level2)
seurat_meta_selected <- seurat_meta[selected_cells, ]
# redefine the levels
seurat_meta_selected$organ_simplified <- factor(seurat_meta_selected$organ_simplified)
seurat_meta_selected$annotation_level2 <- factor(seurat_meta_selected$annotation_level2)

compute_coef_per_factor <- function(factor_id){
  Y <- L_pm[selected_cells, factor_id]
  organ <- seurat_meta_selected$organ_simplified
  level2 <- seurat_meta_selected$annotation_level2
  model_full <- lm(Y ~ organ + level2)
  coefs <- coef(model_full)
  # Extract coefficients for organ
  organ_coefs <- coefs[grep("^organ", names(coefs))]
  # Extract coefficients for annotation_level2
  level2_coefs <- coefs[grep("^level2", names(coefs))]

  return(
    list(
      factor_id = factor_id,
      organ_coefs = condition_coefs,
      level2_coefs = level2_coefs
    )
  )
}

factor_coef_list <- lapply(1:ncol(L_pm), compute_coef_per_factor)
saveRDS(factor_coef_list, file = paste0(data_path, "conditional_analysis_factor_organ_simplified.rds"))
