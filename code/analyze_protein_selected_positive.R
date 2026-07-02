library(flashier)
library(fastTopics)
# library(ggplot2)
# library(ggrepel)
# library(cowplot)
# data_path <- "data/"
data_path <- "/project2/mstephens/immgent/"

# Load data matrix and loading matrix
protein_mat_normalized <- readRDS(file = paste0(data_path, "protein_mat_normalized.rds"))
scRNA_result <- readRDS(file = paste0(data_path, "flashier_snmf_summary.rds"))

# factors to be considered
k <- ncol(scRNA_result$L_pm)
factors_considered <- which(scRNA_result$pve > 1e-4)
L_mat <- scRNA_result$L_pm[rownames(protein_mat_normalized), factors_considered, drop = FALSE]
F_mat_int <- matrix(0, nrow = ncol(protein_mat_normalized), ncol = length(factors_considered))

# remove rows that are entirely zero
zero_rows <- which(rowSums(protein_mat_normalized) == 0)
if(length(zero_rows) > 0){
  protein_mat_normalized <- protein_mat_normalized[-zero_rows, , drop = FALSE]
  L_mat <- L_mat[-zero_rows, , drop = FALSE]
}

# remove columns that are entirely zero
zero_cols <- which(colSums(protein_mat_normalized) == 0)
if(length(zero_cols) > 0){
  protein_mat_normalized <- protein_mat_normalized[, -zero_cols, drop = FALSE]
  F_mat_int <- F_mat_int[-zero_cols, , drop = FALSE]
}

# run flashier fixed loading
flash_fixed_loading <- flash_init(protein_mat_normalized, var_type = 2) |>
  flash_set_verbose(1) |>
  flash_factors_init(list(L_mat,
                          F_mat_int), ebnm_fn = ebnm_point_exponential) |>
  flash_factors_fix(kset = 1:length(factors_considered),
                    which_dim = "loadings") |>
  flash_backfit()

protein_flash_summary_pos <- list(L_pm = flash_fixed_loading$L_pm,
                              F_pm = flash_fixed_loading$F_pm,
                              pve = flash_fixed_loading$pve,
                              elbo = flash_fixed_loading$elbo,
                              residuals_sd = flash_fixed_loading$residuals_sd)

saveRDS(protein_flash_summary_pos, file = paste0(data_path, "protein_flash_selected_summary_pos.rds"))



