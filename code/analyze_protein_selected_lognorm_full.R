library(flashier)
library(fastTopics)
# data_path <- "data/"
data_path <- "/project2/mstephens/immgent/"

# Load data matrix and loading matrix
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
protein_mat_normalized <- readRDS(file = paste0(data_path, "protein_mat_normalized_lognorm.rds"))
protein_mat_normalized <- as.matrix(protein_mat_normalized)
scRNA_result <- readRDS(file = paste0(data_path, "flashier_snmf_summary.rds"))
# cells with CITE-seq measurements
cells_measured <- seurat_meta$cellID[seurat_meta$cite_seq]
# Assign L_mat based on the cells that are measured in the protein matrix
L_mat <- scRNA_result$L_pm[rownames(protein_mat_normalized), , drop = FALSE]
# cells_unmeasured <- seurat_meta$cellID[!seurat_meta$cite_seq]
# cells_unmeasured_in_protein_mat <- cells_unmeasured[cells_unmeasured %in% rownames(protein_mat_normalized)]

# remove cells that are not measured in the protein matrix
L_mat <- L_mat[cells_measured, , drop = FALSE]
protein_mat_normalized <- protein_mat_normalized[cells_measured, , drop = FALSE]


# projection (OLS)
qrL <- qr(L_mat)                         # QR of N x K
U_t <- qr.coef(qrL, protein_mat_normalized)               # K x Q
U <- t(U_t)                          # Q x K
saveRDS(U, file = paste0(data_path, "protein_projection_OLS_lognorm.rds"))

# use OLS projection as initialization for flashier
F_mat_int <- U
# run flashier fixed loading
flash_fixed_loading <- flash_init(protein_mat_normalized, var_type = 2) |>
  flash_set_verbose(1) |>
  flash_factors_init(list(L_mat,
                          F_mat_int), ebnm_fn = ebnm_point_laplace) |>
  flash_factors_fix(kset = 1:ncol(L_mat),
                    which_dim = "loadings")

flash_fixed_loading <- flash_backfit(flash_fixed_loading,extrapolate = FALSE,maxiter = 100,verbose = 2)
flash_fixed_loading <- flash_backfit(flash_fixed_loading,extrapolate = TRUE, maxiter = 40,verbose = 2)
saveRDS(flash_fixed_loading, file = paste0(data_path, "protein_flash_fixed_loading_lognorm.rds"))

protein_flash_summary_lognorm <- list(L_pm = flash_fixed_loading$L_pm,
                                      F_pm = flash_fixed_loading$F_pm,
                                      pve = flash_fixed_loading$pve,
                                      elbo = flash_fixed_loading$elbo,
                                      residuals_sd = flash_fixed_loading$residuals_sd)

saveRDS(protein_flash_summary_lognorm, file = paste0(data_path, "protein_flash_selected_summary_lognorm.rds"))



