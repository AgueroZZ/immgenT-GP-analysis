# Pipeline step 4: project CITE-seq protein measurements onto the fixed
# scRNA-derived cell loadings (L), re-estimating a protein factor matrix U
# via EBMF (Figure 6's caption panel (a) schematic: Y ~= L U^T).
#
# GAP: this produces protein_flash_selected_summary_lognorm.rds, but the
# file actually consumed by Figure6.R/FigureS6.R is
# protein_flash_selected_summary_lognorm_backfit200.rds (a longer-backfit
# variant, presumably maxiter=200 instead of the 100+40 used here). No
# script anywhere in this repository produces that exact backfit200 file --
# it was most likely produced by re-running this same script with a larger
# `maxiter` interactively and saving under a different name. Treat
# ..._backfit200.rds as an existing input until that variant is
# reconstructed.
#
# Source: ported from code/analyze_protein_selected_lognorm_full.R,
# unchanged apart from path variables and this header.

library(flashier)
library(fastTopics)

data_path <- "data/" # original script used a cluster path; adjust as needed

seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
protein_mat_normalized <- readRDS(paste0(data_path, "protein_mat_normalized_lognorm.rds"))
protein_mat_normalized <- as.matrix(protein_mat_normalized)
scRNA_result <- readRDS(paste0(data_path, "flashier_snmf_summary.rds")) # see 01_extract_data.R gap note

cells_measured <- seurat_meta$cellID[seurat_meta$cite_seq]
L_mat <- scRNA_result$L_pm[rownames(protein_mat_normalized), , drop = FALSE]
L_mat <- L_mat[cells_measured, , drop = FALSE]
protein_mat_normalized <- protein_mat_normalized[cells_measured, , drop = FALSE]

# OLS projection, used as the flashier initialization
qrL <- qr(L_mat)
U_t <- qr.coef(qrL, protein_mat_normalized)
U <- t(U_t)
saveRDS(U, file = paste0(data_path, "protein_projection_OLS_lognorm.rds"))

F_mat_int <- U
flash_fixed_loading <- flash_init(protein_mat_normalized, var_type = 2) |>
  flash_set_verbose(1) |>
  flash_factors_init(list(L_mat, F_mat_int), ebnm_fn = ebnm_point_laplace) |>
  flash_factors_fix(kset = 1:ncol(L_mat), which_dim = "loadings")

flash_fixed_loading <- flash_backfit(flash_fixed_loading, extrapolate = FALSE, maxiter = 100, verbose = 2)
flash_fixed_loading <- flash_backfit(flash_fixed_loading, extrapolate = TRUE, maxiter = 40, verbose = 2)
saveRDS(flash_fixed_loading, file = paste0(data_path, "protein_flash_fixed_loading_lognorm.rds"))

protein_flash_summary_lognorm <- list(
  L_pm = flash_fixed_loading$L_pm,
  F_pm = flash_fixed_loading$F_pm,
  pve = flash_fixed_loading$pve,
  elbo = flash_fixed_loading$elbo,
  residuals_sd = flash_fixed_loading$residuals_sd
)
saveRDS(protein_flash_summary_lognorm, file = paste0(data_path, "protein_flash_selected_summary_lognorm.rds"))
