#!/usr/bin/env Rscript

# Fit CITE-seq protein scores while fixing the scRNA-derived GP cell loadings.
# The model is initialized by OLS and refined with flashier checkpoints through
# 200 cumulative backfitting iterations.
#
# Usage:
#   Rscript code/other/fit_citeseq_fixed_loading_ebmf_20260206.R \
#     /project2/mstephens/immgent

suppressPackageStartupMessages({
  library(flashier)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
data_path <- if (length(args) >= 1L) args[[1]] else "/project2/mstephens/immgent"

metadata_file <- file.path(data_path, "seurat_meta_20260206.rds")
protein_file <- file.path(data_path, "protein_mat_normalized_lognorm.rds")
scrna_file <- file.path(data_path, "flashier_snmf_summary.rds")

message("Reading 20260206 metadata: ", metadata_file)
seurat_meta <- readRDS(metadata_file)
message("Reading LogNormalize protein matrix: ", protein_file)
protein_mat_normalized <- readRDS(protein_file)
message("Reading scRNA GP fit summary: ", scrna_file)
scRNA_result <- readRDS(scrna_file)

required_metadata <- c("cellID", "cite_seq")
missing_metadata <- setdiff(required_metadata, colnames(seurat_meta))
if (length(missing_metadata) > 0L) {
  stop(
    "Required metadata columns are missing: ",
    paste(missing_metadata, collapse = ", ")
  )
}
if (is.null(scRNA_result$L_pm) || is.null(rownames(scRNA_result$L_pm))) {
  stop("flashier_snmf_summary.rds does not contain a row-named L_pm matrix.")
}
if (is.null(rownames(protein_mat_normalized))) {
  stop("The protein matrix does not have cell IDs as row names.")
}

# Use the 20260206 metadata order and retain only measured cells that occur in
# both matrices. This explicit intersection avoids silent NA rows when object
# versions differ by a small number of cells.
cells_measured <- as.character(
  seurat_meta$cellID[!is.na(seurat_meta$cite_seq) & seurat_meta$cite_seq]
)
cells_used <- cells_measured[
  cells_measured %in% rownames(protein_mat_normalized) &
    cells_measured %in% rownames(scRNA_result$L_pm)
]

if (length(cells_used) == 0L) {
  stop("No measured CITE-seq cells are shared by the metadata and both matrices.")
}
if (length(cells_used) < length(cells_measured)) {
  warning(
    length(cells_measured) - length(cells_used),
    " measured cells were absent from the protein or scRNA loading matrix."
  )
}

L_mat <- scRNA_result$L_pm[cells_used, , drop = FALSE]
protein_mat_normalized <- as.matrix(
  protein_mat_normalized[cells_used, , drop = FALSE]
)
if (!identical(rownames(L_mat), rownames(protein_mat_normalized))) {
  stop("The aligned scRNA loading and protein matrices have different row order.")
}

message(
  "Aligned matrices: ",
  nrow(L_mat),
  " cells x ",
  ncol(L_mat),
  " fixed GP loadings; ",
  ncol(protein_mat_normalized),
  " proteins."
)

# OLS projection gives a protein x GP score matrix used to initialize flashier.
message("Computing the OLS protein-score initialization...")
qrL <- qr(L_mat)
U_t <- qr.coef(qrL, protein_mat_normalized)
U <- t(U_t)
saveRDS(U, file.path(data_path, "protein_projection_OLS_lognorm.rds"))

message("Initializing flashier and fixing the scRNA-derived loading matrix...")
flash_fixed_loading <- flash_init(protein_mat_normalized, var_type = 2) |>
  flash_set_verbose(1) |>
  flash_factors_init(
    list(L_mat, U),
    ebnm_fn = ebnm_point_laplace
  ) |>
  flash_factors_fix(
    kset = seq_len(ncol(L_mat)),
    which_dim = "loadings"
  )

save_checkpoint <- function(fit, cumulative_iterations, output_dir) {
  fit_file <- file.path(
    output_dir,
    paste0(
      "protein_flash_fixed_loading_lognorm_backfit",
      cumulative_iterations,
      ".rds"
    )
  )
  summary_file <- file.path(
    output_dir,
    paste0(
      "protein_flash_selected_summary_lognorm_backfit",
      cumulative_iterations,
      ".rds"
    )
  )

  saveRDS(fit, fit_file)
  protein_flash_summary <- list(
    L_pm = fit$L_pm,
    F_pm = fit$F_pm,
    pve = fit$pve,
    elbo = fit$elbo,
    residuals_sd = fit$residuals_sd
  )
  saveRDS(protein_flash_summary, summary_file)
  message("Saved cumulative backfit checkpoint ", cumulative_iterations)
  invisible(NULL)
}

checkpoint_schedule <- data.frame(
  cumulative_iterations = c(20L, 40L, 80L, 120L, 160L, 200L),
  additional_iterations = c(20L, 20L, 40L, 40L, 40L, 40L),
  extrapolate = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
)

for (i in seq_len(nrow(checkpoint_schedule))) {
  flash_fixed_loading <- flash_backfit(
    flash_fixed_loading,
    extrapolate = checkpoint_schedule$extrapolate[[i]],
    maxiter = checkpoint_schedule$additional_iterations[[i]],
    verbose = 2
  )
  save_checkpoint(
    flash_fixed_loading,
    checkpoint_schedule$cumulative_iterations[[i]],
    data_path
  )
}

message("Completed fixed-loading CITE-seq EBMF through backfit200.")
