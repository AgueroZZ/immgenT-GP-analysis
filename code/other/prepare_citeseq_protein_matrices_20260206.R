#!/usr/bin/env Rscript

# Prepare reproducible CITE-seq protein matrices from the 20260206 ADT-only
# Seurat object. This is a cluster-scale provenance script.
#
# Usage:
#   Rscript code/other/prepare_citeseq_protein_matrices_20260206.R \
#     /project2/mstephens/immgent
#
# Outputs in data_path:
#   seurat_meta_20260206.rds              current cell metadata
#   protein_mat.rds                       raw ADT counts (cells x proteins)
#   protein_mat_normalized.rds            CLR-normalized ADT (legacy name)
#   protein_mat_normalized_lognorm.rds    LogNormalize ADT used by protein EBMF

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
data_path <- if (length(args) >= 1L) args[[1]] else "/project2/mstephens/immgent"
seurat_file <- file.path(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
)

message("Reading the 20260206 ADT-only Seurat object: ", seurat_file)
seurat_obj <- readRDS(seurat_file)
seurat_meta <- seurat_obj@meta.data

required_metadata <- c("cellID", "cite_seq", "nCount_ADT")
missing_metadata <- setdiff(required_metadata, colnames(seurat_meta))
if (length(missing_metadata) > 0L) {
  stop(
    "Required metadata columns are missing: ",
    paste(missing_metadata, collapse = ", ")
  )
}
if (!"ADT" %in% names(seurat_obj@assays)) {
  stop("The input Seurat object does not contain an ADT assay.")
}
if (!identical(rownames(seurat_meta), as.character(seurat_meta$cellID))) {
  stop("Metadata row names and cellID are not identical and ordered.")
}

metadata_file <- file.path(data_path, "seurat_meta_20260206.rds")
saveRDS(seurat_meta, metadata_file)
message("Saved versioned metadata: ", metadata_file)

expected_cell_ids <- rownames(seurat_meta)
expected_proteins <- rownames(seurat_obj[["ADT"]])
expected_dim <- c(length(expected_cell_ids), length(expected_proteins))

validate_protein_matrix <- function(x, label) {
  if (!identical(dim(x), expected_dim)) {
    stop(label, " has unexpected dimensions: ", paste(dim(x), collapse = " x "))
  }
  if (!identical(rownames(x), expected_cell_ids)) {
    stop(label, " cell IDs do not match the 20260206 metadata.")
  }
  if (!identical(colnames(x), expected_proteins)) {
    stop(label, " protein names do not match the ADT assay.")
  }
  invisible(TRUE)
}

# Raw ADT counts. Matrices are transposed so rows are cells and columns are
# proteins, matching the downstream flashier scripts.
protein_mat <- t(seurat_obj[["ADT"]]$counts)
validate_protein_matrix(protein_mat, "Raw ADT matrix")
saveRDS(protein_mat, file.path(data_path, "protein_mat.rds"))
rm(protein_mat)
invisible(gc())

# CLR normalization is retained because it was part of the original extraction
# workflow. The legacy filename protein_mat_normalized.rds refers specifically
# to this CLR matrix; it is not the matrix used for the fixed-loading EBMF fit.
message("Computing CLR-normalized ADT matrix (margin = 2)...")
seurat_obj <- NormalizeData(
  seurat_obj,
  assay = "ADT",
  normalization.method = "CLR",
  margin = 2,
  verbose = TRUE
)
protein_mat_normalized_clr <- t(seurat_obj[["ADT"]]$data)
validate_protein_matrix(protein_mat_normalized_clr, "CLR-normalized ADT matrix")
saveRDS(
  protein_mat_normalized_clr,
  file.path(data_path, "protein_mat_normalized.rds")
)
rm(protein_mat_normalized_clr)
invisible(gc())

# The cached protein_mat_normalized_lognorm.rds was audited against the raw
# counts and is exactly reproduced by Seurat LogNormalize with this rounded
# dataset-wide mean ADT library size. For the 20260206 object this evaluates to
# 3472. The mean includes every cell, including cells without CITE-seq counts,
# matching the original cached matrix.
adt_scale_factor <- round(mean(seurat_meta$nCount_ADT, na.rm = TRUE))
message("Computing LogNormalize ADT matrix with scale.factor = ", adt_scale_factor)
if (adt_scale_factor != 3472) {
  warning(
    "The audited 20260206 scale factor is 3472, but this object produced ",
    adt_scale_factor,
    ". Check that the intended object version was supplied."
  )
}
seurat_obj <- NormalizeData(
  seurat_obj,
  assay = "ADT",
  normalization.method = "LogNormalize",
  scale.factor = adt_scale_factor,
  verbose = TRUE
)
protein_mat_normalized_lognorm <- t(seurat_obj[["ADT"]]$data)
validate_protein_matrix(
  protein_mat_normalized_lognorm,
  "LogNormalize ADT matrix"
)
saveRDS(
  protein_mat_normalized_lognorm,
  file.path(data_path, "protein_mat_normalized_lognorm.rds")
)

message(
  "Completed CITE-seq matrix preparation: ",
  nrow(protein_mat_normalized_lognorm),
  " cells x ",
  ncol(protein_mat_normalized_lognorm),
  " proteins; ",
  sum(seurat_meta$cite_seq, na.rm = TRUE),
  " cells have cite_seq == TRUE."
)
