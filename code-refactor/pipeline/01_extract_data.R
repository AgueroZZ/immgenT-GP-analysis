# Pipeline step 1: raw data extraction.
#
# From the raw Seurat object, produces:
#   - shifted_log_counts.qs / counts.qs   (RNA counts, gene-filtered)
#   - protein_mat.rds / protein_mat_normalized.rds  (CITE-seq ADT)
#   - flashier_snmf_summary.rds  (L_pm/F_pm/elbo/pve; condensed from the
#     upstream flashier fit)
#
# GAP: this script reads an existing `flashier_snmf.rds` (line "flashier_snmf
# <- readRDS(...)") -- the flashier semi-NMF fit itself. No script anywhere
# in this repository (old code/ or otherwise) produces that file; it was
# fit interactively/on a cluster and the fitting script was not preserved
# here. The repo owner has this script and will add it separately -- until
# then, treat flashier_snmf.rds as a required upstream input to this step,
# not something this pipeline reproduces from scratch.
#
# This is a cluster-scale job (needs Seurat/flashier/qs and the full raw
# Seurat object, ~683k cells) -- ported here for provenance/documentation,
# not intended to be run against the data/ directory in this repo (which
# only contains the already-extracted, filtered outputs).
#
# Source: ported from code/extract_data.R, unchanged apart from path
# variables and this header.

libs <- c("fastTopics", "flashier", "Matrix", "Seurat", "qs")
invisible(sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))))

data_path <- "/project2/mstephens/immgent/" # cluster path; adjust to your raw-data location

seurat_obj <- readRDS(paste0(data_path, "igt1_96_withtotalvi20250710_clean.Rds"))
seurat_meta <- seurat_obj@meta.data
saveRDS(seurat_meta, file = paste0(data_path, "seurat_meta.rds"))

norm_fac <- mean(seurat_meta$nCount_RNA)
seurat_obj <- NormalizeData(seurat_obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = norm_fac)
counts <- seurat_obj[["RNA"]]$counts
shifted_log_counts <- seurat_obj[["RNA"]]$data

# transpose so rows = cells, cols = genes
counts <- t(counts)
shifted_log_counts <- t(shifted_log_counts)

message("Removing TCR, mt, ribo, Gm, and Rik genes and genes that are not expressed")
tcr_genes <- grepl(x = rownames(seurat_obj), pattern = "Trbv|Trbd|Trbj|Trbc|Trav|Traj|Trac|Trgv|Trgd|Trgj|Trgc|Trdv|Trdj|Trdc")
gm_rik_genes <- grepl(x = rownames(seurat_obj), pattern = "Gm|Rik$|\\-ps$")
ribo_genes <- grepl(x = rownames(seurat_obj), pattern = "Rpl|Rps|Mrpl|Mrps|Rsl")
mt_genes <- grepl(x = rownames(seurat_obj), pattern = "^mt-")
genes_not_expressed <- rowSums(seurat_obj[["RNA"]]$counts) == 0
genes_to_keep <- !(tcr_genes | gm_rik_genes | ribo_genes | mt_genes | genes_not_expressed)
cat("Number of genes to keep:", sum(genes_to_keep), "\n")
shifted_log_counts <- shifted_log_counts[, genes_to_keep]
counts <- counts[, genes_to_keep]
qsave(shifted_log_counts, file = file.path(data_path, "shifted_log_counts.qs"), preset = "balanced")
qsave(counts, file = file.path(data_path, "counts.qs"), preset = "balanced")

# CITE-seq (ADT) protein matrix
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "CLR", margin = 2, assay = "ADT")
protein_mat <- seurat_obj[["ADT"]]$counts
protein_mat_normalized <- seurat_obj[["ADT"]]$data
protein_mat <- t(protein_mat)
protein_mat_normalized <- t(protein_mat_normalized)
saveRDS(protein_mat, file = paste0(data_path, "protein_mat.rds"))
saveRDS(protein_mat_normalized, file = paste0(data_path, "protein_mat_normalized.rds"))

# Condense the full flashier fit (see GAP note above) into the summary
# object every figure script reads.
flashier_snmf <- readRDS(paste0(data_path, "flashier_snmf.rds"))
flashier_snmf_summary <- list(
  L_pm = flashier_snmf$L_pm,
  F_pm = flashier_snmf$F_pm,
  elbo = flashier_snmf$elbo,
  residuals_sd = flashier_snmf$residuals_sd,
  pve = flashier_snmf$pve
)
saveRDS(flashier_snmf_summary, file = paste0(data_path, "flashier_snmf_summary.rds"))
