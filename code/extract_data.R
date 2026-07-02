data_path <- "/project2/mstephens/immgent/"
libs = c("fastTopics", "flashier", "Matrix", "Seurat", "qs") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

# read in data:
seurat_obj <- readRDS(paste0(data_path, "igt1_96_withtotalvi20250710_clean.Rds"))
seurat_meta <- seurat_obj@meta.data
saveRDS(seurat_meta, file = paste0(data_path, "seurat_meta.rds"))
seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))

norm_fac = mean(seurat_meta$nCount_RNA)
seurat_obj = NormalizeData(seurat_obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = norm_fac)
counts = seurat_obj[['RNA']]$counts
shifted_log_counts = seurat_obj[['RNA']]$data

#transpose!
counts = t(counts)
shifted_log_counts = t(shifted_log_counts)
message("Removing TCR, mt, ribo, Gm, and Rik genes and genes that are not expressed")
tcr_genes = grepl(x = rownames(seurat_obj), pattern = "Trbv|Trbd|Trbj|Trbc|Trav|Traj|Trac|Trgv|Trgd|Trgj|Trgc|Trdv|Trdj|Trdc")
gm_rik_genes = grepl(x = rownames(seurat_obj), pattern = "Gm|Rik$|\\-ps$")
ribo_genes = grepl(x = rownames(seurat_obj), pattern = "Rpl|Rps|Mrpl|Mrps|Rsl")
mt_genes = grepl(x = rownames(seurat_obj), pattern = "^mt-")
genes_not_expressed = rowSums(seurat_obj[["RNA"]]$counts) == 0
genes_to_keep = !(tcr_genes | gm_rik_genes | ribo_genes | mt_genes | genes_not_expressed)
cat("Number of genes to keep:", sum(genes_to_keep), "\n")
shifted_log_counts = shifted_log_counts[,genes_to_keep]
counts = counts[,genes_to_keep]
qsave(shifted_log_counts, file = file.path(data_path, "shifted_log_counts.qs"), preset = "balanced")
qsave(counts,             file = file.path(data_path, "counts.qs"),             preset = "balanced")

# cite-seq
seurat_obj = NormalizeData(seurat_obj, normalization.method = "CLR", margin = 2, assay = "ADT")
protein_mat = seurat_obj[['ADT']]$counts
protein_mat_normalized = seurat_obj[['ADT']]$data
protein_mat <- t(protein_mat); protein_mat_normalized <- t(protein_mat_normalized)
saveRDS(protein_mat, file = paste0(data_path, "protein_mat.rds"))
saveRDS(protein_mat_normalized, file = paste0(data_path, "protein_mat_normalized.rds"))

# full matrix factorization result
flashier_snmf <- readRDS(paste0(data_path, "flashier_snmf.rds"))
flashier_snmf_summary <- list(L_pm = flashier_snmf$L_pm,
                              F_pm = flashier_snmf$F_pm,
                              elbo = flashier_snmf$elbo,
                              residuals_sd = flashier_snmf$residuals_sd,
                              pve = flashier_snmf$pve)
saveRDS(flashier_snmf_summary, file = paste0(data_path, "flashier_snmf_summary.rds"))
