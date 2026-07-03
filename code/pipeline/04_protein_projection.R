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
# Source: ported from analyze_protein_selected_lognorm_full.R,
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

# ============================================================
# CITEseq_markers_full.rds: per-GP positive/negative protein markers
# (|score| >= 0.5 on the max|score|=1-scaled protein factor matrix).
#
# Source: recovered from live (not commented-out) code in the original
# Figure_CITEseq.R -- the marker-selection logic itself was never removed,
# just fed by a *_backfit200.rds variant of protein_flash_summary_lognorm
# that this pipeline doesn't reproduce exactly (see GAP note above). Using
# this step's own protein_flash_summary_lognorm instead, so the result may
# differ slightly from the cached data/CITEseq_markers_full.rds.
# ============================================================
Protein_F_pm <- protein_flash_summary_lognorm$F_pm
colnames(Protein_F_pm) <- paste0("GP", seq_len(ncol(Protein_F_pm)))

isotype_proteins <- grep("^Isotype", rownames(Protein_F_pm), value = TRUE)
proteins_quality <- read.csv(paste0(data_path, "TableS4_citeseq_qc_20250513.csv"), header = TRUE, stringsAsFactors = FALSE, skip = 1)
good_proteins <- c(proteins_quality$protein[proteins_quality$classification == "good"], "IL2RA.CD25", "ITB7", "CD69")
exclude_proteins <- c("CD19", "CD34", "CD45.1", "CD45.2", "CD138", "TCRVA2", "TER119")
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% isotype_proteins, ]
Protein_F_pm <- Protein_F_pm[rownames(Protein_F_pm) %in% good_proteins, ]
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% exclude_proteins, ]
thy11_proteins <- grep("THY1.1", rownames(Protein_F_pm), value = TRUE)
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% thy11_proteins, ]
D_lognorm <- diag(1 / apply(Protein_F_pm, 2, function(x) max(abs(x))))
Protein_F_pm <- Protein_F_pm %*% D_lognorm
colnames(Protein_F_pm) <- paste0("GP", seq_len(ncol(Protein_F_pm)))
Protein_F_pm[is.na(Protein_F_pm)] <- 0

# For each GP, list marker proteins with |score| >= threshold, split by sign.
marker_list_to_df <- function(gp_marker_sign_list, collapse = ", ") {
  sets <- names(gp_marker_sign_list)
  data.frame(
    Set = sets,
    Positive = vapply(gp_marker_sign_list, function(x) if (length(x$pos) == 0) "" else paste(x$pos, collapse = collapse), character(1)),
    Negative = vapply(gp_marker_sign_list, function(x) if (length(x$neg) == 0) "" else paste(x$neg, collapse = collapse), character(1)),
    stringsAsFactors = FALSE
  )
}
marker_threshold <- 0.5
gp_marker_sign_list <- list()
for (gp_name in colnames(Protein_F_pm)) {
  marker_values <- Protein_F_pm[, gp_name]
  marker_genes <- rownames(Protein_F_pm)[abs(marker_values) >= marker_threshold]
  marker_genes <- marker_genes[order(-abs(marker_values[marker_genes]))]
  gp_marker_sign_list[[gp_name]] <- list(
    pos = marker_genes[marker_values[marker_genes] >= marker_threshold],
    neg = marker_genes[marker_values[marker_genes] <= -marker_threshold]
  )
}
df_markers <- marker_list_to_df(gp_marker_sign_list)
saveRDS(df_markers, file = paste0(data_path, "CITEseq_markers_full.rds"))
