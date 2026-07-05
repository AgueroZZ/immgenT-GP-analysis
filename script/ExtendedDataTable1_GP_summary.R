# Extended Data Table 1: summary of GP characteristics and annotations.
#
# One row per GP (GP1..GP200) with:
#   - Lineage / Cluster / Tissue: the categories each GP *positively* predicts
#     well -- one-vs-rest AUC > 0.8 AND the optimal decision threshold at or
#     above the GP's median loading, so high loading (not low) drives the
#     prediction. Categories are annotation_level1 (lineage), annotation_level2
#     (cluster) and organ_simplified (tissue), all from the healthy
#     non-thymocyte AUC family (*_no_thymocytes_healthy, same population as
#     Figure 2/Figure 4); the median loading is likewise computed on the healthy
#     non-thymocyte cells so the threshold comparison is on the same population.
#   - Signature genes / proteins: the top 5 up- and top 5 down-regulated genes
#     and proteins by factor score, among those with |score| > 0.1 on the
#     max|.|=1-per-GP-scaled gene (F_pm_filtered) and protein (Figure 6's
#     Protein_F_pm) factor matrices.
#
# Reworks the retired Table S1 (script/TableS1.R): drops its Condition column,
# switches the annotations from the non-thymocyte healthy+diseased AUC to the
# healthy non-thymocyte AUC, adds protein signatures, and loosens the gene
# signature cutoff from 0.25 to 0.1.

data_path <- "data/"
output_path <- "figures/generated/"

# ---- Protein factor matrix (Figure 6's Protein_F_pm: filtered + max|.|=1 per
# GP column). Load, extract, and drop the heavy fit before loading L. ----
Protein_flash_result <- readRDS(paste0(data_path, "protein_flash_selected_summary_lognorm_backfit200.rds"))
Protein_F_pm_raw <- Protein_flash_result$F_pm
rm(Protein_flash_result); gc()

isotype_proteins <- grep("^Isotype", rownames(Protein_F_pm_raw), value = TRUE)
proteins_quality <- read.csv(paste0(data_path, "TableS4_citeseq_qc_20250513.csv"), header = TRUE, stringsAsFactors = FALSE, skip = 1)
good_proteins <- c(proteins_quality$protein[proteins_quality$classification == "good"], "IL2RA.CD25", "ITB7", "CD69")
exclude_proteins <- c("CD19", "CD34", "CD45.1", "CD45.2", "CD138", "TCRVA2", "TER119")
thy11_proteins <- grep("THY1.1", rownames(Protein_F_pm_raw), value = TRUE)

Protein_F_pm <- Protein_F_pm_raw[!rownames(Protein_F_pm_raw) %in% isotype_proteins, ]
Protein_F_pm <- Protein_F_pm[rownames(Protein_F_pm) %in% good_proteins, ]
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% exclude_proteins, ]
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% thy11_proteins, ]
D_lognorm <- diag(1 / apply(Protein_F_pm, 2, function(x) max(abs(x))))
Protein_F_pm <- Protein_F_pm %*% D_lognorm
Protein_F_pm[is.na(Protein_F_pm)] <- 0

# ---- Loadings, gene factor matrix, healthy non-thymocyte AUC ----
# L_pm_filtered / F_pm_filtered and the cached AUC matrices all share the raw
# "K1","K2",... column names (from the flashier fit); keep them for matching and
# relabel to "GP1","GP2",... only for output.
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
level_1_AUC_list <- readRDS(paste0(data_path, "level_1_AUC_list_figure_no_thymocytes_healthy.rds"))
level_2_AUC_list <- readRDS(paste0(data_path, "level_2_AUC_list_figure_no_thymocytes_healthy.rds"))
organ_AUC_list   <- readRDS(paste0(data_path, "organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds"))

# Use L_pm_filtered's "K1".."K200" names as the canonical GP key: the AUC
# matrices carry these same names, whereas F_pm_filtered's raw columns are
# "F1".."F200" and the protein factor columns are unnamed -- but all three are
# in the same factor order as L (same flashier fit / fixed-loading projection),
# so relabeling them to L's K-names aligns every signature/annotation lookup.
gps <- colnames(L_pm_filtered)          # "K1".."K200"
gp_labels <- paste0("GP", seq_along(gps))
colnames(Protein_F_pm) <- gps

# Normalize the gene factor matrix so each GP column has max|score| = 1.
F_pm_filtered <- apply(F_pm_filtered, 2, function(x) x / max(abs(x)))
colnames(F_pm_filtered) <- gps

# ---- Median loading on the healthy non-thymocyte cell subset ----
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
rm(seurat_meta); gc()
# which() drops any cells whose condition/lineage metadata is NA (so they don't
# leak in as NA rows and turn every median into NA). seurat_meta_filtered is row-
# aligned to L_pm_filtered, so the integer index selects the matching L rows.
healthy_nonthy_idx <- which(
  seurat_meta_filtered$condition_broad == "healthy" &
    seurat_meta_filtered$annotation_level1 != "thymocyte"
)
L_healthy_nonthy <- L_pm_filtered[healthy_nonthy_idx, , drop = FALSE]

# ---- Positively-predicted categories: AUC > 0.8 and threshold >= median
# loading (so high loading drives the prediction) ----
get_passing_categories <- function(auc_list, gps, L_mat, threshold = 0.8) {
  auc_mat <- auc_list$auc
  thr_mat <- auc_list$threshold
  gp_median_loading <- apply(L_mat, 2, median)
  lapply(gps, function(gp) {
    if (!gp %in% colnames(auc_mat)) return("")
    auc_vals <- auc_mat[, gp]
    thr_vals <- thr_mat[, gp]
    med_load <- gp_median_loading[gp]
    pass <- !is.na(auc_vals) & !is.na(thr_vals) &
      auc_vals > threshold & thr_vals >= med_load
    cats <- rownames(auc_mat)[pass]
    paste(cats, collapse = "; ")
  })
}

lineage_cats <- get_passing_categories(level_1_AUC_list, gps, L_healthy_nonthy)
cluster_cats <- get_passing_categories(level_2_AUC_list, gps, L_healthy_nonthy)
tissue_cats  <- get_passing_categories(organ_AUC_list, gps, L_healthy_nonthy)

# ---- Top 5 up / down signature genes and proteins (|score| > 0.1) ----
top_signatures <- function(score_mat, gps, cutoff = 0.1, n = 5) {
  pos <- lapply(gps, function(gp) {
    vals <- score_mat[, gp]
    cand <- names(vals)[vals > cutoff]
    top <- head(cand[order(vals[cand], decreasing = TRUE)], n)
    paste(top, collapse = "; ")
  })
  neg <- lapply(gps, function(gp) {
    vals <- score_mat[, gp]
    cand <- names(vals)[vals < -cutoff]
    top <- head(cand[order(abs(vals[cand]), decreasing = TRUE)], n)
    paste(top, collapse = "; ")
  })
  list(pos = unlist(pos), neg = unlist(neg))
}

gene_sig <- top_signatures(F_pm_filtered, gps)
prot_sig <- top_signatures(Protein_F_pm, gps)

supp_table <- data.frame(
  GP = gp_labels,
  Lineage = unlist(lineage_cats),
  Cluster = unlist(cluster_cats),
  Tissue = unlist(tissue_cats),
  Signature_Genes_Pos = gene_sig$pos,
  Signature_Genes_Neg = gene_sig$neg,
  Signature_Proteins_Pos = prot_sig$pos,
  Signature_Proteins_Neg = prot_sig$neg,
  stringsAsFactors = FALSE
)

write.csv(
  supp_table,
  file = paste0(output_path, "ExtendedDataTable1_GP_summary.csv"),
  row.names = FALSE
)
