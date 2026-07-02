# Supplementary Table 1: per-GP summary.
#
# One row per GP: which Level-1/Level-2/condition/organ categories it
# predicts well (AUC > 0.8, driven by high loading), plus its top positive
# and negative signature genes.
#
# Source: ported from Supplement_Table1.R (unchanged apart from the
# output path).

data_path <- "data/"
output_path <- "figures/generated/"

L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
level_1_AUC_list <- readRDS(paste0(data_path, "level_1_AUC_list_figure.rds"))
level_2_AUC_list <- readRDS(paste0(data_path, "level_2_AUC_list_figure.rds"))
condition_detailed_AUC_list_figure <- readRDS(paste0(data_path, "condition_detailed_AUC_list_figure.rds"))
organ_simplified_AUC_list_figure <- readRDS(paste0(data_path, "organ_simplified_AUC_list_figure.rds"))

# Normalize F_pm_filtered such that each column has max abs value of 1
F_pm_filtered <- apply(F_pm_filtered, 2, function(x) x / max(abs(x)))
colnames(F_pm_filtered) <- colnames(L_pm_filtered)

gps <- colnames(F_pm_filtered)

get_passing_categories <- function(auc_list, gps, L_mat, threshold = 0.8) {
  auc_mat <- auc_list$auc
  thr_mat <- auc_list$threshold
  gp_median_loading <- apply(L_mat, 2, median)
  lapply(gps, function(gp) {
    if (!gp %in% colnames(auc_mat)) return("")
    auc_vals <- auc_mat[, gp]
    thr_vals <- thr_mat[, gp]
    med_load <- gp_median_loading[gp]
    # Keep only categories predicted by HIGH loading (threshold >= median loading)
    cats <- rownames(auc_mat)[auc_vals > threshold & thr_vals >= med_load]
    paste(cats, collapse = "; ")
  })
}

level1_cats <- get_passing_categories(level_1_AUC_list, gps, L_pm_filtered)
level2_cats <- get_passing_categories(level_2_AUC_list, gps, L_pm_filtered)
condition_cats <- get_passing_categories(condition_detailed_AUC_list_figure, gps, L_pm_filtered)
organ_cats <- get_passing_categories(organ_simplified_AUC_list_figure, gps, L_pm_filtered)

sig_genes_pos <- lapply(gps, function(gp) {
  vals <- F_pm_filtered[, gp]
  candidates <- names(vals)[vals > 0.25]
  top <- head(candidates[order(vals[candidates], decreasing = TRUE)], 5)
  paste(top, collapse = "; ")
})
sig_genes_neg <- lapply(gps, function(gp) {
  vals <- F_pm_filtered[, gp]
  candidates <- names(vals)[vals < -0.25]
  top <- head(candidates[order(abs(vals[candidates]), decreasing = TRUE)], 5)
  paste(top, collapse = "; ")
})

supp_table <- data.frame(
  GP = gps,
  Level1 = unlist(level1_cats),
  Level2 = unlist(level2_cats),
  Condition = unlist(condition_cats),
  Organ = unlist(organ_cats),
  Signature_Genes_Pos = unlist(sig_genes_pos),
  Signature_Genes_Neg = unlist(sig_genes_neg),
  stringsAsFactors = FALSE
)

write.csv(supp_table, file = paste0(output_path, "Supplementary_Table1_GP_summary.csv"), row.names = FALSE)
