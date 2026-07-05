# Extended Data Table 4: GP AUC by cluster.
#
# One row per GP (GP1..GP200), one column per sub-lineage cluster
# (annotation_level2: CD8.A, CD8.B, ..., CD4.A, ...), holding the one-vs-rest
# AUC for predicting membership in that cluster from the GP's loading.
# Restricted to healthy non-thymocyte cells (the *_no_thymocytes_healthy AUC
# family, same as Figure 4), extracted directly from the precomputed AUC
# matrices -- see code/pipeline/02_compute_auc.R.
#
# The AUC file's $auc is a categories x GPs matrix with raw "K##" GP names;
# auc_list_to_gp_table() transposes it to a per-GP table and relabels K -> GP.

data_path <- "data/"
output_path <- "figures/generated/"

source("code/R/roc_auc.R")

level_2_AUC_list <- readRDS(paste0(data_path, "level_2_AUC_list_figure_no_thymocytes_healthy.rds"))

gp_auc_cluster <- auc_list_to_gp_table(level_2_AUC_list)

write.csv(
  gp_auc_cluster,
  file = paste0(output_path, "ExtendedDataTable4_GP_AUC_cluster.csv"),
  row.names = FALSE
)
