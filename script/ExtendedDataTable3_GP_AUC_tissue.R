# Extended Data Table 3: GP AUC by tissue.
#
# One row per GP (GP1..GP200), one column per tissue (organ_simplified: blood,
# spleen, LN, bone marrow, lung, liver, ...), holding the one-vs-rest AUC for
# predicting membership in that tissue from the GP's loading. Restricted to
# healthy non-thymocyte cells (the *_no_thymocytes_healthy AUC family, same as
# Figure 4), extracted directly from the precomputed AUC matrices -- see
# code/pipeline/02_compute_auc.R.
#
# The AUC file's $auc is a categories x GPs matrix with raw "K##" GP names;
# auc_list_to_gp_table() transposes it to a per-GP table and relabels K -> GP.

data_path <- "data/"
output_path <- "figures/generated/"

source("code/R/roc_auc.R")

organ_simplified_AUC_list <- readRDS(paste0(data_path, "organ_simplified_AUC_list_figure_no_thymocytes_healthy.rds"))

gp_auc_tissue <- auc_list_to_gp_table(organ_simplified_AUC_list)

write.csv(
  gp_auc_tissue,
  file = paste0(output_path, "ExtendedDataTable3_GP_AUC_tissue.csv"),
  row.names = FALSE
)
