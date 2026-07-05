# Extended Data Table 2: GP AUC by lineage.
#
# One row per GP (GP1..GP200), one column per major lineage (annotation_level1:
# CD8, CD4, Treg, gdT, CD8aa, Tz, DN, DP), holding the one-vs-rest AUC for
# predicting membership in that lineage from the GP's loading. Restricted to
# healthy non-thymocyte cells (the *_no_thymocytes_healthy AUC family, same as
# Figure 2/Figure 4), extracted directly from the precomputed AUC matrices --
# see code/pipeline/02_compute_auc.R.
#
# The AUC file's $auc is a categories x GPs matrix with raw "K##" GP names;
# auc_list_to_gp_table() transposes it to a per-GP table and relabels K -> GP.

data_path <- "data/"
output_path <- "figures/generated/"

source("code/R/roc_auc.R")

level_1_AUC_list <- readRDS(paste0(data_path, "level_1_AUC_list_figure_no_thymocytes_healthy.rds"))

gp_auc_lineage <- auc_list_to_gp_table(level_1_AUC_list)

write.csv(
  gp_auc_lineage,
  file = paste0(output_path, "ExtendedDataTable2_GP_AUC_lineage.csv"),
  row.names = FALSE
)
