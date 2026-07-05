# Extended Data Table 5: GP during activation.
#
# One row per GP (GP1..GP200) summarizing how each GP's loading changes with
# T-cell activation, computed separately in CD4 and CD8: the mean change in
# loading (activated minus resting), the average loading over activated+resting
# cells (AveExpr), the standardized mean difference d (= change / pooled SD),
# and the CD8/CD4 ratio of the mean loading changes. This reproduces the
# `GP_activation_summary` object built in script/Figure3.R (panel 3a) from
# code/R/activation_shared_setup.R.
#
# Note vs. the published Table S2: the mean_change_loadings / AveExpr / Ratio
# columns reproduce exactly (identical CD4/CD8 activated-vs-resting populations),
# but the published z_CD4 / z_CD8 were a large-sample z-test statistic from an
# older limma-style DGE that the analysis has since replaced with the
# standardized mean difference d_CD4 / d_CD8 (change / pooled loading SD; see
# activation_shared_setup.R::std_mean_diff) -- the same statistic Figure 3a
# plots. This table therefore carries d_CD4 / d_CD8, the current quantity, not
# the retired z.

library(dplyr)

data_path <- "data/"
output_path <- "figures/generated/"

source("code/R/plot_utils.R") # scale_cols(), used by activation_shared_setup.R
source("code/R/setup_data.R")

gp_data <- load_gp_data(data_path = data_path)
L_pm_filtered <- gp_data$L_pm_filtered
F_pm_filtered <- gp_data$F_pm_filtered
seurat_meta_filtered <- gp_data$seurat_meta_filtered

source("code/R/activation_shared_setup.R") # -> diff_factors_merged, d_factors_merged

GP_activation_summary <- diff_factors_merged %>%
  dplyr::inner_join(d_factors_merged %>% dplyr::select(SYMBOL, d_CD4, d_CD8), by = "SYMBOL") %>%
  dplyr::mutate(Ratio_CD8_CD4 = mean_change_loadings_CD8 / mean_change_loadings_CD4) %>%
  dplyr::select(
    GP = SYMBOL,
    mean_change_loadings_CD4, mean_change_loadings_CD8,
    AveExpr_CD4, AveExpr_CD8,
    d_CD4, d_CD8,
    Ratio_CD8_CD4
  )

write.csv(
  GP_activation_summary,
  file = paste0(output_path, "ExtendedDataTable5_GP_during_activation.csv"),
  row.names = FALSE
)
