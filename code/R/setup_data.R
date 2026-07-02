# Consolidated data loading for the figure scripts in script/.
#
# Every Figure_*.R in the old script/ directory repeated the same ~15-line
# block at the top: read L_pm_filtered.rds / F_pm_filtered.rds, rename their
# K## columns to GP##, read the cell metadata, and subset it to the filtered
# cells. This factors that block into one function so every figure script
# just calls load_gp_data() and gets a consistent set of objects back.
#
# L_pm_filtered.rds / F_pm_filtered.rds are themselves produced upstream by
# code/pipeline/01_extract_data.R (cell-loading matrix filtered by
# filter_cells_by_total_membership(), see code/R/plot_utils.R).
#
# NOTE: cell metadata is read from the Seurat object's @meta.data, not from
# the cached data/seurat_meta.rds -- that cached file was found during this
# refactor to be a stale snapshot (682,951 cells / 60 columns) that no
# longer matches the current Seurat object (682,935 cells / 66 columns).
# Reading straight from the Seurat object matches what every original
# Figure_*.R script actually did.

load_gp_data <- function(
    data_path = "data/",
    seurat_file = "igt1_96_withtotalvi20260206_clean_ADTonly.Rds",
    load_mde = FALSE,
    load_prior = FALSE
) {
  L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
  colnames(L_pm_filtered) <- paste0("GP", seq_len(ncol(L_pm_filtered)))

  F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
  colnames(F_pm_filtered) <- paste0("GP", seq_len(ncol(F_pm_filtered)))

  seurat_meta <- readRDS(paste0(data_path, seurat_file))@meta.data
  seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]

  out <- list(
    L_pm_filtered = L_pm_filtered,
    F_pm_filtered = F_pm_filtered,
    seurat_meta = seurat_meta,
    seurat_meta_filtered = seurat_meta_filtered
  )

  if (load_mde) {
    mde_result <- readRDS(paste0(data_path, "umap_result.rds"))
    colnames(mde_result) <- c("MDE_1", "MDE_2")
    out$mde_result <- mde_result[rownames(L_pm_filtered), ]
  }

  if (load_prior) {
    # Loads `flashier_snmf_fitted_prior` into this function's environment and
    # returns it — needed for PVE-type panels (e.g. Figure 1's scatter of
    # expected active-cell/active-gene proportion per GP).
    load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
    out$flashier_snmf_fitted_prior <- flashier_snmf_fitted_prior
  }

  out
}
