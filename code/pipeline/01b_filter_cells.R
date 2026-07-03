# Pipeline step 1b: filter cells by total GP membership.
#
# Unlike 01_extract_data.R (cluster-scale, needs the full raw Seurat
# object), this step only needs flashier_snmf_summary.rds and the Seurat
# metadata, and is runnable against the data/ already in this repo.
#
# Produces L_pm_filtered.rds / F_pm_filtered.rds: restricts
# flashier_snmf_summary.rds's L_pm/F_pm to cells also present in the
# current Seurat object (see note below), normalizes (each column scaled
# to max 1), then iteratively drops cells whose total membership (rowSum
# of the normalized L) exceeds 10, renormalizing after each pass
# (filter_cells_by_total_membership(), code/R/plot_utils.R), and rescales
# F_pm to match. This is the step every figure script's L_pm_filtered.rds /
# F_pm_filtered.rds ultimately comes from.
#
# Source: recovered from a commented-out block in the original
# Figure_Overview.R (the live code had been replaced by
# `readRDS("L_pm_filtered.rds")` once this step's output was cached) --
# confirmed against the same live logic in the original plotUMAPs.R.
#
# Verified byte-identical against the cached data/L_pm_filtered.rds /
# F_pm_filtered.rds (max abs diff 0, same 681,423 cells in the same row
# order). Getting there required the `cells_flashier %in% cells_seurat`
# restriction below: flashier_snmf_summary.rds's L_pm has 682,953 cells,
# but the current Seurat object has only 682,935 -- the same kind of
# data-version drift already flagged for the cached seurat_meta.rds
# (see code/R/setup_data.R), here between flashier_snmf_summary.rds and
# the Seurat object. Omitting this restriction reproduces 18 extra cells
# that survive the membership filter but shouldn't be there.

source("code/R/plot_utils.R") # scale_cols(), filter_cells_by_total_membership()

data_path <- "data/"

flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
cells_flashier <- rownames(flashier_snmf_summary$L_pm)
cells_seurat <- seurat_meta$cellID

L_pm <- flashier_snmf_summary$L_pm[cells_flashier %in% cells_seurat, ]
F_pm <- flashier_snmf_summary$F_pm

D <- diag(1 / apply(L_pm, 2, function(x) max(x)))
L <- L_pm %*% D

cells <- filter_cells_by_total_membership(L, numiter = 12)
L_pm_filtered <- L_pm[cells, ]
d <- apply(L_pm_filtered, 2, max)
L_pm_filtered <- scale_cols(L_pm_filtered, 1 / d)

# Normalize F_pm to match the same per-column scaling applied to L_pm_filtered.
F_pm_filtered <- scale_cols(F_pm, d)

saveRDS(L_pm_filtered, paste0(data_path, "L_pm_filtered.rds"))
saveRDS(F_pm_filtered, paste0(data_path, "F_pm_filtered.rds"))
