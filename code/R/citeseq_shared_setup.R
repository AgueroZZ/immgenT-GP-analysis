# Shared setup for the Figure 6 / Figure S6 scripts (CITE-seq protein
# projection). Both figure scripts subset to CITE-seq-measured cells, build
# the same curated marker table (df_markers2), the same cell exclusions
# (thymocyte/proliferating/miniverse), and the same protein quality filter
# (select_proteins) -- factored here so both source it once instead of
# repeating ~50 lines of setup each.
#
# Requires data_path to already be set. Defines/overwrites: L_pm_filtered,
# F_pm_filtered, seurat_meta_filtered, protein_mat_normalized_lognorm,
# mde_result, Protein_F_pm_raw, select_proteins, threshold_results_subset_manual,
# df_markers2, thymocyte_cells, proliferating_cells, miniverse_cells,
# well_aligned_gps.

L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(L_pm_filtered) <- paste0("GP", 1:ncol(L_pm_filtered))
colnames(F_pm_filtered) <- paste0("GP", 1:ncol(F_pm_filtered))
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
protein_mat_normalized_lognorm <- readRDS(paste0(data_path, "protein_mat_normalized_lognorm.rds"))
mde_result <- readRDS(paste0(data_path, "umap_result.rds"))
colnames(mde_result) <- c("MDE_1", "MDE_2")
mde_result <- mde_result[rownames(L_pm_filtered), ]

Protein_flash_result <- readRDS(paste0(data_path, "protein_flash_selected_summary_lognorm_backfit200.rds"))
Protein_F_pm_raw <- Protein_flash_result$F_pm
colnames(Protein_F_pm_raw) <- paste0("GP", 1:ncol(Protein_F_pm_raw))

cells_citeseq <- seurat_meta_filtered$cellID[seurat_meta_filtered$cite_seq]
L_pm_filtered <- L_pm_filtered[cells_citeseq, , drop = FALSE]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[cells_citeseq, , drop = FALSE]
mde_result <- mde_result[cells_citeseq, , drop = FALSE]
seurat_meta_filtered <- seurat_meta_filtered[cells_citeseq, , drop = FALSE]

isotype_proteins <- grep("^Isotype", rownames(Protein_F_pm_raw), value = TRUE)
proteins_quality <- read.csv(paste0(data_path, "TableS4_citeseq_qc_20250513.csv"), header = TRUE, stringsAsFactors = FALSE, skip = 1)
good_proteins <- c(proteins_quality$protein[proteins_quality$classification == "good"], "IL2RA.CD25", "ITB7", "CD69")
exclude_proteins <- c("CD19", "CD34", "CD45.1", "CD45.2", "CD138", "TCRVA2", "TER119")
thy11_proteins <- grep("THY1.1", rownames(Protein_F_pm_raw), value = TRUE)
select_proteins <- setdiff(good_proteins, c(exclude_proteins, thy11_proteins, isotype_proteins))

threshold_results_subset <- read.csv(paste0(data_path, "Thresholds_Selected_Proteins.csv"), header = TRUE, stringsAsFactors = FALSE)
threshold_results_subset <- threshold_results_subset[, c("Protein", "Threshold", "Threshold_manual")]
threshold_results_subset <- rbind(threshold_results_subset, data.frame(Protein = "CD62L", Threshold = 3, Threshold_manual = 3))
threshold_results_subset_manual <- data.frame(Protein = threshold_results_subset$Protein, Threshold = threshold_results_subset$Threshold_manual)

# Curated marker overrides for specific GPs (manually reviewed), on top of
# the automatically-derived df_markers.
df_markers <- readRDS(paste0(data_path, "CITEseq_markers_full.rds"))
df_markers2 <- df_markers
df_markers2$Positive[8] <- "TCRVG3"
df_markers2$Negative[8] <- ""
df_markers2$Negative[22] <- "CD4, CD8B, TCRGD"
df_markers2$Positive[23] <- "ITB7, CD103, CD44"
df_markers2$Negative[23] <- "CD29, ITA4.CD49D"
df_markers2$Positive[29] <- "CD8A"
df_markers2$Negative[29] <- "CD4, CD8B"
df_markers2$Positive[30] <- "CD11A, CD49B, CD38"
df_markers2$Negative[30] <- "CD8B, CD8A, CD4"
df_markers2$Positive[41] <- "CD55.DAF, CD4, CD44, CD45RB, CD62L, CD31, CD2, IL7RA.CD127, CD27, SCA1, CD5"
df_markers2$Positive[57] <- "CD44, ICOS.CD278"
df_markers2$Negative[57] <- "CD62L"
df_markers2$Positive[68] <- "IL2RA.CD25, FR4, GITR.CD357, NEUROPILIN1.CD304"
df_markers2$Negative[68] <- ""
df_markers2$Positive[80] <- "CD29, CD44, ITA4.CD49D"
df_markers2$Negative[80] <- "ITB7, CD103"
df_markers2$Positive[170] <- "ITB7, CD103, CD4, CD38"
df_markers2$Negative[170] <- "GITR.CD357, CD62L"
df_markers2$Positive[171] <- "CD62L"
df_markers2$Negative[171] <- "CD44"

thymocyte_cells <- seurat_meta_filtered %>% dplyr::filter(annotation_level1 == "thymocyte") %>% dplyr::pull(cellID)
proliferating_cells <- seurat_meta_filtered %>% dplyr::filter(annotation_level2_group == "proliferating") %>% dplyr::pull(cellID)
miniverse_cells <- seurat_meta_filtered %>% dplyr::filter(annotation_level2_group == "miniverse") %>% dplyr::pull(cellID)

# GPs judged well-aligned between protein gating and GP loading (manually
# curated; see data/CITEseq_alignment_scores_manual.csv), used by Figure6.R's
# panel 6b heatmap and FigureS6.R's gallery.
well_aligned_gps <- c(
  "GP10", "GP26", "GP68", "GP171", "GP58", "GP8", "GP30", "GP27", "GP170", "GP80",
  "GP35", "GP12", "GP3", "GP29", "GP77", "GP22", "GP25", "GP41", "GP181", "GP63",
  "GP107", "GP23", "GP126", "GP127", "GP192", "GP159", "GP153"
)

# plot_gated_gp_vs_protein()/get_gated_cell_ids() expect the loading matrix
# with "K##" colnames (they map internally via gp_label()); keep a separate
# copy so the "GP##" renaming above isn't disturbed for callers that want it.
L_pm_for_gating <- L_pm_filtered
colnames(L_pm_for_gating) <- paste0("K", seq_len(ncol(L_pm_for_gating)))
