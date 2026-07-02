# Figure S2. GP30 and GP58 loadings across T-cell subsets.
#
# Panels produced (see figures/final-selected/bits/Figure S2/FigureS2_caption.md
# for the full caption text):
#   S2A  Boxplots of GP30 loading across Tz subsets (iNKT, MAIT, other Tz)
#        vs all other, non-Tz T cells.
#   S2B  Boxplots of GP58 loading: resting CD8, activated CD8, vs the other
#        six T-cell lineages pooled as "Other T cells".
#   S2C  Boxplots of CD8A/CD8B log-normalized CITE-seq protein expression in
#        CD8 cells, resting vs activated.
#
# Source: ported from Figure_Lineage.R (see Figure2.R for the main
# Figure 2 panels from the same file).

library(ggplot2)
library(dplyr)
library(ggrastr)
library(tidyr)
library(Matrix) # protein_mat_normalized_lognorm is a dgCMatrix; rownames()
                # dispatch on it is unreliable unless Matrix is attached
                # (not just loaded as a namespace), which silently intersected
                # to zero cells during this refactor -- keep this library() call.

data_path <- "data/"
figure_path <- "figures/generated/Figure S2/"
source("code/R/plot_utils.R") # tukey_outliers()

# ============================================================
# Load data
# ============================================================
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
colnames(L_pm_filtered) <- gsub("^K", "GP", colnames(L_pm_filtered))
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
selected_lineage_in_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN")

# ============================================================
# S2A: GP30 loading, Tz subsets (iNKT/MAIT/other Tz) vs other T cells
# ============================================================
GP30_df <- seurat_meta_filtered %>%
  select(cellID, annotation_level1, iNKT, MAIT) %>%
  filter(annotation_level1 != "thymocyte") %>%
  mutate(GP30_loading = L_pm_filtered[cellID, "GP30"]) %>%
  mutate(
    Group = case_when(
      annotation_level1 == "Tz" & iNKT ~ "iNKT",
      annotation_level1 == "Tz" & MAIT ~ "MAIT",
      annotation_level1 == "Tz" ~ "Other Tz",
      TRUE ~ "Other T cells"
    ),
    Group = factor(Group, levels = c("iNKT", "MAIT", "Other Tz", "Other T cells"))
  )

p_S2A <- ggplot(GP30_df, aes(x = Group, y = GP30_loading, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  ggrastr::rasterise(
    geom_point(data = tukey_outliers(GP30_df, "GP30_loading", "Group"), size = 0.5, alpha = 0.3, show.legend = FALSE),
    dpi = 300
  ) +
  scale_fill_manual(values = c("iNKT" = "darkgoldenrod2", "MAIT" = "darkgoldenrod3", "Other Tz" = "darkgoldenrod1", "Other T cells" = "grey70")) +
  labs(title = "GP30 loading across Tz subsets and other T cells", x = NULL, y = "GP30 loading") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.major.x = element_blank(), legend.position = "none")
ggsave(filename = paste0(figure_path, "S2A.pdf"), plot = p_S2A, width = 6, height = 5)

# ============================================================
# S2B: GP58 loading, resting/activated CD8 vs other T cells
# ============================================================
non_CD8_lineages <- setdiff(selected_lineage_in_order, "CD8")
GP58_df <- seurat_meta_filtered %>%
  select(cellID, annotation_level1, annotation_level2_group) %>%
  mutate(GP58_loading = L_pm_filtered[cellID, "GP58"]) %>%
  mutate(
    Group = case_when(
      annotation_level1 == "CD8" & annotation_level2_group == "resting" ~ "CD8 (Resting)",
      annotation_level1 == "CD8" & annotation_level2_group == "activated" ~ "CD8 (Activated)",
      annotation_level1 %in% non_CD8_lineages ~ "Other T cells",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%
  mutate(Group = factor(Group, levels = c("CD8 (Resting)", "CD8 (Activated)", "Other T cells")))

p_S2B <- ggplot(GP58_df, aes(x = Group, y = GP58_loading, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  ggrastr::rasterise(
    geom_point(data = tukey_outliers(GP58_df, "GP58_loading", "Group"), size = 0.5, alpha = 0.3, show.legend = FALSE),
    dpi = 300
  ) +
  scale_fill_manual(values = c("CD8 (Resting)" = "darkorange2", "CD8 (Activated)" = "orange", "Other T cells" = "grey70")) +
  labs(title = "GP58 loading across CD8 subsets and other T cells", x = NULL, y = "GP58 loading") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.major.x = element_blank(), legend.position = "none")
ggsave(filename = paste0(figure_path, "S2B.pdf"), plot = p_S2B, width = 5, height = 5)

# ============================================================
# S2C: CD8A/CD8B CITE-seq protein expression, resting vs activated CD8
# ============================================================
protein_mat_normalized_lognorm <- readRDS(paste0(data_path, "protein_mat_normalized_lognorm.rds"))

CD8_citeseq_cells <- seurat_meta_filtered$cellID[
  seurat_meta_filtered$annotation_level1 == "CD8" &
    seurat_meta_filtered$cite_seq &
    seurat_meta_filtered$annotation_level2_group %in% c("resting", "activated")
]
CD8_citeseq_cells <- intersect(CD8_citeseq_cells, rownames(protein_mat_normalized_lognorm))

CD8_protein_df <- data.frame(
  cellID = CD8_citeseq_cells,
  Group = ifelse(seurat_meta_filtered[CD8_citeseq_cells, "annotation_level2_group"] == "resting", "CD8 (Resting)", "CD8 (Activated)"),
  CD8A = protein_mat_normalized_lognorm[CD8_citeseq_cells, "CD8A"],
  CD8B = protein_mat_normalized_lognorm[CD8_citeseq_cells, "CD8B"]
) %>%
  tidyr::pivot_longer(cols = c("CD8A", "CD8B"), names_to = "Protein", values_to = "Expression") %>%
  mutate(Group = factor(Group, levels = c("CD8 (Resting)", "CD8 (Activated)")))

p_S2C <- ggplot(CD8_protein_df, aes(x = Protein, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  ggrastr::rasterise(
    geom_point(
      data = tukey_outliers(CD8_protein_df, "Expression", c("Protein", "Group")),
      aes(group = Group), size = 0.5, alpha = 0.3, position = position_dodge(width = 0.75), show.legend = FALSE
    ),
    dpi = 300
  ) +
  scale_fill_manual(values = c("CD8 (Resting)" = "darkorange2", "CD8 (Activated)" = "orange")) +
  labs(title = "CD8A / CD8B protein expression in CD8 cells", x = NULL, y = "Log-normalized protein expression", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.major.x = element_blank(), legend.position = "top")
ggsave(filename = paste0(figure_path, "S2C.pdf"), plot = p_S2C, width = 5, height = 5)
