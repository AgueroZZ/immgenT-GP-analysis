# Figure 2. Active cells and genes in gene programs.
#
# Panels produced (final files figures/generated/Figure 2/2A.pdf .. 2F.pdf):
#   2A  "Signature volcano" for GP1: each gene's normalized weight vs. mean
#       shifted-log expression, top genes labeled (via code/R/volcano_helpers.R;
#       ported from the "immgen-signature" Shiny app -- see that header for the
#       caveat about n_label/threshold not being exactly recoverable).
#   2B  Histogram of the proportion of highly-active cells per GP.
#   2C  Histogram of the number of highly-active genes per GP.
#   2D  Boxplot of active-GP-count per cell, by major lineage.
#   2E  Boxplot of active-GP-count per cell, activated vs resting.
#   2F  Scatter of CD44 protein level vs. number of active GPs per cell.
# Two extra non-lettered diagnostics (raw active-cell-count histogram, active-gene
# proportion histogram, and a PVE-colored active-cell-vs-active-gene scatter) are
# saved under descriptive names. Split out of the former Figure 1 (panels 1D-1I).
#
# Required inputs (data/): L_pm_filtered.rds, F_pm_filtered.rds,
#   igt1_96_..._ADTonly.Rds, protein_mat_normalized_lognorm.rds,
#   mean_shifted_log_expr.rds, flashier_snmf_summary.rds,
#   flashier_snmf_fitted_prior.rda.

library(ggplot2)
library(dplyr)
library(Matrix)  # protein matrix is a dgCMatrix; attach for `[` dispatch

data_path <- "data/"
figure_path <- "figures/generated/Figure 2/"
source("code/R/volcano_helpers.R") # plot_gp_signature_volcano(), normalize_maxabs()

# ============================================================
# Load data (non-thymocyte cells only, matching the former Figure 1)
# ============================================================
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
protein_mat_normalized_lognorm <- readRDS(paste0(data_path, "protein_mat_normalized_lognorm.rds"))
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[rownames(L_pm_filtered), "CD44"]
level1_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN")

non_thymo_cells <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 != "thymocyte"]
L_pm_filtered <- L_pm_filtered[non_thymo_cells, ]
seurat_meta_filtered <- seurat_meta_filtered[non_thymo_cells, ]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[non_thymo_cells]

# ============================================================
# 2A: GP1 signature volcano (see code/R/volcano_helpers.R header)
# ============================================================
F_pm_filtered_1d <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F_pm_filtered_1d) <- paste0("GP", seq_len(ncol(F_pm_filtered_1d)))
F_pm_normalized_1d <- normalize_maxabs(F_pm_filtered_1d)
mean_shifted_log_expr <- readRDS(paste0(data_path, "mean_shifted_log_expr.rds"))
p_2A <- plot_gp_signature_volcano("GP1", F_pm_normalized_1d, mean_shifted_log_expr, threshold = 0.1, n_label = 47, bg_alpha = 0.2)
ggsave(filename = paste0(figure_path, "2A.pdf"), plot = p_2A, width = 7, height = 5.5)

# ============================================================
# 2B: histogram of the proportion of highly-active cells per GP
# ============================================================
L_pm_norm_col <- L_pm_filtered / matrix(apply(L_pm_filtered, 2, function(x) max(x)), nrow = nrow(L_pm_filtered), ncol = ncol(L_pm_filtered), byrow = TRUE)
gp_active_cell_counts <- colSums((L_pm_norm_col) > 1e-1)
gp_active_cell_prop <- gp_active_cell_counts / nrow(L_pm_norm_col)
p_2B <- ggplot(data.frame(prop = gp_active_cell_prop), aes(x = prop)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  scale_x_log10(labels = scales::label_percent()) +
  annotation_logticks(sides = "b") +
  labs(x = "Proportion of highly active cells per GP (log scale)", y = "Count",
       title = "Histogram of highly active cells per GP (proportion)") +
  theme_minimal(base_size = 13)
ggsave(filename = paste0(figure_path, "2B.pdf"), plot = p_2B, width = 6, height = 4, dpi = 300)

# ============================================================
# 2C: histogram of the number of highly-active genes per GP
# ============================================================
F_pm_norm_col <- F_pm_filtered / matrix(apply(F_pm_filtered, 2, function(x) max(abs(x))), nrow = nrow(F_pm_filtered), ncol = ncol(F_pm_filtered), byrow = TRUE)
gp_active_gene_counts <- colSums(abs(F_pm_norm_col) > 0.25)
p_2C <- ggplot(data.frame(count = gp_active_gene_counts), aes(x = count)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  scale_x_log10(labels = scales::label_comma()) +
  annotation_logticks(sides = "b") +
  labs(x = "Number of highly active genes per GP (log scale)", y = "Count",
       title = "Histogram of highly active genes per GP") +
  theme_minimal(base_size = 13)
ggsave(filename = paste0(figure_path, "2C.pdf"), plot = p_2C, width = 6, height = 4, dpi = 300)

# ============================================================
# 2D: boxplot of active-GP-count per cell, by major lineage
# ============================================================
gp_active_cell_counts_level1 <- dplyr::bind_rows(lapply(level1_order, function(grp) {
  cells_grp <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 == grp]
  L_grp <- L_pm_filtered[seurat_meta_filtered$cellID %in% cells_grp, , drop = FALSE]
  data.frame(Group = grp, Active_Cell_Counts = rowSums(L_grp > 1e-1))
}))
group_counts <- gp_active_cell_counts_level1 %>%
  group_by(Group) %>% summarise(n = n(), .groups = "drop") %>%
  mutate(Group_Label = paste0(Group, "\n(n=", n, ")")) %>%
  mutate(Group = factor(Group, levels = level1_order)) %>%
  arrange(Group) %>% mutate(Group_Label = factor(Group_Label, levels = Group_Label))
plot_df_g <- gp_active_cell_counts_level1 %>% left_join(group_counts, by = "Group")
p_2D <- ggplot(plot_df_g, aes(x = Group_Label, y = Active_Cell_Counts, fill = Group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.6, alpha = 0.8, color = "gray40") +
  scale_fill_manual(values = ZemmourLib::immgent_colors$level1) +
  labs(title = "Active Gene Programs per Group", x = "Cell Group (Annotation Level 1)", y = "Number of highly active GPs") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5), axis.text.x = element_text(size = 11, angle = 45, hjust = 1), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 13, face = "bold"), legend.position = "none", panel.grid.minor = element_blank())
ggsave(filename = paste0(figure_path, "2D.pdf"), plot = p_2D, width = 6, height = 4, dpi = 300)

# ============================================================
# 2E: boxplot of active-GP-count per cell, activated vs resting
# ============================================================
cells_activated <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level2_group == "activated"]
gp_active_cell_counts_activated <- rowSums(L_pm_filtered[seurat_meta_filtered$cellID %in% cells_activated, ] > 1e-1)
cells_resting <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level2_group == "resting"]
gp_active_cell_counts_resting <- rowSums(L_pm_filtered[seurat_meta_filtered$cellID %in% cells_resting, ] > 1e-1)
gp_active_cell_counts_df <- data.frame(
  Group = c(rep("Activated", length(gp_active_cell_counts_activated)), rep("Resting", length(gp_active_cell_counts_resting))),
  Active_Cell_Counts = c(gp_active_cell_counts_activated, gp_active_cell_counts_resting))
p_2E <- ggplot(gp_active_cell_counts_df, aes(x = Group, y = Active_Cell_Counts, fill = Group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.6, alpha = 0.8, color = "gray40") +
  scale_fill_manual(values = c("Activated" = "#1f78b4", "Resting" = "#e31a1c")) +
  labs(title = "", x = "", y = "Number of highly active GPs") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 13, face = "bold"), legend.position = "none")
ggsave(filename = paste0(figure_path, "2E.pdf"), plot = p_2E, width = 6, height = 4, dpi = 300)

# ============================================================
# 2F: scatter of CD44 protein level vs. number of active GPs per cell
# ============================================================
gp_active_cell_counts_all <- rowSums(L_pm_filtered > 1e-1)
gp_cd44_df <- data.frame(CD44_Protein_Level = protein_mat_normalized_lognorm, Active_GP_Counts = gp_active_cell_counts_all)
set.seed(123)
df_nz <- gp_cd44_df %>% dplyr::filter(CD44_Protein_Level > 0) %>% sample_n(min(10000, nrow(.)))
R <- cor(df_nz$CD44_Protein_Level, df_nz$Active_GP_Counts, use = "complete.obs")
p_2F <- ggplot(df_nz, aes(CD44_Protein_Level, Active_GP_Counts)) +
  geom_point(alpha = 0.1, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "CD44 Protein Level vs Number of Active GPs", x = "CD44 Protein Level (log-normalized)", y = "Number of Active GPs") +
  annotate("text", x = min(df_nz$CD44_Protein_Level, na.rm = TRUE) + 0.5, y = max(df_nz$Active_GP_Counts, na.rm = TRUE) - 1, label = paste0("R = ", round(R, 2)), size = 4) +
  theme_minimal(base_size = 13)
ggsave(filename = paste0(figure_path, "2F.pdf"), plot = p_2F, width = 6, height = 4, dpi = 300)

# ============================================================
# Extra (non-lettered) diagnostics kept under descriptive names
# ============================================================
pdf(paste0(figure_path, "hist_active_cells_per_GP.pdf"), width = 6, height = 4, useDingbats = FALSE)
hist(gp_active_cell_counts, breaks = 100, xlab = "Number of highly active cells per GP", main = "Histogram of highly active cells per GP", freq = TRUE)
dev.off()

gp_active_gene_prop <- gp_active_gene_counts / nrow(F_pm_norm_col)
pdf(paste0(figure_path, "hist_active_genes_prop_per_GP.pdf"), width = 6, height = 4, useDingbats = FALSE)
hist(gp_active_gene_prop, breaks = 100, xlab = "Proportion of highly active genes per GP", main = "Histogram of highly active genes per GP (proportion)", freq = TRUE)
dev.off()

pve <- flashier_snmf_summary$pve
pve_log <- log10(pve)
pve_log_min <- quantile(pve_log, 0.02, na.rm = TRUE); pve_log_max <- quantile(pve_log, 0.98, na.rm = TRUE)
pve_log_clipped <- pmin(pmax(pve_log, pve_log_min), pve_log_max)
pve_scaled <- (pve_log_clipped - min(pve_log_clipped)) / (max(pve_log_clipped) - min(pve_log_clipped))
n_col <- 100; pal <- colorRampPalette(c("lightgray", "skyblue", "darkblue"))(n_col)
col_points <- pal[pmin(pmax(1, floor(pve_scaled * (n_col - 1)) + 1), n_col)]
load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
p_cells <- 1 - sapply(1:200, function(i) flashier_snmf_fitted_prior$L_ghat[[i]]$pi[1])
p_genes <- 1 - sapply(1:200, function(i) flashier_snmf_fitted_prior$F_ghat[[i]]$pi[1])
pdf(paste0(figure_path, "scatter_e_active_cells_vs_genes.pdf"), width = 6, height = 6)
plot(p_cells, p_genes, xlab = "Expected proportion of active cells", ylab = "Expected proportion of active genes",
     log = "xy", col = col_points, pch = 19, xaxt = "n", yaxt = "n",
     main = "Expected proportion of active cells vs genes per GP\n(colored by PVE)")
ticks <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 1)
axis(1, at = ticks, labels = paste0(ticks * 100, "%")); axis(2, at = ticks, labels = paste0(ticks * 100, "%"))
dev.off()
