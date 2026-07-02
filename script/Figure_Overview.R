library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(fastTopics)
library(dplyr)
library(pheatmap)
library(scales)
library(scattermore)

#####################################################
#####################################################
####################################################
##### Defining directory and loading functions
####################################################
#####################################################
# data_path <- "../data/"
# code_path <- "../code/"
data_path <- "data/"
code_path <- "code/"
figure_path <- "figures/Figure1_Overview/"
source(paste0(code_path, "ROC.R"))
source(paste0(code_path, "filtering_membership.R"))


#####################################################
#####################################################
#####################################################
### Pre-process data and save them (only need once)
#####################################################
#####################################################
# seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
# flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
# cells_flashier <- rownames(flashier_snmf_summary$L_pm)
# cells_seurat <- seurat_meta$cellID
# # check if all cells in flashier are in seurat
# if (!all(cells_flashier %in% cells_seurat)) {
#   stop("Not all cells in flashier are in seurat")
# }
# L_pm <- flashier_snmf_summary$L_pm[cells_flashier %in% cells_seurat, ]
# D <- apply(L_pm,2,max)
# L_pm_norm <- scale_cols(L_pm,1/D)
# F_pm <- flashier_snmf_summary$F_pm
# F_pm_norm <- scale_cols(F_pm,D)
# colnames(L_pm) <- paste0("K", seq_len(ncol(L_pm)))
# colnames(F_pm) <- paste0("F", seq_len(ncol(F_pm)))
# D <- diag(1 / apply(L_pm, 2, function(x) max(x)))
# L <- L_pm %*% D
# cells <- filter_cells_by_total_membership(L,numiter = 12)
# seurat_meta_filtered <- seurat_meta[cells,]
# L_pm_filtered <- L_pm[cells,]
# d <- apply(L_pm_filtered,2,max)
# L_pm_filtered <- scale_cols(L_pm_filtered,1/d)
# # normalize F correspondingly
# F_pm_filtered <- scale_cols(F_pm,d)
# L_pm_no_thymocytes_normalized <- L_pm_filtered[seurat_meta_filtered$annotation_level1 != "thymocyte", ]
# seurat_meta_filtered_no_thymocytes <- seurat_meta_filtered[seurat_meta_filtered$annotation_level1 != "thymocyte", ]
# saveRDS(L_pm_filtered, paste0(data_path, "L_pm_filtered.rds"))
# saveRDS(F_pm_filtered, paste0(data_path, "F_pm_filtered.rds"))

#####################################################
#####################################################
#####################################################
### Loading Data
#####################################################
#####################################################
load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
seurat_meta <- readRDS(paste0(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
protein_mat_normalized_lognorm <- readRDS(paste0(
  data_path,
  "protein_mat_normalized_lognorm.rds"
))
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[
  rownames(L_pm_filtered),
  "CD44"
]

# Drop thymocytes from all downstream cell-level visualizations.
non_thymo_cells <- seurat_meta_filtered$cellID[
  seurat_meta_filtered$annotation_level1 != "thymocyte"
]
L_pm_filtered <- L_pm_filtered[non_thymo_cells, ]
seurat_meta_filtered <- seurat_meta_filtered[non_thymo_cells, ]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[
  non_thymo_cells
]


#####################################################
#####################################################
#####################################################
### A giant heatmap showing loadings of 200 GPs
#####################################################
#####################################################
# Rows: stratified sample from top organs per level1, unioned with top-K loading anchor cells
#        per GP (to ensure sparse GPs are not all-white). Ordered biologically, not clustered.
# Columns: GPs clustered by ward.D2 to reveal modular structure.
library(ComplexHeatmap)
library(circlize)
library(dplyr)

set.seed(6173)
MIN_CELLS <- 20 # minimum cells in a level1 x organ combo to include
N_SAMPLE <- 80 # cells to sample per combo
TOP_ORGANS <- 5 # top organs per level1 by cell count
K_ANCHOR <- 5 # top-loading anchor cells per GP (ensures sparse GPs are visible)

level1_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN")

# Organ ordering: ensure spleen and LN are adjacent; other organ order is flexible.
organs_all <- as.character(unique(seurat_meta_filtered$organ_simplified))
ln_match <- organs_all[grepl("^LN$|lymph", organs_all, ignore.case = TRUE)]
spleen_match <- organs_all[grepl("spleen", organs_all, ignore.case = TRUE)]
other_organs <- sort(setdiff(organs_all, c(ln_match, spleen_match)))
organ_order <- c(spleen_match, ln_match, other_organs)

# Step 1: Top organs per level1 (exclude levels not in level1_order, e.g. thymocyte)
top_organ_combos <- seurat_meta_filtered |>
  dplyr::filter(annotation_level1 %in% level1_order) |>
  dplyr::count(annotation_level1, organ_simplified) |>
  dplyr::group_by(annotation_level1) |>
  dplyr::slice_max(n, n = TOP_ORGANS, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(annotation_level1, organ_simplified)

# Step 2: Stratified random sampling restricted to top organs
sampled_random <- seurat_meta_filtered |>
  dplyr::filter(annotation_level1 %in% level1_order) |>
  dplyr::inner_join(
    top_organ_combos,
    by = c("annotation_level1", "organ_simplified")
  ) |>
  dplyr::group_by(annotation_level1, organ_simplified) |>
  dplyr::filter(dplyr::n() >= MIN_CELLS) |>
  dplyr::slice_sample(n = N_SAMPLE) |>
  dplyr::ungroup()

# Step 3: Anchor cells — top-K loading per GP, restricted to top organ combos
anchor_cellids <- apply(L_pm_filtered, 2, function(x) {
  rownames(L_pm_filtered)[order(x, decreasing = TRUE)[seq_len(K_ANCHOR)]]
}) |>
  as.vector() |>
  unique()

anchor_meta <- seurat_meta_filtered |>
  dplyr::filter(
    cellID %in% anchor_cellids,
    annotation_level1 %in% level1_order
  ) |>
  dplyr::inner_join(
    top_organ_combos,
    by = c("annotation_level1", "organ_simplified")
  )

# Step 4: Union and order by level1 -> organ
all_meta <- dplyr::bind_rows(sampled_random, anchor_meta) |>
  dplyr::distinct(cellID, .keep_all = TRUE) |>
  dplyr::arrange(
    factor(annotation_level1, levels = level1_order),
    factor(organ_simplified, levels = organ_order)
  )

L_sampled <- L_pm_filtered[all_meta$cellID, ]
clip_val <- quantile(L_sampled, 0.99)
L_display <- pmin(L_sampled, clip_val)
colnames(L_display) <- gsub("^K", "GP", colnames(L_display))

col_fun <- colorRamp2(
  c(0, clip_val / 2, clip_val),
  c("white", "#4393c3", "#08306b")
)

level1_colors <- ZemmourLib::immgent_colors$level1
organ_colors <- ZemmourLib::immgent_colors$organ_simplified
# Order legend entries to match level1_order / organ_order (factor levels drive legend order).
level1_present <- intersect(level1_order, as.character(unique(all_meta$annotation_level1)))
organ_present  <- intersect(organ_order,  as.character(unique(all_meta$organ_simplified)))

row_ann <- rowAnnotation(
  Cell_Type = factor(as.character(all_meta$annotation_level1), levels = level1_present),
  Organ     = factor(as.character(all_meta$organ_simplified),  levels = organ_present),
  col = list(
    Cell_Type = level1_colors[level1_present],
    Organ = organ_colors[organ_present]
  ),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    Cell_Type = list(title = "Cell Type"),
    Organ = list(title = "Organ")
  )
)

ht <- Heatmap(
  L_display,
  name = "Loading",
  col = col_fun,
  left_annotation = row_ann,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 4),
  column_title = "Gene Programs (GPs)",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  use_raster = TRUE,
  raster_quality = 3,
  border = FALSE,
  heatmap_legend_param = list(title = "Loading", direction = "vertical")
)

pdf(
  paste0(figure_path, "heatmap_loading_overview.pdf"),
  width = 15,
  height = 20,
  useDingbats = FALSE
)
draw(ht, merge_legend = TRUE)
dev.off()


#####################################################
#####################################################
#####################################################
### Heatmap showing loadings by level1 then level2
#####################################################
#####################################################
# Same logic as above but rows grouped by level1 x level2 (cell subtype),
# ordered by level1 canonical order then by descending cell count within each level1.
set.seed(6173)
MIN_CELLS_L2 <- 20 # minimum cells in a level2 group to sample from
N_SAMPLE_L2 <- 40 # cells per level2 group
K_ANCHOR_L2 <- 5 # top-loading anchor cells per GP

# Order level2 within each level1 by descending cell count
level2_order <- seurat_meta_filtered |>
  dplyr::filter(annotation_level1 %in% level1_order) |>
  dplyr::count(annotation_level1, annotation_level2) |>
  dplyr::group_by(annotation_level1) |>
  dplyr::arrange(dplyr::desc(n), .by_group = TRUE) |>
  dplyr::ungroup() |>
  dplyr::pull(annotation_level2)

sampled_random_l2 <- seurat_meta_filtered |>
  dplyr::filter(annotation_level1 %in% level1_order) |>
  dplyr::group_by(annotation_level2) |>
  dplyr::filter(dplyr::n() >= MIN_CELLS_L2) |>
  dplyr::slice_sample(n = N_SAMPLE_L2) |>
  dplyr::ungroup()

anchor_cellids_l2 <- apply(L_pm_filtered, 2, function(x) {
  rownames(L_pm_filtered)[order(x, decreasing = TRUE)[seq_len(K_ANCHOR_L2)]]
}) |>
  as.vector() |>
  unique()

anchor_meta_l2 <- seurat_meta_filtered |>
  dplyr::filter(
    cellID %in% anchor_cellids_l2,
    annotation_level1 %in% level1_order
  )

all_meta_l2 <- dplyr::bind_rows(sampled_random_l2, anchor_meta_l2) |>
  dplyr::distinct(cellID, .keep_all = TRUE) |>
  dplyr::arrange(
    factor(annotation_level1, levels = level1_order),
    factor(annotation_level2, levels = level2_order)
  )

L_sampled_l2 <- L_pm_filtered[all_meta_l2$cellID, ]
clip_val_l2 <- quantile(L_sampled_l2, 0.99)
L_display_l2 <- pmin(L_sampled_l2, clip_val_l2)
colnames(L_display_l2) <- gsub("^K", "GP", colnames(L_display_l2))

col_fun_l2 <- colorRamp2(
  c(0, clip_val_l2 / 2, clip_val_l2),
  c("white", "#4393c3", "#08306b")
)

level2_colors <- ZemmourLib::immgent_colors$level2
level1_present_l2 <- intersect(level1_order, as.character(unique(all_meta_l2$annotation_level1)))
level2_present    <- intersect(level2_order, as.character(unique(all_meta_l2$annotation_level2)))

row_ann_l2 <- rowAnnotation(
  Cell_Type    = factor(as.character(all_meta_l2$annotation_level1), levels = level1_present_l2),
  Cell_Subtype = factor(as.character(all_meta_l2$annotation_level2), levels = level2_present),
  col = list(
    Cell_Type = level1_colors[level1_present_l2],
    Cell_Subtype = level2_colors[level2_present]
  ),
  annotation_name_gp = gpar(fontsize = 8),
  show_legend = c(Cell_Type = TRUE, Cell_Subtype = FALSE)
)

ht_l2 <- Heatmap(
  L_display_l2,
  name = "Loading",
  col = col_fun_l2,
  left_annotation = row_ann_l2,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 4),
  column_title = "Gene Programs (GPs)",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  use_raster = TRUE,
  raster_quality = 3,
  border = FALSE,
  heatmap_legend_param = list(title = "Loading", direction = "vertical")
)

pdf(
  paste0(figure_path, "heatmap_loading_overview_level1_level2.pdf"),
  width = 10,
  height = 16,
  useDingbats = FALSE
)
draw(ht_l2, merge_legend = TRUE)
dev.off()


#####################################################
#####################################################
#####################################################
### Histograms showing highly active cells per GP
#####################################################
#####################################################
options(scipen = 999)
L_pm_norm_col <- L_pm_filtered /
  matrix(
    apply(L_pm_filtered, 2, function(x) max(x)),
    nrow = nrow(L_pm_filtered),
    ncol = ncol(L_pm_filtered),
    byrow = TRUE
  )
gp_active_cell_counts <- colSums((L_pm_norm_col) > 1e-1)
pdf(
  paste0(figure_path, "hist_active_cells_per_GP.pdf"),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
hist(
  gp_active_cell_counts,
  breaks = 100,
  xlab = "Number of highly active cells per GP",
  main = "Histogram of highly active cells per GP",
  freq = T
)
dev.off()

gp_active_cell_prop <- gp_active_cell_counts / nrow(L_pm_norm_col)
pdf(
  paste0(figure_path, "hist_active_cells_prop_per_GP.pdf"),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
hist(
  gp_active_cell_prop,
  breaks = 100,
  xlab = "Proportion of highly active cells per GP",
  main = "Histogram of highly active cells per GP (proportion)",
  freq = T
)
dev.off()


#####################################################
#####################################################
#####################################################
### Boxplot showing highly active cells per GP in activated vs resting cells
#####################################################
#####################################################
cells_activated <- seurat_meta_filtered$cellID[
  seurat_meta_filtered$annotation_level2_group == "activated"
]
L_pm_activated <- L_pm_filtered[
  seurat_meta_filtered$cellID %in% cells_activated,
]
gp_active_cell_counts_activated <- rowSums(L_pm_activated > 1e-1)
cells_resting <- seurat_meta_filtered$cellID[
  seurat_meta_filtered$annotation_level2_group == "resting"
]
L_pm_resting <- L_pm_filtered[seurat_meta_filtered$cellID %in% cells_resting, ]
gp_active_cell_counts_resting <- rowSums(L_pm_resting > 1e-1)
gp_active_cell_counts_df <- data.frame(
  Group = c(
    rep("Activated", length(gp_active_cell_counts_activated)),
    rep("Resting", length(gp_active_cell_counts_resting))
  ),
  Active_Cell_Counts = c(
    gp_active_cell_counts_activated,
    gp_active_cell_counts_resting
  )
)
p <- ggplot(
  gp_active_cell_counts_df,
  aes(x = Group, y = Active_Cell_Counts, fill = Group)
) +
  geom_boxplot(outlier.size = 0.4, width = 0.6, alpha = 0.8, color = "gray40") +
  scale_fill_manual(
    values = c("Activated" = "#1f78b4", "Resting" = "#e31a1c")
  ) +
  labs(
    title = "",
    x = "",
    y = "Number of highly active GPs"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 13, face = "bold"),
    legend.position = "none"
  )
ggsave(
  filename = paste0(figure_path, "boxplot_active_cells_per_GP.pdf"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)


#####################################################
#####################################################
#####################################################
### Boxplot showing highly active cells per GP in different level 1 groups
#####################################################
#####################################################
# Build GP activity counts by annotation_level1
gp_active_cell_counts_level1 <- dplyr::bind_rows(lapply(
  level1_order,
  function(grp) {
    cells_grp <- seurat_meta_filtered$cellID[
      seurat_meta_filtered$annotation_level1 == grp
    ]
    L_grp <- L_pm_filtered[
      seurat_meta_filtered$cellID %in% cells_grp,
      ,
      drop = FALSE
    ]

    data.frame(
      Group = grp,
      Active_Cell_Counts = rowSums(L_grp > 1e-1)
    )
  }
))

group_counts <- gp_active_cell_counts_level1 %>%
  group_by(Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Group_Label = paste0(Group, "\n(n=", n, ")")) %>%
  mutate(Group = factor(Group, levels = level1_order)) %>%
  arrange(Group) %>%
  mutate(Group_Label = factor(Group_Label, levels = Group_Label))

plot_df <- gp_active_cell_counts_level1 %>%
  left_join(group_counts, by = "Group")
p <- ggplot(
  plot_df,
  aes(x = Group_Label, y = Active_Cell_Counts, fill = Group)
) +
  geom_boxplot(outlier.size = 0.4, width = 0.6, alpha = 0.8, color = "gray40") +
  scale_fill_manual(values = ZemmourLib::immgent_colors$level1) +
  labs(
    title = "Active Gene Programs per Group",
    x = "Cell Group (Annotation Level 1)",
    y = "Number of highly active GPs"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "gray30", hjust = 0.5),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 13, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )
ggsave(
  filename = paste0(figure_path, "boxplot_active_cells_per_GP_level1.pdf"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)


#####################################################
#####################################################
#####################################################
### Scatter plot showing the relationship between CD44 protein level and
### number of active GPs
#####################################################
#####################################################
gp_active_cell_counts <- rowSums(L_pm_filtered > 1e-1)
gp_active_cell_counts_df <- data.frame(
  CD44_Protein_Level = protein_mat_normalized_lognorm,
  Active_GP_Counts = gp_active_cell_counts
)
set.seed(123)
df_nz <- gp_active_cell_counts_df %>% dplyr::filter(CD44_Protein_Level > 0)
df_nz <- df_nz %>% sample_n(min(10000, nrow(df_nz)))
R <- cor(df_nz$CD44_Protein_Level, df_nz$Active_GP_Counts, use = "complete.obs")
p <- ggplot(df_nz, aes(CD44_Protein_Level, Active_GP_Counts)) +
  geom_point(alpha = 0.1, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "CD44 Protein Level vs Number of Active GPs",
    x = "CD44 Protein Level (log-normalized)",
    y = "Number of Active GPs"
  ) +
  annotate(
    "text",
    x = min(df_nz$CD44_Protein_Level, na.rm = TRUE) + 0.5,
    y = max(df_nz$Active_GP_Counts, na.rm = TRUE) - 1,
    label = paste0("R = ", round(R, 2)),
    size = 4
  ) +
  theme_minimal(base_size = 13)
ggsave(
  filename = paste0(figure_path, "scatterplot_CD44_vs_active_GPs.pdf"),
  plot = p,
  width = 6,
  height = 4,
  dpi = 300
)


#####################################################
#####################################################
#####################################################
### Histogram showing highly active genes per GP
#####################################################
#####################################################
F_pm_norm_col <- F_pm_filtered /
  matrix(
    apply(F_pm_filtered, 2, function(x) max(abs(x))),
    nrow = nrow(F_pm_filtered),
    ncol = ncol(F_pm_filtered),
    byrow = TRUE
  )
# a gene is active if its abs factor value is larger than 0.25 of the max value in that GP
gp_active_gene_counts <- colSums(abs(F_pm_norm_col) > 0.25)
pdf(
  paste0(figure_path, "hist_active_genes_per_GP.pdf"),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
hist(
  gp_active_gene_counts,
  breaks = 100,
  xlab = "Number of highly active genes per GP",
  main = "Histogram of highly active genes per GP",
  freq = T
)
dev.off()

gp_active_gene_prop <- gp_active_gene_counts / nrow(F_pm_norm_col)
pdf(
  paste0(figure_path, "hist_active_genes_prop_per_GP.pdf"),
  width = 6,
  height = 4,
  useDingbats = FALSE
)
hist(
  gp_active_gene_prop,
  breaks = 100,
  xlab = "Proportion of highly active genes per GP",
  main = "Histogram of highly active genes per GP (proportion)",
  freq = T
)
dev.off()


#####################################################
#####################################################
#####################################################
### Scatter plot showing the relationship between expected proportion of
### active cells and expected proportion of active genes per GP, colored by PVE
#####################################################
#####################################################
pve <- flashier_snmf_summary$pve
eps <- 0
pve_log <- log10(pve + eps)
pve_log_min <- quantile(pve_log, 0.02, na.rm = TRUE)
pve_log_max <- quantile(pve_log, 0.98, na.rm = TRUE)
pve_log_clipped <- pmin(pmax(pve_log, pve_log_min), pve_log_max)
pve_scaled <- (pve_log_clipped - min(pve_log_clipped)) /
  (max(pve_log_clipped) - min(pve_log_clipped))
n_col <- 100
pal <- colorRampPalette(c("lightgray", "skyblue", "darkblue"))(n_col)
idx_col <- pmin(pmax(1, floor(pve_scaled * (n_col - 1)) + 1), n_col)
col_points <- pal[idx_col]
load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
l_pi_vec <- c()
for (i in 1:200) {
  l_pi_vec <- c(l_pi_vec, flashier_snmf_fitted_prior$L_ghat[[i]]$pi[1])
}
f_pi_vec <- c()
for (i in 1:200) {
  f_pi_vec <- c(f_pi_vec, flashier_snmf_fitted_prior$F_ghat[[i]]$pi[1])
}
p_cells <- (1 - l_pi_vec)
p_genes <- (1 - f_pi_vec)
pdf(
  paste0(figure_path, "scatter_e_active_cells_vs_genes.pdf"),
  width = 6,
  height = 6
)
plot(
  p_cells,
  p_genes,
  xlab = "Expected proportion of active cells",
  ylab = "Expected proportion of active genes",
  log = "xy",
  col = col_points,
  pch = 19,
  xaxt = "n",
  yaxt = "n",
  main = "Expected proportion of active cells vs genes per GP\n(colored by PVE)"
)
ticks <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 1)
axis(1, at = ticks, labels = paste0(ticks * 100, "%"))
axis(2, at = ticks, labels = paste0(ticks * 100, "%"))
brks_pve <- c(1e-7, 1e-5, 1e-3, 1e-2, 0.1, 0.2)
brks_pve <- brks_pve[brks_pve >= min(pve) & brks_pve <= max(pve)]
brks_log <- log10(brks_pve + eps)
brks_log_clipped <- pmin(pmax(brks_log, pve_log_min), pve_log_max)
brks_scaled <- (brks_log_clipped - min(pve_log_clipped)) /
  (max(pve_log_clipped) - min(pve_log_clipped))
brks_idx <- pmin(pmax(1, floor(brks_scaled * (n_col - 1)) + 1), n_col)
brks_cols <- pal[brks_idx]
legend(
  "bottomright",
  legend = format(brks_pve, scientific = TRUE, digits = 2),
  col = brks_cols,
  pch = 15,
  pt.cex = 1.5,
  title = "PVE",
  bty = "n"
)
dev.off()
