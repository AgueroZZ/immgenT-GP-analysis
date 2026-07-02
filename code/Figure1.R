library(ggplot2)
library(dplyr)
library(fastTopics)
library(qs)
library(cowplot)
library(ggrepel)
library(ZemmourLib)

data_path <- "data/"
code_path <- "code/"
figure_path <- "figures/"
seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))

### Alignment of cells
cells_flashier <- rownames(flashier_snmf_summary$L_pm)
cells_seurat <- seurat_meta$cellID
which(!cells_flashier %in% cells_seurat)
L_pm <- flashier_snmf_summary$L_pm[cells_flashier %in% cells_seurat, ]
all(rownames(L_pm) == cells_seurat)

### Filter cells with high total membership
filter_cells_by_total_membership <- function (L, max_val = 10, numiter = 10) {
  n <- nrow(L)
  rows <- 1:n
  for (iter in 1:numiter) {
    x <- rowSums(L)
    cat(sprintf("%d. Filtered out %d cells.\n",iter,sum(x > max_val)))
    i <- which(x <= max_val)
    L <- L[i,]
    rows <- rows[i]
    d <- apply(L,2,max)
    L <- scale_cols(L,1/d)
  }
  return(rows)
}
scale_cols <- function (A, b) t(t(A) * b)
D <- diag(1 / apply(L_pm, 2, function(x) max(x)))
L <- L_pm %*% D
cells <- filter_cells_by_total_membership(L,numiter = 12)
seurat_meta <- seurat_meta[cells,]
L_pm <- L_pm[cells,]

### Normalizing loadings
d <- apply(L_pm,2,max)
L_pm <- scale_cols(L_pm,1/d)
D <- (1/d) * D

### Remove unhealthy cells and thymocytes
cells_thymocyte <- which(seurat_meta$annotation_level1 == "thymocyte")
cells_unhealthy_thymocyte <- which(seurat_meta$annotation_level1 == "thymocyte" | seurat_meta$condition_broad != "healthy")
L_pm_no_unhealthy_thymocytes <- L_pm[-cells_unhealthy_thymocyte, ]


### Highlighted cell types
set.seed(12345)
### These are all factors that passed the validation test.
highlights <- c("K3","K22","K30","K58","K68","K171")
color_coding <- ZemmourLib::immgent_colors$level1
Level1_match <- c("CD8aa", "DN", "nonconv", "CD8", "Treg", "gdT")
cell_type_factors <- as.numeric(gsub("K", highlights, replacement = ""))
fit2 <- L_pm_no_unhealthy_thymocytes[, cell_type_factors, drop = FALSE]
colnames(fit2) <- paste0("k", cell_type_factors)

cell_type <- seurat_meta$annotation_level1[-cells_unhealthy_thymocyte]
cells <- which(cell_type == "CD4" | cell_type == "CD8")
cells <- sample(cells,5e4)
cells <- sort(c(cells,which(cell_type != "CD4" & cell_type != "CD8")))



### Figure 1(a): Structure plot
structure_plot(fit2[cells,],gap = 20,n = 3000,
               colors = color_coding[Level1_match],
               grouping = cell_type[cells]) +
  labs(y = "membership",color = "",fill = "") +
  guides(
    fill = guide_legend(nrow = 1),
    color = guide_legend(nrow = 1)
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )
ggsave(paste0(figure_path, "Figure1a_structure_plot.pdf"),
       width = 8, height = 5)


### Figure 1(b): Gene Tile plot
F_pm <- flashier_snmf_summary$F_pm
DF <- apply(F_pm,2,max)
F_pm <- F_pm %*% solve(diag(DF))
colnames(F_pm) <- paste0("K", 1:ncol(F_pm))
p <- MyGeneTilePlot(F_pm,
                    highlights,
                    alpha_range = NULL,
                    n_genes_per_topic = 10) +
  theme(aspect.ratio = 1.5) +
  # Remove legend
  guides(alpha = "none") +
  labs(
  title = "Top Expressed Genes per Factor",
  fill = "Relative Expression")
p
ggsave(paste0(figure_path, "Figure1b_gene_tile_plot.pdf"),
       plot = p,
       width = 6, height = 9)
# annotation_heatmap(F_pm[,highlights],n = 2,
#                    dims = highlights,
#                    compare_dims = highlights,
#                          select_features = "distinctive",
#                          font_size = 9) +
#   theme(plot.title = element_text(face = "plain",size = 9))




### Figure 1(c): Protein Annotation Plot
protein_flash_summary_lognorm <- readRDS(file = paste0(data_path, "protein_flash_selected_summary_lognorm.rds"))
# Normalize protein loadings and factors for visualization
D_lognorm <- diag(1 / apply(protein_flash_summary_lognorm$F_pm, 2, function(x) (max(x) - min(x))))
# D_lognorm <- diag(1 / apply(protein_flash_summary_lognorm$F_pm, 2, function(x) max(x)))
protein_flash_summary_lognorm$F_pm <- protein_flash_summary_lognorm$F_pm %*% D_lognorm
protein_flash_summary_lognorm$L_pm <- protein_flash_summary_lognorm$L_pm %*% solve(D_lognorm)

F_pm <- protein_flash_summary_lognorm$F_pm
colnames(F_pm) <- paste0("F", 1:ncol(flashier_snmf_summary$F_pm))[ flashier_snmf_summary$pve > 1e-4]

# kset <- paste0("F", c(3, 22, 30, 58, 68, 171))
kset <- paste0("F", c(68, 58, 30, 22, 3, 171))
F_pm <- F_pm[, kset]
annotation_heatmap(F_pm,n = 3,
                   dim = kset,
                   show_dims = kset,
                   select_features = "largest",
                   feature_sign = "both",
                   # zero_value = 0.1,
                   compare_dims = kset,
                   font_size = 15) +
  theme(plot.title = element_text(face = "plain",size = 10))
ggsave(paste0(figure_path, "Figure1c_protein_annotation_plot.pdf"),
       width = 6, height = 5)

