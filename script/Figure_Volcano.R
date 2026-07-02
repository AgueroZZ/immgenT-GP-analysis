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
data_path <- "data/"
code_path <- "code/"
figure_path <- "figures/"
source(paste0(code_path, "ROC.R"))
source(paste0(code_path, "filtering_membership.R"))

#####################################################
#####################################################
#####################################################
### Loading Data
#####################################################
#####################################################
load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
# L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
# F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered_500.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered_500.rds"))
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered),]
umap_result <- readRDS(paste0(data_path, "umap_result.rds"))
colnames(umap_result) <- c("UMAP_1", "UMAP_2")
umap_result <- umap_result[rownames(L_pm_filtered),]
df_umap <- as.data.frame(umap_result)
mean_shifted_log_expr <- readRDS(paste0(data_path, "mean_shifted_log_expr.rds"))




#####################################################
#####################################################
#####################################################
### Defining functions for plotting
#####################################################
#####################################################
plot_volcano_gp <- function(factor_index,
                            factor_matrix,
                            mean_expr_vector,
                            n_highlight = 10,
                            normalize = "max_abs") {

  # 1. Align genes
  gene_names <- rownames(factor_matrix)
  common_genes <- intersect(gene_names, names(mean_expr_vector))

  if(length(common_genes) == 0) {
    stop("Error: No matching gene names found between factor_matrix and mean_expr_vector.")
  }

  # 2. Extract values
  f_values <- factor_matrix[common_genes, factor_index]

  # Logic for Normalization
  if (!is.null(normalize) && normalize == "max_abs") {
    norm_val <- max(abs(f_values), na.rm = TRUE)
    # Avoid division by zero if the factor is empty/constant
    if(norm_val != 0) {
      f_values <- f_values / norm_val
    }
  }

  df <- data.frame(
    Gene = common_genes,
    FactorValue = f_values,
    MeanLogExpr = mean_expr_vector[common_genes]
  )

  # 3. Identify Top Genes
  top_genes <- df %>%
    slice_max(order_by = abs(FactorValue), n = n_highlight)

  # 4. Construct Plot
  p <- ggplot(df, aes(x = FactorValue, y = MeanLogExpr)) +
    # RASTER LAYER
    scattermore::geom_scattermore(
      pointsize = 2,
      alpha = 0.2,
      color = "grey75"
    ) +
    # VECTOR LAYERS (Illustrator friendly)
    geom_point(data = top_genes, color = "dodgerblue", size = 1.5) +
    ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = Gene),
      color = "black",
      size = 3,
      segment.color = "grey50",
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = Inf,
      min.segment.length = 0
    ) +
    labs(
      title = paste0("Factor: ", colnames(factor_matrix)[factor_index]),
      x = "Increase in Shifted Log Expression (Normalized)",
      y = "Mean Shifted Log Expression"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  return(p)
}




#####################################################
#####################################################
#####################################################
### Plotting Volcano Plot for GP
#####################################################
#####################################################
# Parameters
selected_factors <- 1:ncol(F_pm_filtered)
output_file <- paste0(figure_path, "volcano_GPs_backfitting.pdf")
output_file <- paste0(figure_path, "volcano_GPs_backfitting_500.pdf")
# Open PDF
pdf(output_file, width = 7, height = 7)
# Loop through and plot
for (i in selected_factors) {
  # Generate the plot
  # Passing F_pm_filtered directly as the factor_matrix
  p <- plot_volcano_gp(
    factor_index = i,
    factor_matrix = F_pm_filtered,
    mean_expr_vector = mean_shifted_log_expr,
    n_highlight = 10
  )

  # Explicitly print to the PDF device
  print(p)

  # Optional: Progress tracker in console
  if(i %% 50 == 0) message("Finished factor ", i)
}
dev.off()



