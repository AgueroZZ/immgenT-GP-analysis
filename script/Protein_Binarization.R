library(ggplot2)
library(dplyr)
library(tidyr)
data_path <- "data/"
figure_path <- "figures/"
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered),]
protein_mat_normalized_lognorm <- readRDS(paste0(data_path, "protein_mat_normalized_lognorm.rds"))
proteins_quality <- read.csv(paste0(data_path, "TableS4_citeseq_qc_20250513.csv"), header = TRUE, stringsAsFactors = FALSE, skip = 1)
cells_citeseq <- seurat_meta_filtered$cellID[seurat_meta_filtered$cite_seq]




MyFeatureScatter_df <- function(
    x,
    y,
    highlight,
    split = NULL,
    feature1 = "feature1",
    feature2 = "feature2",
    raster = TRUE,
    color_backgroud = "grey",
    cols = rev(rainbow(10, end = 4/6)),
    highlight_size = 1,
    highlight_alpha = 1,
    background_alpha = 1,
    base_pixels = c(512, 512),
    highlight_pixels = c(216, 216),
    nbin = 500
) {
  requireNamespace("scattermore", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)

  # ---- checks ----
  if (missing(x) || missing(y) || missing(highlight)) {
    stop("Please provide x, y, and highlight.")
  }
  if (length(x) != length(y) || length(x) != length(highlight)) {
    stop("x, y, and highlight must have the same length.")
  }
  if (!is.null(split) && length(split) != length(x)) {
    stop("split must have the same length as x/y.")
  }

  # coerce highlight to logical safely
  highlight <- as.logical(highlight)

  df <- data.frame(
    feature1 = as.numeric(x),
    feature2 = as.numeric(y),
    highlight = highlight,
    stringsAsFactors = FALSE
  )
  if (!is.null(split)) df$split <- split

  # drop rows with NA in x/y (and split if present)
  keep <- !is.na(df$feature1) & !is.na(df$feature2)
  if (!is.null(split)) keep <- keep & !is.na(df$split)
  df <- df[keep, , drop = FALSE]

  df2 <- df[!is.na(df$highlight) & df$highlight, , drop = FALSE]

  # If nothing to highlight, still return the base plot
  p1 <- ggplot2::ggplot(df) +
    scattermore::geom_scattermore(
      ggplot2::aes(feature1, feature2),
      color = color_backgroud,
      alpha = background_alpha,
      pixels = base_pixels
    )

  if (nrow(df2) > 0) {
    df2$density_col <- grDevices::densCols(
      df2$feature1, df2$feature2,
      colramp = grDevices::colorRampPalette(cols),
      nbin = nbin
    )

    if (isTRUE(raster)) {
      p2 <- scattermore::geom_scattermore(
        data = df2,
        ggplot2::aes(feature1, feature2, color = density_col),
        pointsize = highlight_size,
        pixels = highlight_pixels
      )
    } else {
      p2 <- ggplot2::geom_point(
        data = df2,
        ggplot2::aes(feature1, feature2, color = density_col),
        size = highlight_size,
        alpha = highlight_alpha
      )
    }
  } else {
    p2 <- NULL
  }

  p <- p1 + p2 +
    ggplot2::xlab(feature1) + ggplot2::ylab(feature2) +
    ggplot2::scale_colour_identity() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 10),
      axis.title.x = ggplot2::element_text(size = 20),
      axis.title.y = ggplot2::element_text(size = 20),
      legend.title = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(split)) {
    p <- p + ggplot2::facet_wrap(~split)
  }

  return(p)
}


# keep only cells in cells_citeseq
L_pm_filtered <- L_pm_filtered[intersect(rownames(L_pm_filtered), cells_citeseq),]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[intersect(rownames(protein_mat_normalized_lognorm), cells_citeseq),]

# remove proteins with poor quality
poor_proteins <- proteins_quality$protein[proteins_quality$classification == "poor"]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[, !colnames(protein_mat_normalized_lognorm) %in% poor_proteins]


# take out thymocytes
thymocyte_cells <- seurat_meta_filtered$cellID[seurat_meta_filtered$celltype == "thymocyte"]
# take out thymocyte_cells cells from L_pm_filtered and protein_mat_normalized_lognorm
L_pm_filtered <- L_pm_filtered[!rownames(L_pm_filtered) %in% thymocyte_cells,]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[!rownames(protein_mat_normalized_lognorm) %in% thymocyte_cells,]







# 1. Define and create the output directory
output_dir <- paste0(figure_path, "protein_protein/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Identify all proteins to plot against CD62L
all_proteins <- colnames(protein_mat_normalized_lognorm)
# Exclude CD62L to avoid plotting it against itself
target_proteins <- setdiff(all_proteins, "CD62L")

message(sprintf("Starting batch plotting for %d proteins...", length(target_proteins)))

# 3. Loop through each protein and generate/save the scatter plot
for (prot in target_proteins) {

  # Create a filesystem-safe filename by replacing special characters with underscores
  safe_prot_name <- gsub("[^A-Za-z0-9]", "_", prot)
  file_path <- paste0(output_dir, "CD62L_vs_", safe_prot_name, ".png")

  # Generate the plot using your custom function
  p <- MyFeatureScatter_df(
    x = protein_mat_normalized_lognorm[, "CD62L"],
    y = protein_mat_normalized_lognorm[, prot],
    color_backgroud = "black",
    background_alpha = 0.5,
    # No highlighting, just show density
    highlight = rep(FALSE, nrow(protein_mat_normalized_lognorm)),
    cols = grDevices::colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(10)
  ) +
    xlab("CD62L (log-normalized)") +
    ylab(paste0(prot, " (log-normalized)")) +
    ggtitle(paste0("Protein Expression: CD62L vs ", prot)) +
    theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      legend.text = element_text(size = 10),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  # Save the plot
  # Adjust width/height as needed for your research report
  ggsave(filename = file_path, plot = p, width = 8, height = 7, dpi = 300)
}

message("Batch plotting complete. Figures saved in: ", output_dir)






library(mclust)
library(ggplot2)
library(dplyr)

# 1. Define and create the output directory
output_dir <- paste0(figure_path, "protein_density/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Identify all proteins to process
all_proteins <- colnames(protein_mat_normalized_lognorm)

# Initialize threshold results
threshold_results <- data.frame(Protein = character(), Threshold = numeric(), stringsAsFactors = FALSE)

message("Starting GMM-based thresholding for counts > 0.5...")

for (prot in all_proteins) {

  # Get log-normalized values
  prot_vals <- protein_mat_normalized_lognorm[, prot]

  # SUBSET: Only values > 0.5 for GMM fitting
  non_zero_vals <- prot_vals[prot_vals > 0.5]

  if (length(non_zero_vals) < 50) {
    message(sprintf("Skipping %s: too few entries > 0.5.", prot))
    next
  }

  # 3. Apply Gaussian Mixture Modeling (GMM)
  set.seed(42)
  gmm_fit <- Mclust(non_zero_vals, G = 2, verbose = FALSE)

  # Logic: threshold is the boundary of the 'negative' cluster within the subset
  means <- gmm_fit$parameters$mean
  neg_cluster <- which.min(means)
  threshold <- max(non_zero_vals[gmm_fit$classification == neg_cluster], na.rm = TRUE)

  threshold_results <- rbind(threshold_results, data.frame(Protein = prot, Threshold = threshold))

  # 4. Create Histogram + Density Overlay
  safe_prot_name <- gsub("[^A-Za-z0-9]", "_", prot)

  # We create a dataframe of ONLY the values > 0.5 to match your GMM logic
  plot_df <- data.frame(val = non_zero_vals)

  p <- ggplot(plot_df, aes(x = val)) +
    # Histogram using density scale on Y axis
    geom_histogram(aes(y = after_stat(density)), bins = 50,
                   fill = "grey90", color = "grey60", alpha = 0.8) +
    # Density curve overlay
    geom_density(fill = "steelblue", alpha = 0.3, color = "steelblue", linewidth = 1) +
    # Threshold line
    geom_vline(xintercept = threshold, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = threshold, y = Inf, label = paste("Threshold =", round(threshold, 3)),
             vjust = 2, hjust = -0.1, color = "red", fontface = "bold") +
    labs(
      title = paste("GMM Thresholding:", prot),
      subtitle = "Histogram & Density (Filtered for values > 0.5)",
      x = paste(prot, "(log-normalized)"),
      y = "Density"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )

  # 5. Save
  ggsave(filename = paste0(output_dir, "Density_Hist_", safe_prot_name, ".png"),
         plot = p, width = 7, height = 5, dpi = 300)
}

write.csv(threshold_results, paste0(output_dir, "GMM_Thresholds_Summary.csv"), row.names = FALSE)
message("GMM Thresholding complete.")
