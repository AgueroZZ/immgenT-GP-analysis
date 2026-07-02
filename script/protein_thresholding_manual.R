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
Protein_flash_result <- readRDS(file = paste0(data_path, "protein_flash_selected_summary_lognorm_backfit200.rds"))
Protein_F_pm <- Protein_flash_result$F_pm
colnames(Protein_F_pm) <- paste0("K", 1:ncol(Protein_F_pm))

good_proteins <- proteins_quality$protein[proteins_quality$classification == "good"]
good_proteins <- c(good_proteins, "IL2RA.CD25", "ITB7", "CD69")
isotype_proteins <- grep("^Isotype", rownames(Protein_F_pm), value = TRUE)
exclude_proteins <- c("CD19", "CD34", "CD45.1","CD45.2", "CD138", "TCRVA2", "TER119")
thy11_proteins <- grep("THY1.1", rownames(Protein_F_pm), value = TRUE)
select_proteins <- setdiff(good_proteins, c(exclude_proteins, thy11_proteins, isotype_proteins))


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

# load "GMM_Thresholds_Summary.csv"
threshold_results <- read.csv(paste0(data_path, "GMM_Thresholds_Summary.csv"), header = TRUE, stringsAsFactors = FALSE)
# subset for selected proteins
threshold_results_subset <- threshold_results[threshold_results$Protein %in% target_proteins,]
rownames(threshold_results_subset) <- 1:nrow(threshold_results_subset)
threshold_results_subset_manual <- c(
  "B220" = 2.75,
  "BTLA.CD272" = 3.25,
  "CD2" = 5.75,
  "CD4" = 3.85,
  "CD5" = 3.3,
  "CD11A" = 5,
  "CD24" = 4,
  "CD27" = 4,
  "CD29" = 3.5,
  "CD31" = 4.5,
  "CD38" = 4,
  "CD39" = 3,
  "CD44" = 5,
  "CD49B" = 3,
  "CD103" = 3.35,
  "CD155.PVR" = 3.25,
  "CD160" = 3,
  "CD55.DAF" = 4.25,
  "CD73.5NTD" = 5,
  "CD80" = 3,
  "CD86" = 4,
  "CD8A" = 3.95,
  "CD8B" = 4,
  "FR4" = 5,
  "GITR.CD357" = 3,
  "GR1-LY6G-LY6C1-LY6C2" = 2.5,
  "ICAM1" = 3,
  "ICOS.CD278" = 3,
  "IL7RA.CD127" = 4,
  "ITA4.CD49D" = 3.95,
  "ITAM.CD11B" = 2.25,
  "ITAX.CD11C" = 2,
  "KLRG1" = 2.85,
  "LY49A" = 3,
  "LY108" = 3,
  "CD45RB" = 7,
  "NEUROPILIN1.CD304" = 3,
  "SCA1" = 5,
  "SLAM.CD150" = 3.75,
  "TCRGD" = 3.1,
  "TCRVG2" = 2.3,
  "TCRVG3" = 3.25,
  "THY1.2" = 7,
  "IL2RA.CD25" = 2.15,
  "ITB7" = 5,
  "CD69" = 4
)
# match this to threshold_results_subset_manual to threshold_results_subset$Protein and add a new column "Threshold_manual" to threshold_results_subset
threshold_results_subset$Threshold_manual <- threshold_results_subset_manual[match(threshold_results_subset$Protein, names(threshold_results_subset_manual))]
# save threshold_results_subset as "Thresholds_Selected_Proteins.csv"
write.csv(threshold_results_subset, paste0(data_path, "Thresholds_Selected_Proteins.csv"), row.names = FALSE)



# keep only cells in cells_citeseq
L_pm_filtered <- L_pm_filtered[intersect(rownames(L_pm_filtered), cells_citeseq),]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[intersect(rownames(protein_mat_normalized_lognorm), cells_citeseq),]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[, select_proteins]


# take out thymocytes
thymocyte_cells <- seurat_meta_filtered$cellID[seurat_meta_filtered$celltype == "thymocyte"]
# take out thymocyte_cells cells from L_pm_filtered and protein_mat_normalized_lognorm
L_pm_filtered <- L_pm_filtered[!rownames(L_pm_filtered) %in% thymocyte_cells,]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[!rownames(protein_mat_normalized_lognorm) %in% thymocyte_cells,]


# 1. Define and create the output directory
output_dir <- paste0(figure_path, "protein_protein_selected/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Identify all proteins to plot against CD62L
all_proteins <- colnames(protein_mat_normalized_lognorm)
# Exclude CD62L to avoid plotting it against itself
target_proteins <- setdiff(all_proteins, "CD62L")
message(sprintf("Starting batch plotting for %d proteins...", length(target_proteins)))
# 3. Loop through each protein and generate/save the scatter plot
# Define CD62L GMM threshold for the vertical line
cd62l_gmm <- threshold_results$Threshold[threshold_results$Protein == "CD62L"][1]
for (prot in target_proteins) {

  safe_prot_name <- gsub("[^A-Za-z0-9]", "_", prot)
  file_path <- paste0(output_dir, "CD62L_vs_", safe_prot_name, ".png")

  # Extract thresholds
  gmm_val <- threshold_results_subset$Threshold[threshold_results_subset$Protein == prot][1]
  man_val <- threshold_results_subset$Threshold_manual[threshold_results_subset$Protein == prot][1]

  # 1. Generate the base plot
  p <- MyFeatureScatter_df(
    x = protein_mat_normalized_lognorm[, "CD62L"],
    y = protein_mat_normalized_lognorm[, prot],
    color_backgroud = "grey90",
    background_alpha = 0.2,
    highlight = rep(TRUE, nrow(protein_mat_normalized_lognorm)), # MUST be TRUE for density colors
    cols = grDevices::colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(10)
  ) +
    xlab("CD62L (log-normalized)") +
    ylab(paste0(prot, " (log-normalized)")) +
    ggtitle(paste0("CD62L vs ", prot)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      panel.grid = element_blank()
    )

  # 2. Add Vertical CD62L Threshold (Optional but helpful for gating context)
  if (!is.na(cd62l_gmm)) {
    p <- p + geom_vline(xintercept = cd62l_gmm, linetype = "dashed", color = "grey50", alpha = 0.6)
  }

  # 3. Add Horizontal Thresholds with Anti-Overlap Logic
  # Check if they are dangerously close (e.g., within 0.3 units)
  vjust_gmm <- -0.5
  vjust_man <- 1.5

  if (!is.na(gmm_val) && !is.na(man_val)) {
    if (abs(gmm_val - man_val) < 0.3) {
      # Push them further apart if they overlap
      vjust_gmm <- -1.2
      vjust_man <- 2.2
    }
  }

  if (!is.na(gmm_val)) {
    p <- p +
      geom_hline(yintercept = gmm_val, linetype = "dashed", color = "red", size = 0.8) +
      annotate("text", x = Inf, y = gmm_val, label = paste0("GMM: ", round(gmm_val, 2)),
               color = "red", size = 4, hjust = 1.05, vjust = vjust_gmm, fontface = "bold")
  }

  if (!is.na(man_val)) {
    p <- p +
      geom_hline(yintercept = man_val, linetype = "dotted", color = "blue", size = 0.8) +
      annotate("text", x = Inf, y = man_val, label = paste0("Manual: ", round(man_val, 2)),
               color = "blue", size = 4, hjust = 1.05, vjust = vjust_man, fontface = "bold")
  }

  ggsave(filename = file_path, plot = p, width = 8, height = 7, dpi = 300)
}
message("Batch plotting complete with dual thresholds. Figures saved in: ", output_dir)




