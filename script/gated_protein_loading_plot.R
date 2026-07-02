library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

data_path <- "data/"
figure_path <- "figures/"
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta <- readRDS(paste0(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
protein_mat_normalized_lognorm <- readRDS(paste0(
  data_path,
  "protein_mat_normalized_lognorm.rds"
))
if (isS4(protein_mat_normalized_lognorm)) {
  protein_mat_normalized_lognorm <- as.matrix(protein_mat_normalized_lognorm)
}
proteins_quality <- read.csv(
  paste0(data_path, "TableS4_citeseq_qc_20250513.csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  skip = 1
)
cells_citeseq <- seurat_meta_filtered$cellID[seurat_meta_filtered$cite_seq]
Protein_flash_result <- readRDS(
  file = paste0(
    data_path,
    "protein_flash_selected_summary_lognorm_backfit200.rds"
  )
)
Protein_F_pm <- Protein_flash_result$F_pm
colnames(Protein_F_pm) <- paste0("K", 1:ncol(Protein_F_pm))

# Keep only cells in cells_citeseq
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[
  rownames(protein_mat_normalized_lognorm) %in% cells_citeseq,
]
seurat_meta_filtered <- seurat_meta_filtered %>%
  filter(cellID %in% cells_citeseq)
L_pm_filtered <- L_pm_filtered[rownames(L_pm_filtered) %in% cells_citeseq, ]


good_proteins <- proteins_quality$protein[
  proteins_quality$classification == "good"
]
good_proteins <- c(good_proteins, "IL2RA.CD25", "ITB7", "CD69")
isotype_proteins <- grep("^Isotype", rownames(Protein_F_pm), value = TRUE)
exclude_proteins <- c(
  "CD19",
  "CD34",
  "CD45.1",
  "CD45.2",
  "CD138",
  "TCRVA2",
  "TER119"
)
thy11_proteins <- grep("THY1.1", rownames(Protein_F_pm), value = TRUE)
select_proteins <- setdiff(
  good_proteins,
  c(exclude_proteins, thy11_proteins, isotype_proteins)
)

# threshold_results <- read.csv(paste0(data_path, "GMM_Thresholds_Summary.csv"), header = TRUE, stringsAsFactors = FALSE)
threshold_results_subset <- read.csv(
  paste0(data_path, "Thresholds_Selected_Proteins.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)
threshold_results_subset <- threshold_results_subset[, c(
  "Protein",
  "Threshold",
  "Threshold_manual"
)]
threshold_results_subset <- rbind(
  threshold_results_subset,
  data.frame(Protein = "CD62L", Threshold = 3, Threshold_manual = 3)
)
threshold_results_subset_manual <- data.frame(
  Protein = threshold_results_subset$Protein,
  Threshold = threshold_results_subset$Threshold_manual
)
df_markers <- readRDS(paste0(data_path, "CITEseq_markers_full.rds"))
umap_result <- readRDS(paste0(data_path, "umap_result.rds"))
colnames(umap_result) <- c("UMAP_1", "UMAP_2")
umap_result <- umap_result[rownames(L_pm_filtered), ]


MyDimPlotHighlightDensity_df <- function(
  emb, # n x 2 matrix/data.frame: columns are dim1, dim2
  highlight, # length n logical (or coercible): which cells to highlight
  split = NULL, # optional length n vector for facet
  dim_names = c("Dim1", "Dim2"),
  raster = TRUE,
  highlight_size = 0.5,
  highlight_alpha = 0.5,
  highlight_pointsize = 0L, # scattermore pointsize for raster highlight layer
  base_pixels = c(512, 512),
  highlight_pixels = c(512, 512),
  cols = rev(rainbow(10, end = 4 / 6)),
  nbin = 500
) {
  requireNamespace("scattermore", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)

  # ---- checks ----
  if (missing(emb) || missing(highlight)) {
    stop("Please provide emb (n x 2) and highlight.")
  }

  emb <- as.data.frame(emb)
  if (ncol(emb) < 2) {
    stop("emb must have at least 2 columns (dim1, dim2).")
  }
  emb <- emb[, 1:2, drop = FALSE]

  n <- nrow(emb)
  if (length(highlight) != n) {
    stop("highlight must have length nrow(emb).")
  }
  if (!is.null(split) && length(split) != n) {
    stop("split must have length nrow(emb).")
  }

  highlight <- as.logical(highlight)

  df <- data.frame(
    feature1 = as.numeric(emb[[1]]),
    feature2 = as.numeric(emb[[2]]),
    highlight = highlight,
    stringsAsFactors = FALSE
  )
  if (!is.null(split)) {
    df$split <- split
  }

  # drop NA coords (and split if present)
  keep <- !is.na(df$feature1) & !is.na(df$feature2)
  if (!is.null(split)) {
    keep <- keep & !is.na(df$split)
  }
  df <- df[keep, , drop = FALSE]

  df2 <- df[!is.na(df$highlight) & df$highlight, , drop = FALSE]

  # ---- base plot (all cells in grey) ----
  p1 <- ggplot2::ggplot(df) +
    scattermore::geom_scattermore(
      ggplot2::aes(feature1, feature2),
      color = "grey",
      pixels = base_pixels
    )

  # ---- highlight layer (density-colored) ----
  if (nrow(df2) > 0) {
    # densCols needs enough unique points to form bin breaks; fall back to a
    # solid colour (top of the ramp) when there are too few highlighted cells.
    df2$density_col <- tryCatch(
      grDevices::densCols(
        df2$feature1,
        df2$feature2,
        colramp = grDevices::colorRampPalette(cols),
        nbin = nbin
      ),
      error = function(e) rep(tail(cols, 1), nrow(df2))
    )

    if (isTRUE(raster)) {
      p2 <- scattermore::geom_scattermore(
        data = df2,
        ggplot2::aes(feature1, feature2, color = density_col),
        pixels = highlight_pixels,
        pointsize = highlight_pointsize
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

  p <- p1 +
    p2 +
    ggplot2::xlab(dim_names[1]) +
    ggplot2::ylab(dim_names[2]) +
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

gp_label <- function(x) sub("^K", "GP", x)

plot_gated_gp_vs_protein <- function(
  gp_name,
  df_markers,
  protein_mat,
  loading_mat,
  umap_emb,
  threshold_df,
  selected_proteins = NULL,
  exclude_cells = NULL,
  loading_q = 0.9,
  min_cells = 2,
  missing_threshold_action = "median", # NEW: "median" or "skip"
  min_pointsize = 0L, # floor on raster highlight pointsize (overrides count-based scaling)
  save_path = NULL
) {
  # 1. Cell Exclusion Logic
  if (!is.null(exclude_cells)) {
    cells_to_keep <- setdiff(rownames(loading_mat), exclude_cells)
    loading_mat <- loading_mat[cells_to_keep, , drop = FALSE]
    common <- intersect(cells_to_keep, rownames(protein_mat))
    protein_mat <- protein_mat[common, , drop = FALSE]
    common_umap <- intersect(common, rownames(umap_emb))
    umap_emb <- umap_emb[common_umap, , drop = FALSE]
  }

  # 2. Extract markers (df_markers is keyed by "GP..."; gp_name is "K...")
  row_idx <- which(rownames(df_markers) == gp_label(gp_name))
  if (length(row_idx) == 0) {
    stop("GP name not found in df_markers.")
  }

  pos_markers <- strsplit(df_markers[row_idx, "Positive"], ", ")[[1]]
  neg_markers <- strsplit(df_markers[row_idx, "Negative"], ", ")[[1]]

  if (!is.null(selected_proteins)) {
    pos_markers <- intersect(pos_markers, selected_proteins)
    neg_markers <- intersect(neg_markers, selected_proteins)
  }
  pos_markers <- pos_markers[pos_markers != ""]
  neg_markers <- neg_markers[neg_markers != ""]

  # 3. Protein Gating Logic
  is_protein_gated <- rep(TRUE, nrow(protein_mat))
  names(is_protein_gated) <- rownames(protein_mat)

  # Internal helper to handle the logic for each marker
  apply_gate <- function(marker_list, is_positive = TRUE) {
    for (p in marker_list) {
      if (p %in% colnames(protein_mat)) {
        # Check if threshold exists
        if (p %in% threshold_df$Protein) {
          thresh <- threshold_df$Threshold[threshold_df$Protein == p]
          if (is_positive) {
            is_protein_gated <<- is_protein_gated & (protein_mat[, p] > thresh)
          } else {
            is_protein_gated <<- is_protein_gated & (protein_mat[, p] <= thresh)
          }
        } else {
          # Handle missing threshold
          if (missing_threshold_action == "skip") {
            message(sprintf(
              "[%s] Warning: %s threshold missing. Skipping this marker.",
              gp_name,
              p
            ))
          } else {
            thresh <- median(protein_mat[, p], na.rm = TRUE)
            message(sprintf(
              "[%s] Warning: %s threshold missing. Using median (%.2f).",
              gp_name,
              p,
              thresh
            ))
            if (is_positive) {
              is_protein_gated <<- is_protein_gated &
                (protein_mat[, p] > thresh)
            } else {
              is_protein_gated <<- is_protein_gated &
                (protein_mat[, p] <= thresh)
            }
          }
        }
      }
    }
  }

  apply_gate(pos_markers, is_positive = TRUE)
  apply_gate(neg_markers, is_positive = FALSE)

  n_prot <- sum(is_protein_gated)

  # 4. GP Loading Gating Logic
  loadings <- loading_mat[, gp_name]
  if (is.null(loading_q)) {
    if (n_prot <= 1) {
      loading_q_val <- 0.999
    } else {
      loading_q_val <- 1 - (n_prot / length(loadings))
      loading_q_val <- max(0, min(0.9999, loading_q_val))
    }
    q_label <- paste0("Matched n=", n_prot)
  } else {
    loading_q_val <- loading_q
    q_label <- paste0("Top ", round((1 - loading_q_val) * 100), "%")
  }

  loading_cutoff <- quantile(loadings, loading_q_val)
  is_loading_gated <- loadings >= loading_cutoff
  n_load <- sum(is_loading_gated)

  message(sprintf(
    "[%s] Final Gate: Protein=%d, Loading=%d",
    gp_name,
    n_prot,
    n_load
  ))

  # 5. Plot Preparation (Same as before)
  common_cells <- intersect(rownames(loading_mat), rownames(umap_emb))
  emb_subset <- umap_emb[common_cells, ]
  pos_subtitle <- if (length(pos_markers) > 0) {
    paste0(paste(pos_markers, collapse = "+ "), "+")
  } else {
    "None+"
  }
  neg_subtitle <- if (length(neg_markers) > 0) {
    paste0(paste(neg_markers, collapse = "- "), "-")
  } else {
    "None-"
  }

  render_plot <- function(highlight_vec, title_text, cell_count) {
    if (cell_count < min_cells) {
      label_text <- if (cell_count == 0) {
        "Zero cells in gate"
      } else {
        paste0("Too few cells (n=", cell_count, ")")
      }
      return(
        ggplot() +
          annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = label_text,
            size = 5,
            fontface = "italic"
          ) +
          ggtitle(title_text) +
          theme_void() +
          theme(plot.title = element_text(size = 9, face = "bold"))
      )
    }

    # Scale raster pointsize by number of gated cells; background stays rasterized
    hl_pointsize <- if (cell_count < 100) {
      8L
    } else if (cell_count < 1000) {
      3L
    } else {
      0L
    }
    hl_pointsize <- max(hl_pointsize, min_pointsize)

    MyDimPlotHighlightDensity_df(
      emb = emb_subset,
      highlight = highlight_vec[common_cells],
      dim_names = c("MDE1", "MDE2"),
      cols = grDevices::colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(
        10
      ),
      highlight_pointsize = hl_pointsize
    ) +
      ggtitle(title_text) +
      theme(
        plot.title = element_text(size = 9, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()
      )
  }

  gp_disp <- gp_label(gp_name)
  p1 <- render_plot(
    is_protein_gated,
    paste0(gp_disp, " (", q_label, ")\n", pos_subtitle, "\n", neg_subtitle),
    n_prot
  )
  p2 <- render_plot(
    is_loading_gated,
    "",
    n_load
  )

  combined_plot <- p1 + p2 + plot_layout(ncol = 2)

  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = combined_plot, width = 12, height = 6)
  }

  return(combined_plot)
}

get_gated_cell_ids <- function(
  gp_name,
  df_markers,
  protein_mat,
  threshold_df,
  selected_proteins = NULL,
  exclude_cells = NULL,
  missing_threshold_action = "skip"
) {
  # 1. Filter cells initially
  all_cells <- rownames(protein_mat)
  if (!is.null(exclude_cells)) {
    all_cells <- setdiff(all_cells, exclude_cells)
  }

  # 2. Extract markers (df_markers is keyed by "GP..."; gp_name is "K...")
  row_idx <- which(rownames(df_markers) == gp_label(gp_name))
  if (length(row_idx) == 0) {
    warning(sprintf("GP %s not found. Returning empty vector.", gp_name))
    return(character(0))
  }

  pos_markers <- strsplit(df_markers[row_idx, "Positive"], ", ")[[1]]
  neg_markers <- strsplit(df_markers[row_idx, "Negative"], ", ")[[1]]

  if (!is.null(selected_proteins)) {
    pos_markers <- intersect(pos_markers, selected_proteins)
    neg_markers <- intersect(neg_markers, selected_proteins)
  }
  pos_markers <- pos_markers[pos_markers != ""]
  neg_markers <- neg_markers[neg_markers != ""]

  # 3. Gating Logic
  # Start with a logical vector for the subset of cells we care about
  sub_mat <- protein_mat[all_cells, , drop = FALSE]
  keep_vec <- rep(TRUE, nrow(sub_mat))
  names(keep_vec) <- all_cells

  process_markers <- function(markers, is_pos) {
    for (p in markers) {
      if (p %in% colnames(sub_mat)) {
        thresh <- NULL
        if (p %in% threshold_df$Protein) {
          thresh <- threshold_df$Threshold[threshold_df$Protein == p]
        } else if (missing_threshold_action == "median") {
          thresh <- median(sub_mat[, p], na.rm = TRUE)
        }

        if (!is.null(thresh)) {
          if (is_pos) {
            keep_vec <<- keep_vec & (sub_mat[, p] > thresh)
          } else {
            keep_vec <<- keep_vec & (sub_mat[, p] <= thresh)
          }
        }
      }
    }
  }

  process_markers(pos_markers, is_pos = TRUE)
  process_markers(neg_markers, is_pos = FALSE)

  # Return the names (Cell IDs) that are TRUE
  return(names(keep_vec)[keep_vec])
}

thymocyte_cells <- seurat_meta_filtered %>%
  filter(annotation_level1 == "thymocyte") %>%
  pull(cellID)
proliferating_cells <- seurat_meta_filtered %>%
  filter(annotation_level2_group == "proliferating") %>%
  pull(cellID)
miniverse_cells <- seurat_meta_filtered %>%
  filter(annotation_level2_group == "miniverse") %>%
  pull(cellID)


# ### Generate gated plot for every GP (using default markers from df_markers)
# all_gp_names <- paste0("K", seq_len(nrow(df_markers)))
# pdf_path_all <- paste0(
#   figure_path,
#   "gated_umap/raw_markers/gated_umap_all_GPs.pdf"
# )
# if (!dir.exists(dirname(pdf_path_all))) {
#   dir.create(dirname(pdf_path_all), recursive = TRUE)
# }

# message("Starting PDF: all GPs (", length(all_gp_names), " plots)...")
# cairo_pdf(pdf_path_all, width = 12, height = 6, onefile = TRUE, bg = "white")
# for (gp in all_gp_names) {
#   message("Plotting ", gp, "...")
#   tryCatch(
#     {
#       p <- plot_gated_gp_vs_protein(
#         gp_name = gp,
#         df_markers = df_markers,
#         protein_mat = protein_mat_normalized_lognorm,
#         loading_mat = L_pm_filtered,
#         umap_emb = umap_result,
#         threshold_df = threshold_results_subset_manual,
#         missing_threshold_action = "skip",
#         exclude_cells = c(
#           thymocyte_cells,
#           proliferating_cells,
#           miniverse_cells
#         ),
#         selected_proteins = select_proteins,
#         loading_q = NULL
#       )
#       print(p)
#     },
#     error = function(e) {
#       message("Error plotting ", gp, ": ", e$message)
#     }
#   )
# }
# dev.off()
# message("Done. Saved: ", pdf_path_all)

#### Manually curate markers for specific GPs and generate individual plots for them

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

df_markers2$Positive[
  41
] <- "CD55.DAF, CD4, CD44, CD45RB, CD62L, CD31, CD2, IL7RA.CD127, CD27, SCA1, CD5"


df_markers2$Positive[57] <- c("CD44, ICOS.CD278")
df_markers2$Negative[57] <- c("CD62L")

df_markers2$Positive[68] <- "IL2RA.CD25, FR4, GITR.CD357, NEUROPILIN1.CD304"
df_markers2$Negative[68] <- ""

df_markers2$Positive[80] <- "CD29, CD44, ITA4.CD49D"
df_markers2$Negative[80] <- "ITB7, CD103"

df_markers2$Positive[170] <- "ITB7, CD103, CD4, CD38"
df_markers2$Negative[170] <- "GITR.CD357, CD62L"

df_markers2$Positive[171] <- "CD62L"
df_markers2$Negative[171] <- "CD44"


# #############################################
# #############################################
# #############################################
# #############################################
# p1 <- plot_gated_gp_vs_protein(
#   gp_name = "K3",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p1,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP3_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# plot_gated_gp_vs_protein(
#   gp_name = "K41",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# p1 <- plot_gated_gp_vs_protein(
#   gp_name = "K6",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p1,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP6_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# p1 <- plot_gated_gp_vs_protein(
#   gp_name = "K8",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p1,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP8_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# p2 <- plot_gated_gp_vs_protein(
#   gp_name = "K22",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p2,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP22_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# p3 <- plot_gated_gp_vs_protein(
#   gp_name = "K23",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p3,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP23_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# p4 <- plot_gated_gp_vs_protein(
#   gp_name = "K26",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p4,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP26_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# p5 <- plot_gated_gp_vs_protein(
#   gp_name = "K27",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p5,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP27_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# p6 <- plot_gated_gp_vs_protein(
#   gp_name = "K29",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p6,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP29_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# p7 <- plot_gated_gp_vs_protein(
#   gp_name = "K30",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p7,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP30_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# p8 <- plot_gated_gp_vs_protein(
#   gp_name = "K57",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p8,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP57_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# p8 <- plot_gated_gp_vs_protein(
#   gp_name = "K58",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p8,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP58_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# p9 <- plot_gated_gp_vs_protein(
#   gp_name = "K68",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p9,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP68_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# p10 <- plot_gated_gp_vs_protein(
#   gp_name = "K80",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p10,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP80_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# p11 <- plot_gated_gp_vs_protein(
#   gp_name = "K170",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p11,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP170_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

# high_cells_k170 <- get_gated_cell_ids(
#   gp_name = "K170",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   threshold_df = threshold_results_subset_manual,
#   selected_proteins = select_proteins,
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells),
#   missing_threshold_action = "skip"
# )

# seurat_meta_filtered %>%
#   filter(cellID %in% high_cells_k170) %>%
#   group_by(annotation_level1) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))
# seurat_meta_filtered %>%
#   filter(cellID %in% high_cells_k170) %>%
#   group_by(annotation_level2) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))

# p12 <- plot_gated_gp_vs_protein(
#   gp_name = "K171",
#   df_markers = df_markers2,
#   protein_mat = protein_mat_normalized_lognorm,
#   loading_mat = L_pm_filtered,
#   umap_emb = umap_result,
#   missing_threshold_action = "skip", # Skip markers with missing thresholds
#   threshold_df = threshold_results_subset_manual, # use manual thresholds
#   exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells), # Exclude thymocyte, proliferating, and miniverse cells
#   selected_proteins = select_proteins, # Use curated list of good proteins
#   loading_q = NULL
# )
# ggsave(
#   p12,
#   filename = paste0(
#     figure_path,
#     "gated_umap/GP171_protein_loading_gated_umap.pdf"
#   ),
#   width = 12,
#   height = 6
# )

##### Produce a combined figure just show the well-aligned gating
##### from CITEseq_alignment_scores_manual.csv
alignment_scores <- read.csv(
  paste0(data_path, "CITEseq_alignment_scores_manual.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

# produce a subset of plots, just for GPs with alignment_scores$Well.aligned = TRUE
well_aligned_gps <- alignment_scores$GP[alignment_scores$Well.aligned == TRUE]
colnames(L_pm_filtered) <- sub("^K", "GP", colnames(L_pm_filtered))

# exclude GP117, GP163, and GP32, GP166, GP9, GP174, GP101, GP79
# these are either too few/many cells, or having bad visualization on UMAP
well_aligned_gps <- setdiff(well_aligned_gps, "GP117")
well_aligned_gps <- setdiff(well_aligned_gps, c("GP163", "GP32"))
well_aligned_gps <- setdiff(
  well_aligned_gps,
  c("GP166", "GP9", "GP174", "GP101", "GP79")
)

# some GPs to consider enlarging the highlighted points for better visibility
enlarge_gps <- c("GP8", "GP30", "GP170", "GP107")

# pdf(
#   paste0(figure_path, "gated_umap/well_aligned_gps_all.pdf"),
#   width = 12,
#   height = 6
# )
# for (gp in well_aligned_gps) {
#   p <- plot_gated_gp_vs_protein(
#     gp_name = gp,
#     df_markers = df_markers2,
#     protein_mat = protein_mat_normalized_lognorm,
#     loading_mat = L_pm_filtered,
#     umap_emb = umap_result,
#     missing_threshold_action = "skip",
#     threshold_df = threshold_results_subset_manual,
#     exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells),
#     selected_proteins = select_proteins,
#     loading_q = NULL,
#     min_pointsize = if (gp %in% enlarge_gps) 3L else 0L
#   )
#   print(p)
# }
# dev.off()

# ### Produce a full GPs PDF, where GPs listed in order of alignment score
# full_GPs <- alignment_scores$GP
# pdf(
#   paste0(figure_path, "gated_umap/all_gps_sorted.pdf"),
#   width = 12,
#   height = 6
# )
# for (gp in full_GPs) {
#   p <- plot_gated_gp_vs_protein(
#     gp_name = gp,
#     df_markers = df_markers2,
#     protein_mat = protein_mat_normalized_lognorm,
#     loading_mat = L_pm_filtered,
#     umap_emb = umap_result,
#     missing_threshold_action = "skip",
#     threshold_df = threshold_results_subset_manual,
#     exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells),
#     selected_proteins = select_proteins,
#     loading_q = NULL,
#     min_pointsize = if (gp %in% enlarge_gps) 3L else 0L
#   )
#   print(p)
# }
# dev.off()

##### 2-page version: 12 GPs per page, 4 columns x 3 rows
# Pre-generate all plots
# GPs whose highlighted points should be slightly larger even at high cell counts
# take out some GPs that are already presented

# GPs in main Fig 6
GPs_fig6 <- c("GP171", "GP23", "GP12", "GP80")
other_GPs <- setdiff(well_aligned_gps, c(GPs_fig6))
plots_all <- lapply(other_GPs, function(gp) {
  plot_gated_gp_vs_protein(
    gp_name = gp,
    df_markers = df_markers2,
    protein_mat = protein_mat_normalized_lognorm,
    loading_mat = L_pm_filtered,
    umap_emb = umap_result,
    missing_threshold_action = "skip",
    threshold_df = threshold_results_subset_manual,
    exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells),
    selected_proteins = select_proteins,
    loading_q = NULL,
    min_pointsize = if (gp %in% enlarge_gps) 3L else 0L
  )
})

# Each GP unit = 2 side-by-side panels (protein gating | loading); 2 GPs per row
# separated by a thin vertical bar, 6 rows per page = 12 GPs per page.
# 13 x 18 gives ~3 in per UMAP panel; scattermore raster detail is pixel-controlled.
out_2page <- paste0(figure_path, "gated_umap/gallery_2page.pdf")
dir.create(dirname(out_2page), recursive = TRUE, showWarnings = FALSE)
graphics.off()
v_bar <- ggplot() +
  theme_void() +
  theme(plot.background = element_rect(fill = "grey60", color = NA))
cairo_pdf(
  out_2page,
  width = 13,
  height = 18,
  onefile = TRUE
)
showtext::showtext_begin()
n_pages <- ceiling(length(plots_all) / 12)
for (page_idx in seq_len(n_pages)) {
  idx_start <- (page_idx - 1) * 12 + 1
  page_plots <- plots_all[idx_start:min(page_idx * 12, length(plots_all))]
  rows_list <- lapply(1:6, function(r) {
    i_left <- (r - 1) * 2 + 1
    i_right <- (r - 1) * 2 + 2
    gp_left <- if (i_left <= length(page_plots)) {
      page_plots[[i_left]]
    } else {
      plot_spacer()
    }
    gp_right <- if (i_right <= length(page_plots)) {
      page_plots[[i_right]]
    } else {
      plot_spacer()
    }
    (gp_left | v_bar | gp_right) + plot_layout(widths = c(1, 0.03, 1))
  })
  combined <- wrap_plots(rows_list, ncol = 1)
  print(combined)
}
showtext::showtext_end()
dev.off()

##### Sorted gallery: GPs in numerical order, one PDF per page, panels labeled a/b/c...
sort_idx <- order(as.integer(sub("^GP", "", other_GPs)))
plots_sorted <- plots_all[sort_idx]

out_sorted_dir <- paste0(figure_path, "gated_umap/gallery_sorted_pages/")
dir.create(out_sorted_dir, recursive = TRUE, showWarnings = FALSE)

n_pages_sorted <- ceiling(length(plots_sorted) / 12)

for (page_idx in seq_len(n_pages_sorted)) {
  idx_start <- (page_idx - 1) * 12 + 1
  idx_end <- min(page_idx * 12, length(plots_sorted))
  page_plots <- plots_sorted[idx_start:idx_end]
  n_on_page <- length(page_plots)

  # Add panel letter (a, b, c...) to the protein-gating (left) sub-panel of each GP unit
  labeled_page_plots <- lapply(seq_len(n_on_page), function(i) {
    unit <- page_plots[[i]]
    p1_labeled <- unit[[1]] +
      labs(tag = letters[i]) +
      theme(plot.tag = element_text(size = 14, face = "bold"))
    p1_labeled + unit[[2]] + plot_layout(ncol = 2)
  })

  rows_list <- lapply(1:6, function(r) {
    i_left <- (r - 1) * 2 + 1
    i_right <- (r - 1) * 2 + 2
    gp_left <- if (i_left <= n_on_page) labeled_page_plots[[i_left]] else plot_spacer()
    gp_right <- if (i_right <= n_on_page) labeled_page_plots[[i_right]] else plot_spacer()
    (gp_left | v_bar | gp_right) + plot_layout(widths = c(1, 0.03, 1))
  })

  combined <- wrap_plots(rows_list, ncol = 1)
  out_page <- sprintf("%sgallery_sorted_page%02d.pdf", out_sorted_dir, page_idx)
  graphics.off()
  cairo_pdf(out_page, width = 13, height = 18)
  showtext::showtext_begin()
  print(combined)
  showtext::showtext_end()
  dev.off()
}
