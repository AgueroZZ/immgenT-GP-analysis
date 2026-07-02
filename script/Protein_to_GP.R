library(ggplot2)
library(dplyr)
library(tidyr)
data_path <- "data/"
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered),]
protein_mat_normalized_lognorm <- readRDS(paste0(data_path, "protein_mat_normalized_lognorm.rds"))
proteins_quality <- read.csv(paste0(data_path, "TableS4_citeseq_qc_20250513.csv"), header = TRUE, stringsAsFactors = FALSE, skip = 1)
cells_citeseq <- seurat_meta_filtered$cellID[seurat_meta_filtered$cite_seq]

# keep only cells in cells_citeseq
L_pm_filtered <- L_pm_filtered[intersect(rownames(L_pm_filtered), cells_citeseq),]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[intersect(rownames(protein_mat_normalized_lognorm), cells_citeseq),]

# remove proteins with poor quality
poor_proteins <- proteins_quality$protein[proteins_quality$classification == "poor"]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[, !colnames(protein_mat_normalized_lognorm) %in% poor_proteins]







# --- Configuration ---
working_dir <- "figures/Figures_Protein_GP/"
n_bins      <-  4
threshold   <- 0.1
setwd(working_dir)

# Create subdirectories for organized storage
dir.create("Boxplots_PVE", showWarnings = FALSE)
dir.create("Barplots_ChiSq", showWarnings = FALSE)

# --- Analysis Functions ---

# 1. Quantile Binning + ANOVA (PVE)
get_quantile_pve <- function(p_vec, l_vec, p_name, g_name, bins) {
  common <- intersect(names(p_vec), names(l_vec))
  df <- data.frame(
    p = as.numeric(p_vec[common]),
    l = as.numeric(l_vec[common])
  ) %>% filter(!is.na(p) & !is.na(l))

  df$bin <- factor(ntile(df$p, bins), levels = 1:bins)

  fit_aov     <- aov(l ~ bin, data = df)
  sum_sq      <- summary(fit_aov)[[1]][,"Sum Sq"]
  pve_val     <- sum_sq[1] / sum(sum_sq)

  summ <- df %>% group_by(bin) %>% summarise(m = median(l), .groups = 'drop')
  p <- ggplot(df, aes(x = bin, y = l)) +
    geom_boxplot(fill = "steelblue", alpha = 0.3, outlier.shape = NA) +
    geom_line(data = summ, aes(x = as.numeric(bin), y = m), color = "red", size = 1, group = 1) +
    labs(title = paste(p_name, "predicting", g_name),
         subtitle = sprintf("PVE (Eta-squared): %.4f", pve_val),
         x = paste(p_name, "Quantile"), y = "GP Loading") +
    theme_classic()

  return(list(stat = pve_val, plot = p))
}

# 2. Active Proportion + Chi-Square
get_active_chi <- function(p_vec, l_vec, p_name, g_name, bins, thresh) {
  common <- intersect(names(p_vec), names(l_vec))
  df <- data.frame(
    p = as.numeric(p_vec[common]),
    l = as.numeric(l_vec[common])
  ) %>% filter(!is.na(p) & !is.na(l))

  curr_thresh <- if(is.null(thresh)) quantile(df$l, 0.8) else thresh
  df$active   <- ifelse(df$l >= curr_thresh, 1, 0)
  df$bin      <- factor(ntile(df$p, bins), levels = 1:bins)

  tab <- table(df$bin, df$active)
  chi_stat <- if(ncol(tab) < 2) 0 else chisq.test(tab)$statistic[[1]]

  summ <- df %>% group_by(bin) %>%
    summarise(prop = mean(active), n = n(), .groups = 'drop') %>%
    mutate(se = sqrt(prop*(1-prop)/n))

  p <- ggplot(summ, aes(x = bin, y = prop)) +
    geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.6) +
    geom_errorbar(aes(ymin = prop-se, ymax = prop+se), width = 0.2) +
    labs(title = paste(p_name, "Activity Prediction:", g_name),
         subtitle = sprintf("Chi-Sq Stat: %.2f | Threshold: %.2f", chi_stat, curr_thresh),
         x = paste(p_name, "Quantile"), y = "Proportion Active Cells") +
    theme_classic()

  return(list(stat = chi_stat, plot = p))
}

# --- Batch Execution ---

proteins <- colnames(protein_mat_normalized_lognorm)
gps      <- colnames(L_pm_filtered)

pve_matrix <- matrix(NA, nrow = length(proteins), ncol = length(gps), dimnames = list(proteins, gps))
chi_matrix <- matrix(NA, nrow = length(proteins), ncol = length(gps), dimnames = list(proteins, gps))

for (p_name in proteins) {
  message("Processing protein: ", p_name)
  p_vec <- protein_mat_normalized_lognorm[, p_name]

  # Open one PDF per protein for each analysis type
  pdf(file.path("Boxplots_PVE", paste0(p_name, "_vs_all_GPs_boxplot.pdf")), width = 8, height = 6)
  dev_box <- dev.cur()

  pdf(file.path("Barplots_ChiSq", paste0(p_name, "_vs_all_GPs_barplot.pdf")), width = 8, height = 6)
  dev_bar <- dev.cur()

  for (g_name in gps) {
    l_vec <- L_pm_filtered[, g_name]

    # Run PVE
    res_pve <- get_quantile_pve(p_vec, l_vec, p_name, g_name, n_bins)
    pve_matrix[p_name, g_name] <- res_pve$stat
    dev.set(dev_box); print(res_pve$plot)

    # Run Chi-Square
    res_chi <- get_active_chi(p_vec, l_vec, p_name, g_name, n_bins, threshold)
    chi_matrix[p_name, g_name] <- res_chi$stat
    dev.set(dev_bar); print(res_chi$plot)
  }

  dev.off(dev_box)
  dev.off(dev_bar)
}

# Save final summary matrices
saveRDS(pve_matrix, "pve_results_matrix.rds")
saveRDS(chi_matrix, "chisq_results_matrix.rds")

message("Done! Check 'Boxplots_PVE' and 'Barplots_ChiSq' folders.")









library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel) # For labeling top GPs

selected_proteins <- proteins_quality$protein[proteins_quality$classification != "poor"]
exclude_proteins <- c("CD19", "CD34", "CD45.1","CD45.2", "CD138", "TCRVA2", "TER119")
thy11_proteins <- grep("THY1.1", selected_proteins, value = TRUE)
isotype_proteins <- grep("^Isotype", selected_proteins, value = TRUE)
selected_proteins <- setdiff(selected_proteins, c(exclude_proteins, thy11_proteins, isotype_proteins))
chi_matrix <- chi_matrix[selected_proteins,]

plot_manhattan_swarm <- function(input_matrix,
                                 title = "Predictive Power across all GPs",
                                 y_label = "Statistic Value",
                                 label_threshold_quantile = 0.99,
                                 global_threshold = NULL,
                                 manual_colors = c("steelblue", "darkorange", "forestgreen", "purple"),
                                 n_colors = 2) {

  # 1. Prepare Data
  unique_features <- rownames(input_matrix)
  plot_df <- as.data.frame(input_matrix) %>%
    mutate(Feature = rownames(.)) %>%
    pivot_longer(cols = -Feature, names_to = "GP", values_to = "Value")

  # 2. Add Color Groups and Threshold Flags
  feature_info <- data.frame(
    Feature = unique_features,
    ColorGroup = as.factor(rep(1:n_colors, length.out = length(unique_features)))
  )

  # Map the actual colors to the features for axis coloring later
  color_map <- manual_colors[1:n_colors]
  feature_info$ActualColor <- color_map[feature_info$ColorGroup]

  plot_df <- plot_df %>%
    left_join(feature_info, by = "Feature") %>%
    mutate(Feature = factor(Feature, levels = unique_features))

  # 3. Handle Highlighting and Labeling
  # Logic: Label if it's a top quantile hit OR over the global threshold
  top_hits <- plot_df %>%
    group_by(Feature) %>%
    filter(Value == max(Value)) %>% # Only the top GP per protein
    ungroup()

  if (!is.null(label_threshold_quantile)) {
    q_val <- quantile(plot_df$Value, label_threshold_quantile, na.rm = TRUE)
    top_hits <- top_hits %>% filter(Value >= q_val)
  }

  if (!is.null(global_threshold)) {
    # If global threshold is provided, ensure those are also labeled
    global_hits <- plot_df %>% filter(Value >= global_threshold)
    top_hits <- bind_rows(top_hits, global_hits) %>% distinct()

    # Add a column for highlighting points above global threshold
    plot_df$is_high <- plot_df$Value >= global_threshold
  } else {
    plot_df$is_high <- FALSE
  }

  # 4. The Plot
  p <- ggplot(plot_df, aes(x = Feature, y = Value, color = ColorGroup)) +
    # Use different sizes/alphas for points above global threshold
    geom_jitter(aes(alpha = ifelse(is_high, 1, 0.4),
                    size = ifelse(is_high, 1.2, 0.7)),
                width = 0.3) +

    geom_text_repel(data = top_hits, aes(label = GP),
                    size = 3, max.overlaps = 20,
                    segment.color = "grey50", show.legend = FALSE) +

    scale_color_manual(values = manual_colors[1:n_colors]) +
    scale_size_identity() +
    scale_alpha_identity() +

    labs(title = title,
         x = "Protein Marker",
         y = y_label) +

    theme_minimal() +
    theme(
      # Color the X-axis text based on the Feature's assigned color
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7,
                                 color = feature_info$ActualColor),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  # If a global threshold is set, add a horizontal line for it
  if (!is.null(global_threshold)) {
    p <- p + geom_hline(yintercept = global_threshold, linetype = "dotted", color = "red", alpha = 0.5)
  }

  return(p)
}

p1 <- plot_manhattan_swarm(chi_matrix,
                           title = "Chi-Square Statistics: Protein vs GP",
                           y_label = "Chi-Square Stat",
                           # global_threshold = 50000,
                           label_threshold_quantile = 0.99,
                           n_colors = 3)
print(p1)
ggsave("chi_square_manhattan_swarm.pdf",
       p1,
       width = 10, height = 6)


pve_matrix <- pve_matrix[selected_proteins,]
p2 <- plot_manhattan_swarm(pve_matrix,
                           title = "PvE Statistics: Protein vs GP",
                           y_label = "PvE Stat",
                           # global_threshold = 50000,
                           label_threshold_quantile = 0.99,
                           n_colors = 3)
print(p2)
ggsave("pve_manhattan_swarm.pdf",
       p2,
       width = 10, height = 6)











