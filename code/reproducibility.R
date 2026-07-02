library(ggplot2)
library(dplyr)
library(fastTopics)
library(qs)
library(cowplot)
library(ggrepel)
data_path <- "./data/"
code_path <- "./code/"
source(paste0(code_path, "compute_cosine_sim.R"))
seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
F_all <- flashier_snmf_summary$F_pm

flashier_snmf_summary1 <- qs::qread(paste0(data_path, "fit1_summary_noback.qs"))
flashier_snmf_summary2 <- qs::qread(paste0(data_path, "fit2_summary_noback.qs"))
F_1 <- flashier_snmf_summary1$F_pm; F_2 <- flashier_snmf_summary2$F_pm
score_50_percent_max <- match_factors_max_cosine(F_all, F_1, F_2)
score_50_percent_min <- match_factors_min_cosine(F_all, F_1, F_2)
score_50_percent_between <- match_factors_cosine_F1_F2(F_all, F_1, F_2)

flashier_snmf_summary1 <- qs::qread(paste0(data_path, "fit1_summary_25percent.qs"))
flashier_snmf_summary2 <- qs::qread(paste0(data_path, "fit2_summary_25percent.qs"))
F_1 <- flashier_snmf_summary1$F_pm; F_2 <- flashier_snmf_summary2$F_pm
score_25_percent_max <- match_factors_max_cosine(F_all, F_1, F_2)
score_25_percent_min <- match_factors_min_cosine(F_all, F_1, F_2)
score_25_percent_between <- match_factors_cosine_F1_F2(F_all, F_1, F_2)

flashier_snmf_summary1 <- qs::qread(paste0(data_path, "fit1_summary_12percent.qs"))
flashier_snmf_summary2 <- qs::qread(paste0(data_path, "fit2_summary_12percent.qs"))
F_1 <- flashier_snmf_summary1$F_pm; F_2 <- flashier_snmf_summary2$F_pm
score_12_percent_max <- match_factors_max_cosine(F_all, F_1, F_2)
score_12_percent_min <- match_factors_min_cosine(F_all, F_1, F_2)
score_12_percent_between <- match_factors_cosine_F1_F2(F_all, F_1, F_2)

library(dplyr)
library(tidyr)
library(ggplot2)

threshold_grid <- seq(0.2, 0.8, by = 0.05)

scores_long <- bind_rows(
  score_12_percent_max %>%
    transmute(Factor_All, proportion = "12.5%", score = max_cosine_similarity),
  score_25_percent_max %>%
    transmute(Factor_All, proportion = "25%", score = max_cosine_similarity),
  score_50_percent_max %>%
    transmute(Factor_All, proportion = "50%", score = max_cosine_similarity)
  # score_25_percent_min %>%
  #   transmute(Factor_All, proportion = "25%", score = min_cosine_similarity),
  # score_50_percent_min %>%
  #   transmute(Factor_All, proportion = "50%", score = min_cosine_similarity)
  # score_12_percent_between %>%
  #   transmute(Factor_All, proportion = "12.5%", score = cosine_F1_F2),
  # score_25_percent_between %>%
  #   transmute(Factor_All, proportion = "25%", score = cosine_F1_F2),
  # score_50_percent_between %>%
  #   transmute(Factor_All, proportion = "50%", score = cosine_F1_F2)
)

count_df <- tidyr::crossing(
  threshold = threshold_grid,
  proportion = unique(scores_long$proportion)
) %>%
  left_join(scores_long, by = "proportion") %>%
  group_by(threshold, proportion) %>%
  summarize(
    n_pass = sum(score >= threshold, na.rm = TRUE),
    .groups = "drop"
  )

p <- ggplot(count_df, aes(x = threshold, y = n_pass, color = proportion)) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.8) +
  coord_cartesian(ylim = c(0,200)) +
  # scale_x_reverse(breaks = threshold_grid) +
  labs(
    x = "Reproducibility threshold (max cosine similarity)",
    # x = "Reproducibility threshold (min cosine similarity)",
    # x = "Reproducibility threshold (cosine similarity between runs)",
    y = "Number of factors passing threshold",
    color = "Subset proportion",
    title = "Reproducibility curves across thresholds"
  ) +
  theme_bw()

print(p)

