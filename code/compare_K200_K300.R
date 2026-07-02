### Compare the rank K = 200 and K = 300 factorization
library(ggplot2)
library(dplyr)
library(fastTopics)
library(qs)
library(cowplot)
library(ggrepel)
data_path <- "./data/"
code_path <- "./code/"
source(paste0(code_path, "compute_cosine_sim.R"))
source(paste0(code_path, "filtering_membership.R"))
seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))
rank_k200_result_back <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
rank_k200_result <- qs::qread(paste0(data_path, "fit_replicate_summary_noback.qs"))
rank_k300_result <- qs::qread(paste0(data_path, "fit_K300_noback.qs"))

# check if the result without backfitting is close to the result with backfitting
F_K200_back <- rank_k200_result_back$F_pm
F_K200 <- rank_k200_result$F_pm
similarity <- diag(compute_cosine_sim_matrix(F_K200, F_K200_back))
plot(similarity, xlab = "Factor", ylab = "Cosine Similarity", main = "Similarity between K=200 with and without backfitting")
lines(lowess(similarity), col = "red")


# compare the K=200 and K=300 results
F_K300 <- rank_k300_result$F_pm
similarity_K200_K300 <- compute_cosine_sim_matrix(F_K200_back, F_K300)
# For each factor in K=200, find the most similar factor in K=300
max_similarity <- apply(similarity_K200_K300, 1, max)
plot(max_similarity, xlab = "Factor (K=200)",
     ylab = "Max Cosine Similarity with K=300",
     ylim = c(0,1),
     main = "Similarity between K=200 (backfitting) and K=300 (greedy)")
lines(lowess(max_similarity), col = "red")


# apply hungarian algorithm to find the best matching between K=200 and K=300 factors
matched_similiarity <- match_factors_min_cosine(F_K200_back, F_K300, F_K300)
plot(matched_similiarity$min_cosine_similarity, xlab = "Factor (K=200)",
     ylab = "Cosine Similarity with Matched Factor in K=300",
     ylim = c(0,1),
     main = "Similarity between K=200 (backfitting) and K=300 (greedy) after matching")

# take a look at the PVE:
pve_K200_back <- rank_k200_result_back$pve
plot(pve_K200_back, log = "y")
lines(lowess(pve_K200_back), col = "red")

pve_K300 <- rank_k300_result$pve
plot(pve_K300, log = "y")
lines(lowess(pve_K300), col = "red")






