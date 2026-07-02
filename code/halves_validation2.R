# sinteractive -p mstephens --account=pi-mstephens --time=96:00:00 \
#   --mem=120G -c 8 (for nmf = FALSE)
# sinteractive -p mstephens --account=pi-mstephens --time=24:00:00 \
#   --mem=36G -c 8 (for nmf = TRUE)
# module load R/4.2.0
# export OMP_NUM_THREADS=8
# > .libPaths()[1]
# [1] "/home/pcarbo/R_libs_4_20"
nmf       <- TRUE # FALSE
backfit   <- TRUE # TRUE
outfile   <- "fit2_nmf_k=20.RData" # "fit2_summary.qs"
data_path <- "/project2/mstephens/immgent"

library(tools)
library(qs)
library(Matrix)
library(fastTopics)
library(flashier)
shifted_log_counts <- qread(file.path(data_path,"flashier_snmf_matrix.qs"))

# cell-names
set.seed(1234)
all_cells <- rownames(shifted_log_counts)
selected_indices <- sort(sample(1:length(all_cells),
                                ceiling(length(all_cells)/2),
                                replace = FALSE))
second_half <- shifted_log_counts[-selected_indices,]

# remove genes that are entirely zero
genes_removed_2 <- colSums(second_half) == 0
second_half <- second_half[,!genes_removed_2]

## Model fitting
n2 <- nrow(second_half)
x2 <- rpois(1e7, 1/n2)
s2 <- sd(log(x2 + 1))
if (nmf) {

  # Fit an NMF to the shifted-log counts.
  fl_nmf <- flash(second_half,
                  ebnm_fn = ebnm_point_exponential,
                  greedy_Kmax = 20,var_type = 2,S = s2,
                  backfit = FALSE,verbose = 2)
  fl_nmf <- flash_backfit(fl_nmf,extrapolate = FALSE,maxiter = 100,verbose = 2)
  fl_nmf <- flash_backfit(fl_nmf,extrapolate = TRUE,maxiter = 100,verbose = 2)
  session_info <- sessionInfo()
  fl_nmf_ldf <- ldf(fl_nmf,type = "i")
  save(list = c("fl_nmf_ldf","session_info"),
       file = outfile)
  resaveRdaFiles(outfile)
} else {

  # Fit a semi-NMF to the shifted-log counts.
  fit2 <- flash(second_half,
                ebnm_fn = c(ebnm_point_exponential, ebnm_point_laplace),
                var_type = 2,
                S = s2,
                backfit = FALSE,
                greedy_Kmax = 200)
  if (backfit) {
    fit2 <- flash_backfit(fit2,extrapolate = FALSE,maxiter = 100,verbose = 2)
    fit2 <- flash_backfit(fit2,extrapolate = TRUE,maxiter = 40,verbose = 2)
  }
  fit2_summary <- list(L_pm = fit2$L_pm,
                       F_pm = fit2$F_pm,
                       elbo = fit2$elbo,
                       residuals_sd = fit2$residuals_sd,
                       pve = fit2$pve)
  qsave(fit2_summary,file = outfile,preset = "balanced")
}
