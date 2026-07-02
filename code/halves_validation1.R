# sinteractive -p mstephens --account=pi-mstephens --time=96:00:00 \
#   --mem=120G -c 8 (for nmf = FALSE)
# sinteractive -p mstephens --account=pi-mstephens --time=24:00:00 \
#   --mem=36G -c 8 (for nmf = TRUE)
# module load R/4.2.0
# export OMP_NUM_THREADS=8
# > .libPaths()[1]
# [1] "/home/pcarbo/R_libs_4_20"
nmf       <- FALSE
backfit   <- FALSE
outfile   <- "fit1_nmf_k=20.RData" # "fit1_summary.qs"
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
first_half <- shifted_log_counts[selected_indices,]

# remove genes that are entirely zero
genes_removed_1 <- colSums(first_half) == 0
first_half <- first_half[,!genes_removed_1]

## Model fitting
n1 <- nrow(first_half)
x1 <- rpois(1e7, 1/n1)
s1 <- sd(log(x1 + 1))
if (nmf) {

  # Fit an NMF to the shifted-log counts.
  fl_nmf <- flash(first_half,
                  ebnm_fn = ebnm_point_exponential,
                  greedy_Kmax = 20,var_type = 2,S = s1,
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
  fit1 <- flash(first_half,
                ebnm_fn = c(ebnm_point_exponential, ebnm_point_laplace),
                var_type = 2,
                S = s1,
                backfit = FALSE,
                greedy_Kmax = 200)
  if (backfit) {
    fit1 <- flash_backfit(fit1,extrapolate = FALSE,maxiter = 100,verbose = 2)
    fit1 <- flash_backfit(fit1,extrapolate = TRUE,maxiter = 100,verbose = 2)
  }
  fit1_summary <- list(L_pm = fit1$L_pm,
                       F_pm = fit1$F_pm,
                       elbo = fit1$elbo,
                       residuals_sd = fit1$residuals_sd,
                       pve = fit1$pve)
  qsave(fit1_summary,file = outfile,preset = "balanced")
}
