library(flashier)
library(fastTopics)
library(Matrix)

nmf       <- FALSE
backfit   <- TRUE
outfile   <- "fit_K200_longer.qs" 
data_path <- "/project2/mstephens/immgent"

library(tools)
library(qs)
library(Matrix)
library(fastTopics)
library(flashier)
shifted_log_counts <- qread(file.path(data_path,"flashier_snmf_matrix.qs"))

# cell-names
set.seed(1234)
## Model fitting
n <- nrow(shifted_log_counts)
x <- rpois(1e7, 1/n)
s <- sd(log(x + 1))
if (nmf) {
  
  # Fit an NMF to the shifted-log counts.
  fl_nmf <- flash(shifted_log_counts,
                  ebnm_fn = ebnm_point_exponential,
                  greedy_Kmax = 20,var_type = 2,S = s,
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
  fit <- flash(shifted_log_counts,
                ebnm_fn = c(ebnm_point_exponential, ebnm_point_laplace),
                var_type = 2,
                S = s,
                backfit = FALSE,
                greedy_Kmax = 200)
  if (backfit) {
    fit <- flash_backfit(fit, extrapolate = FALSE, maxiter = 250, verbose = 2)
    fit <- flash_backfit(fit, extrapolate = TRUE, maxiter = 250, verbose = 2)
  }
  fit_summary <- list(L_pm = fit$L_pm,
                       F_pm = fit$F_pm,
                       elbo = fit$elbo,
                       residuals_sd = fit$residuals_sd,
                       pve = fit$pve)
  qsave(fit_summary,file = outfile,preset = "balanced")
}