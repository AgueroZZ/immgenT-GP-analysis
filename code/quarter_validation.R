nmf       <- FALSE
backfit   <- FALSE
outfile1   <- "fit1_summary_25percent.qs"
outfile2   <- "fit2_summary_25percent.qs"
data_path <- "/project2/mstephens/immgent"

# consider the rest two quarters as well
outfile3   <- "fit3_summary_25percent.qs"
outfile4   <- "fit4_summary_25percent.qs"

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")


library(tools)
library(qs)
library(Matrix)
library(fastTopics)
library(flashier)
# shifted_log_counts <- qread(file.path(data_path,"flashier_snmf_matrix.qs"))


## ---------- helpers ----------

drop_all_zero_genes <- function(X) {
  keep <- Matrix::colSums(X != 0) > 0
  X[, keep, drop = FALSE]
}

estimate_S_from_n <- function(n_rows, seed = 1L, n_sim = 1e7) {
  set.seed(seed)
  x <- rpois(n_sim, 1 / n_rows)
  sd(log(x + 1))
}

get_unused_cells <- function(all_cells, fit_summaries) {
  used <- unique(unlist(lapply(fit_summaries, function(s) rownames(s$L_pm))))
  setdiff(all_cells, used)
}

sample_two_quarters <- function(cells_pool, total_cells, frac = 0.25, seed = 1234) {
  k_target <- ceiling(total_cells * frac)
  k_max <- floor(length(cells_pool) / 2)

  k <- min(k_target, k_max)
  stopifnot(k >= 1)

  set.seed(seed)
  cells3 <- sample(cells_pool, size = k, replace = FALSE)
  remaining <- setdiff(cells_pool, cells3)
  cells4 <- sample(remaining, size = k, replace = FALSE)

  list(q3 = sort(cells3), q4 = sort(cells4), k = k, k_target = k_target)
}

fit_flash_and_save <- function(X, outfile, nmf = FALSE, backfit = FALSE,
                               greedy_Kmax = 200, S_seed = 1L,
                               qsave_preset = "balanced",
                               verbose = 2) {
  X <- drop_all_zero_genes(X)

  n <- nrow(X)
  S <- estimate_S_from_n(n_rows = n, seed = S_seed)

  if (nmf) {
    # NMF-like
    fit <- flash(
      X,
      ebnm_fn = ebnm_point_exponential,
      greedy_Kmax = min(greedy_Kmax, 20),
      var_type = 2,
      S = S,
      backfit = FALSE,
      verbose = verbose
    )
    if (backfit) {
      fit <- flash_backfit(fit, extrapolate = FALSE, maxiter = 100, verbose = verbose)
      fit <- flash_backfit(fit, extrapolate = TRUE,  maxiter = 100, verbose = verbose)
    }

    summary <- list(
      L_pm = fit$L_pm,
      F_pm = fit$F_pm,
      elbo = fit$elbo,
      residuals_sd = fit$residuals_sd,
      pve = fit$pve,
      S = S,
      nmf = TRUE
    )
    qsave(summary, file = outfile, preset = qsave_preset)
    return(invisible(summary))
  } else {
    # semi-NMF
    fit <- flash(
      X,
      ebnm_fn = c(ebnm_point_exponential, ebnm_point_laplace),
      var_type = 2,
      S = S,
      backfit = FALSE,
      greedy_Kmax = greedy_Kmax
    )
    if (backfit) {
      fit <- flash_backfit(fit, extrapolate = FALSE, maxiter = 100, verbose = verbose)
      fit <- flash_backfit(fit, extrapolate = TRUE,  maxiter = 100, verbose = verbose)
    }

    summary <- list(
      L_pm = fit$L_pm,
      F_pm = fit$F_pm,
      elbo = fit$elbo,
      residuals_sd = fit$residuals_sd,
      pve = fit$pve,
      S = S,
      nmf = FALSE
    )
    qsave(summary, file = outfile, preset = qsave_preset)
    return(invisible(summary))
  }
}



# # cell-names
# set.seed(1234)
# all_cells <- rownames(shifted_log_counts)
# selected_indices <- sort(sample(1:length(all_cells),
#                                 ceiling(length(all_cells)/4),
#                                 replace = FALSE))
# first_half <- shifted_log_counts[selected_indices,]
# # remove genes that are entirely zero
# genes_removed_1 <- colSums(first_half) == 0
# first_half <- first_half[,!genes_removed_1]
#
# # ## Model fitting
# # n1 <- nrow(first_half)
# # x1 <- rpois(1e7, 1/n1)
# # s1 <- sd(log(x1 + 1))
# # if (nmf) {
# #
# #   # Fit an NMF to the shifted-log counts.
# #   fl_nmf <- flash(first_half,
# #                   ebnm_fn = ebnm_point_exponential,
# #                   greedy_Kmax = 20,var_type = 2,S = s1,
# #                   backfit = FALSE,verbose = 2)
# #   fl_nmf <- flash_backfit(fl_nmf,extrapolate = FALSE,maxiter = 100,verbose = 2)
# #   fl_nmf <- flash_backfit(fl_nmf,extrapolate = TRUE,maxiter = 100,verbose = 2)
# #   session_info <- sessionInfo()
# #   fl_nmf_ldf <- ldf(fl_nmf,type = "i")
# #   save(list = c("fl_nmf_ldf","session_info"),
# #        file = outfile)
# #   resaveRdaFiles(outfile)
# # } else {
# #
# #   # Fit a semi-NMF to the shifted-log counts.
# #   fit1 <- flash(first_half,
# #                 ebnm_fn = c(ebnm_point_exponential, ebnm_point_laplace),
# #                 var_type = 2,
# #                 S = s1,
# #                 backfit = FALSE,
# #                 greedy_Kmax = 200)
# #   if (backfit) {
# #     fit1 <- flash_backfit(fit1,extrapolate = FALSE,maxiter = 100,verbose = 2)
# #     fit1 <- flash_backfit(fit1,extrapolate = TRUE,maxiter = 100,verbose = 2)
# #   }
# #   fit1_summary <- list(L_pm = fit1$L_pm,
# #                        F_pm = fit1$F_pm,
# #                        elbo = fit1$elbo,
# #                        residuals_sd = fit1$residuals_sd,
# #                        pve = fit1$pve)
# #   qsave(fit1_summary,file = outfile1,preset = "balanced")
# # }
#
# # second half, sample a quarter of cells again, but not overlapping with first half
# set.seed(5678)
# remaining_indices <- setdiff(1:length(all_cells), selected_indices)
# selected_indices_2 <- sort(sample(remaining_indices,
#                                   ceiling(length(all_cells)/4),
#                                   replace = FALSE))
# second_half <- shifted_log_counts[selected_indices_2,]
# # remove genes that are entirely zero
# genes_removed_2 <- colSums(second_half) == 0
# second_half <- second_half[,!genes_removed_2]
#
# ## Model fitting
# n2 <- nrow(second_half)
# x2 <- rpois(1e7, 1/n2)
# s2 <- sd(log(x2 + 1))
# if (nmf) {
#
#   # Fit an NMF to the shifted-log counts.
#   fl_nmf2 <- flash(second_half,
#                   ebnm_fn = ebnm_point_exponential,
#                   greedy_Kmax = 20,var_type = 2,S = s2,
#                   backfit = FALSE,verbose = 2)
#   fl_nmf2 <- flash_backfit(fl_nmf2,extrapolate = FALSE,maxiter = 100,verbose = 2)
#   fl_nmf2 <- flash_backfit(fl_nmf2,extrapolate = TRUE,maxiter = 100,verbose = 2)
#   session_info <- sessionInfo()
#   fl_nmf2_ldf <- ldf(fl_nmf2,type = "i")
#   save(list = c("fl_nmf2_ldf","session_info"),
#        file = outfile2)
#   resaveRdaFiles(outfile2)
# } else {
#   # Fit a semi-NMF to the shifted-log counts.
#   fit2 <- flash(second_half,
#                 ebnm_fn = c(ebnm_point_exponential, ebnm_point_laplace),
#                 var_type = 2,
#                 S = s2,
#                 backfit = FALSE,
#                 greedy_Kmax = 200)
#   if (backfit) {
#     fit2 <- flash_backfit(fit2,extrapolate = FALSE,maxiter = 100,verbose = 2)
#     fit2 <- flash_backfit(fit2,extrapolate = TRUE,maxiter = 100,verbose = 2)
#   }
#   fit2_summary <- list(L_pm = fit2$L_pm,
#                        F_pm = fit2$F_pm,
#                        elbo = fit2$elbo,
#                        residuals_sd = fit2$residuals_sd,
#                        pve = fit2$pve)
#   qsave(fit2_summary,file = outfile2,preset = "balanced")
# }

library(parallel)

fit_two_quarters_mclapply <- function(X3, X4,
                                      outfile3, outfile4,
                                      nmf, backfit, greedy_Kmax,
                                      preset, verbose,
                                      S_seed_q3 = 3001, S_seed_q4 = 4001,
                                      mc.cores = 2) {

  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

  tasks <- list(
    list(X = X3, outfile = outfile3, S_seed = S_seed_q3),
    list(X = X4, outfile = outfile4, S_seed = S_seed_q4)
  )

  res <- mclapply(tasks, function(t) {
    Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

    fit_flash_and_save(
      t$X, outfile = t$outfile,
      nmf = nmf, backfit = backfit,
      greedy_Kmax = greedy_Kmax,
      S_seed = t$S_seed,
      qsave_preset = preset,
      verbose = verbose
    )
  }, mc.cores = mc.cores, mc.set.seed = TRUE)

  list(fit3_summary = res[[1]], fit4_summary = res[[2]])
}


# ---- build X3 / X4 from existing fit1 + fit2 ----
Xall <- qread(file.path(data_path, "flashier_snmf_matrix.qs"))
all_cells <- rownames(Xall)
total_cells <- length(all_cells)

fit1 <- qread(outfile1)
fit2 <- qread(outfile2)

unused   <- get_unused_cells(all_cells, list(fit1, fit2))
quarters <- sample_two_quarters(unused, total_cells = total_cells, frac = 0.25, seed = 1234)

X3 <- Xall[quarters$q3, , drop = FALSE]
X4 <- Xall[quarters$q4, , drop = FALSE]

message(sprintf("Q3 cells: %d; Q4 cells: %d; unused pool: %d",
                length(quarters$q3), length(quarters$q4), length(unused)))


out <- fit_two_quarters_mclapply(
  X3, X4,
  outfile3, outfile4,
  nmf = nmf, backfit = backfit,
  greedy_Kmax = 200,
  preset = "balanced",
  verbose = 2,
  mc.cores = 2
)




