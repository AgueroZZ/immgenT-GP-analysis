nmf       <- FALSE
backfit   <- FALSE
outfile1   <- "fit1_summary_12percent.qs"
outfile2   <- "fit2_summary_12percent.qs"
data_path <- "/project2/mstephens/immgent"

library(tools)
library(qs)
library(Matrix)
library(fastTopics)
library(flashier)
library(parallel)



## ---------- helpers ----------
sample_k_subsets <- function(cells_pool, total_cells, frac = 1/8, K = 6, seed = 1234) {
  ksize_target <- ceiling(total_cells * frac)
  ksize_max    <- floor(length(cells_pool) / K)
  ksize        <- min(ksize_target, ksize_max)

  if (ksize < 1) stop("Not enough unused cells to form subsets (ksize < 1).")

  if (ksize < ksize_target) {
    message(sprintf(
      "Not enough unused cells for target size. target ksize=%d, actual ksize=%d (pool=%d, K=%d).",
      ksize_target, ksize, length(cells_pool), K
    ))
  }

  set.seed(seed)
  perm <- sample(cells_pool, size = K * ksize, replace = FALSE)

  out <- split(perm, rep(seq_len(K), each = ksize))
  names(out) <- paste0("subset_", seq_len(K))
  out
}
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


fit_remaining_eighths <- function(data_path,
                                  outfile1, outfile2,
                                  outfiles_rest,   # length 6: fit3..fit8
                                  nmf = FALSE, backfit = FALSE,
                                  frac = 1/8,
                                  seed_cells = 1234,
                                  S_seeds = NULL,  # length 6, optional
                                  greedy_Kmax = 200,
                                  preset = "balanced",
                                  verbose = 2,
                                  mc.cores = 2,
                                  skip_existing = TRUE) {

  stopifnot(length(outfiles_rest) == 6)

  # load big matrix once
  Xall <- qread(file.path(data_path, "flashier_snmf_matrix.qs"))
  all_cells <- rownames(Xall)
  total_cells <- length(all_cells)

  # load existing fits (first two eighths)
  fit1 <- qread(outfile1)
  fit2 <- qread(outfile2)

  # unused cells
  unused <- get_unused_cells(all_cells, list(fit1, fit2))

  # build 6 disjoint subsets from unused
  subsets <- sample_k_subsets(unused, total_cells = total_cells, frac = frac, K = 6, seed = seed_cells)

  # default seeds for S if not provided
  if (is.null(S_seeds)) S_seeds <- 3001 + 0:5 * 1000
  stopifnot(length(S_seeds) == 6)

  message(sprintf(
    "Total cells: %d | Used (fit1+fit2): %d | Unused pool: %d | Each subset size: %d",
    total_cells,
    length(union(rownames(fit1$L_pm), rownames(fit2$L_pm))),
    length(unused),
    length(subsets[[1]])
  ))

  # assemble tasks
  tasks <- lapply(seq_len(6), function(i) {
    cells_i <- subsets[[i]]
    list(
      i = i,
      cells = cells_i,
      X = Xall[cells_i, , drop = FALSE],
      outfile = outfiles_rest[[i]],
      S_seed = S_seeds[[i]]
    )
  })

  # optionally skip ones already done
  if (skip_existing) {
    keep <- vapply(tasks, function(t) !file.exists(t$outfile), logical(1))
    if (!all(keep)) {
      message("Skipping existing outputs: ",
              paste(basename(vapply(tasks[!keep], `[[`, "", "outfile")), collapse = ", "))
    }
    tasks <- tasks[keep]
  }

  # nothing to do
  if (length(tasks) == 0) {
    message("All remaining 6 outputs already exist. Nothing to do.")
    return(invisible(NULL))
  }

  # run in parallel (Linux compute node: mclapply OK)
  res <- mclapply(tasks, function(t) {
    fit_flash_and_save(
      t$X, outfile = t$outfile,
      nmf = nmf, backfit = backfit,
      greedy_Kmax = greedy_Kmax,
      S_seed = t$S_seed,
      qsave_preset = preset,
      verbose = verbose
    )
  }, mc.cores = mc.cores, mc.set.seed = TRUE)

  # return summaries + cell lists (for reproducibility)
  invisible(list(
    results = res,
    subsets = lapply(tasks, `[[`, "cells"),
    outfiles = lapply(tasks, `[[`, "outfile")
  ))
}


# shifted_log_counts <- qread(file.path(data_path,"flashier_snmf_matrix.qs"))
#
# # cell-names
# set.seed(1234)
# all_cells <- rownames(shifted_log_counts)
# selected_indices <- sort(sample(1:length(all_cells),
#                                 ceiling(length(all_cells)/8),
#                                 replace = FALSE))
# first_half <- shifted_log_counts[selected_indices,]
# # remove genes that are entirely zero
# genes_removed_1 <- colSums(first_half) == 0
# first_half <- first_half[,!genes_removed_1]
#
# ## Model fitting
# n1 <- nrow(first_half)
# x1 <- rpois(1e7, 1/n1)
# s1 <- sd(log(x1 + 1))
# if (nmf) {
#
#   # Fit an NMF to the shifted-log counts.
#   fl_nmf <- flash(first_half,
#                   ebnm_fn = ebnm_point_exponential,
#                   greedy_Kmax = 20,var_type = 2,S = s1,
#                   backfit = FALSE,verbose = 2)
#   fl_nmf <- flash_backfit(fl_nmf,extrapolate = FALSE,maxiter = 100,verbose = 2)
#   fl_nmf <- flash_backfit(fl_nmf,extrapolate = TRUE,maxiter = 100,verbose = 2)
#   session_info <- sessionInfo()
#   fl_nmf_ldf <- ldf(fl_nmf,type = "i")
#   save(list = c("fl_nmf_ldf","session_info"),
#        file = outfile)
#   resaveRdaFiles(outfile)
# } else {
#
#   # Fit a semi-NMF to the shifted-log counts.
#   fit1 <- flash(first_half,
#                 ebnm_fn = c(ebnm_point_exponential, ebnm_point_laplace),
#                 var_type = 2,
#                 S = s1,
#                 backfit = FALSE,
#                 greedy_Kmax = 200)
#   if (backfit) {
#     fit1 <- flash_backfit(fit1,extrapolate = FALSE,maxiter = 100,verbose = 2)
#     fit1 <- flash_backfit(fit1,extrapolate = TRUE,maxiter = 100,verbose = 2)
#   }
#   fit1_summary <- list(L_pm = fit1$L_pm,
#                        F_pm = fit1$F_pm,
#                        elbo = fit1$elbo,
#                        residuals_sd = fit1$residuals_sd,
#                        pve = fit1$pve)
#   qsave(fit1_summary,file = outfile1,preset = "balanced")
# }
#
# # second half, sample a quarter of cells again, but not overlapping with first half
# set.seed(5678)
# remaining_indices <- setdiff(1:length(all_cells), selected_indices)
# selected_indices_2 <- sort(sample(remaining_indices,
#                                   ceiling(length(all_cells)/8),
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


outfile1 <- "fit1_summary_12percent.qs"
outfile2 <- "fit2_summary_12percent.qs"

outfiles_rest <- sprintf("fit%d_summary_12percent.qs", 3:8)

Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

out <- fit_remaining_eighths(
  data_path = data_path,
  outfile1 = outfile1,
  outfile2 = outfile2,
  outfiles_rest = outfiles_rest,
  nmf = nmf,
  backfit = backfit,
  greedy_Kmax = 200,
  mc.cores = 2,
  seed_cells = 1234
)


