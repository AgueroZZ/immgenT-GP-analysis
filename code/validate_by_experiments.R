## find all the experiment
library(ggplot2)
library(dplyr)
library(fastTopics)
library(qs)
library(cowplot)
library(ggrepel)
data_path <- "/project2/mstephens/immgent"

## igt_factorization.R

Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

library(qs)
library(Matrix)
library(flashier)
library(parallel)

data_path <- "/project2/mstephens/immgent"
X_path    <- file.path(data_path, "flashier_snmf_matrix.qs")
seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))

out_dir <- "igt_specific"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- helpers ---
drop_all_zero_genes <- function(X) {
  keep <- Matrix::colSums(X != 0) > 0
  X[, keep, drop = FALSE]
}

estimate_S_from_n <- function(n_rows, seed = 1L, n_sim = 1e7) {
  set.seed(seed)
  x <- rpois(n_sim, 1 / n_rows)
  sd(log(x + 1))
}

fit_flash_summary <- function(X, nmf = FALSE, backfit = FALSE,
                              greedy_Kmax = 200, S_seed = 1L, verbose = 2) {
  X <- drop_all_zero_genes(X)
  n <- nrow(X)
  S <- estimate_S_from_n(n_rows = n, seed = S_seed)

  if (nmf) {
    fit <- flash(
      X,
      ebnm_fn = ebnm_point_exponential,
      greedy_Kmax = min(greedy_Kmax, 20),
      var_type = 2,
      S = S,
      backfit = FALSE,
      verbose = verbose
    )
  } else {
    fit <- flash(
      X,
      ebnm_fn = c(ebnm_point_exponential, ebnm_point_laplace),
      var_type = 2,
      S = S,
      backfit = FALSE,
      greedy_Kmax = greedy_Kmax
    )
  }

  if (backfit) {
    fit <- flash_backfit(fit, extrapolate = FALSE, maxiter = 100, verbose = verbose)
    fit <- flash_backfit(fit, extrapolate = TRUE,  maxiter = 100, verbose = verbose)
  }

  list(
    L_pm = fit$L_pm,
    F_pm = fit$F_pm,
    elbo = fit$elbo,
    residuals_sd = fit$residuals_sd,
    pve = fit$pve,
    S = S,
    nmf = nmf,
    cells = rownames(fit$L_pm)
  )
}

Xall <- qread(X_path)

# align cells between seurat_meta and Xall
common_cells <- intersect(rownames(seurat_meta), rownames(Xall))
if (length(common_cells) == 0) stop("No overlapping cell names between seurat_meta and Xall.")
seurat_meta <- seurat_meta[common_cells, , drop = FALSE]
Xall <- Xall[common_cells, , drop = FALSE]

igt_levels <- sort(unique(as.character(seurat_meta$IGT)))
tasks <- lapply(igt_levels, function(igt) {
  cells <- rownames(seurat_meta)[as.character(seurat_meta$IGT) == igt]
  outfile <- file.path(out_dir, sprintf("fit_%s_summary.qs", igt))
  list(igt = igt, cells = cells, outfile = outfile)
})

nmf <- FALSE
backfit <- FALSE
greedy_Kmax <- 200
verbose <- 2
min_cells <- 0            # let's not filter by cell number for now
mc.cores <- 4
skip_existing <- TRUE

# skip if it already exists or too few cells
tasks <- Filter(function(t) {
  if (length(t$cells) < min_cells) return(FALSE)
  if (skip_existing && file.exists(t$outfile)) return(FALSE)
  TRUE
}, tasks)

message(sprintf("Total IGTs: %d; to run now: %d", length(igt_levels), length(tasks)))

res <- mclapply(tasks, function(t) {
  Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

  # subset matrix
  X <- Xall[t$cells, , drop = FALSE]

  # use IGT number to set seed (assuming IGT names like "IGT1", "IGT2", ...)
  igt_num <- suppressWarnings(as.integer(gsub("\\D+", "", t$igt)))
  if (is.na(igt_num)) igt_num <- 1L
  S_seed <- 1000L + igt_num

  message(sprintf("[%s] n_cells=%d -> %s", t$igt, nrow(X), t$outfile))

  summ <- fit_flash_summary(
    X,
    nmf = nmf, backfit = backfit,
    greedy_Kmax = greedy_Kmax,
    S_seed = S_seed,
    verbose = verbose
  )

  qsave(summ, file = t$outfile, preset = "balanced")
  return(list(igt = t$igt, outfile = t$outfile, n_cells = length(t$cells)))
}, mc.cores = mc.cores, mc.set.seed = TRUE)

# summary
done <- do.call(rbind, lapply(res, as.data.frame))
print(done)



