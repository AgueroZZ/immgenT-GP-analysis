# Pipeline step 5: per-IGT (per-batch) GP reproducibility validation.
#
# Two stages:
#   Stage A (not runnable here -- see GAP note): for each IGT batch,
#     independently re-fit the flashier semi-NMF model on just that batch's
#     cells, producing data/igt_specific/fit_<IGT>_summary.qs (one per IGT,
#     ~80 files already present in this repo's data/).
#   Stage B: Hungarian-match each IGT-specific fit's factors against the
#     full-data GPs by cosine similarity of gene scores, producing the score
#     matrix consumed by FigureS1.R (panels S1A/S1B).
#
# GAP (Stage A): ported from code/validate_by_experiments.R, this stage
# needs data/flashier_snmf_matrix.qs, which (like flashier_snmf.rds itself,
# see 01_extract_data.R) is not produced by any script in this repository --
# presumably a lightly-filtered copy of shifted_log_counts.qs created
# interactively. The data/igt_specific/*.qs outputs of this stage already
# exist in this repo, so Stage B can run without re-doing Stage A.
#
# Source: Stage B ported from code/summarize_validation_by_experiment.R
# (near-identical to the cosine-matching code embedded in
# script/Figure_Saturation.R, which FigureS1.R's plotting logic was ported
# from) -- kept separate here as the "produce the score matrix" step,
# distinct from FigureS1.R's "plot from the cached score matrix" step.

# ============================================================
# Stage A (cluster-scale; not run here, see GAP note above)
# ============================================================
# libs <- c("qs", "Matrix", "flashier", "parallel")
# invisible(sapply(libs, library, character.only = TRUE))
# data_path <- "/project2/mstephens/immgent"
# X_path <- file.path(data_path, "flashier_snmf_matrix.qs")  # GAP: not produced anywhere in this repo
# seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))
# out_dir <- "igt_specific"
# ... (per-IGT flash() fit + qsave(), see code/validate_by_experiments.R for the full script)

# ============================================================
# Stage B: Hungarian-match each IGT-specific fit against the full model
# ============================================================
library(qs)
library(dplyr)

data_path <- "data/"
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))

extract_igt_id_from_name <- function(path) {
  nm <- basename(path)
  m <- regexec("IGT([0-9]+)", nm)
  reg <- regmatches(nm, m)[[1]]
  if (length(reg) < 2) return(NA_integer_)
  as.integer(reg[2])
}

extract_Fpm <- function(obj) {
  if (is.matrix(obj)) return(obj)
  if (is.list(obj) && !is.null(obj$F_pm)) return(obj$F_pm)
  if (is.list(obj) && !is.null(obj$flashier_snmf_summary$F_pm)) return(obj$flashier_snmf_summary$F_pm)
  stop("Cannot find F_pm in this object. Please edit extract_Fpm().")
}

compute_cosine_sim_matrix <- function(L1, L2) {
  norms1 <- apply(L1, 2, function(x) sqrt(sum(x^2)))
  norms1[norms1 == 0] <- Inf
  norms2 <- apply(L2, 2, function(x) sqrt(sum(x^2)))
  norms2[norms2 == 0] <- Inf
  L1_normalized <- t(t(L1) / norms1)
  L2_normalized <- t(t(L2) / norms2)
  crossprod(L1_normalized, L2_normalized)
}

hungarian_match_with_dummy <- function(cos_mat, dummy_cos = 0) {
  K_all <- nrow(cos_mat)
  K_igt <- ncol(cos_mat)
  if (K_igt > K_all) stop("K_igt > K_all: handle separately if needed.")

  if (K_igt == K_all) {
    assignment <- RcppHungarian::HungarianSolver(-1 * cos_mat)
    pairs <- assignment$pairs
    matched_col <- pairs[, 2]
    matched_cos <- cos_mat[pairs]
    unmatched <- rep(FALSE, K_all)
    return(list(matched_col = matched_col, matched_cos = matched_cos, unmatched = unmatched))
  }

  S <- cbind(cos_mat, matrix(dummy_cos, nrow = K_all, ncol = K_all - K_igt))
  assignment <- RcppHungarian::HungarianSolver(-1 * S)
  pairs <- assignment$pairs
  matched_col <- pairs[, 2]
  unmatched <- matched_col > K_igt
  matched_cos <- rep(dummy_cos, K_all)
  real_rows <- which(!unmatched)
  matched_cos[real_rows] <- cos_mat[cbind(real_rows, matched_col[real_rows])]
  matched_col_out <- matched_col
  matched_col_out[unmatched] <- NA_integer_
  list(matched_col = matched_col_out, matched_cos = matched_cos, unmatched = unmatched)
}

match_factors_hungarian_allow_unmatched <- function(F_all, F_igt, scale_cols = TRUE, align_genes = TRUE, dummy_cos = 0) {
  F_all <- as.matrix(F_all)
  F_igt <- as.matrix(F_igt)
  if (is.null(rownames(F_all)) || is.null(rownames(F_igt))) stop("All inputs must have rownames (gene IDs).")

  if (align_genes) {
    common_genes <- intersect(rownames(F_all), rownames(F_igt))
    if (length(common_genes) == 0) stop("No common genes across F_all/F_igt.")
    F_all <- F_all[common_genes, , drop = FALSE]
    F_igt <- F_igt[common_genes, , drop = FALSE]
  }
  if (scale_cols) {
    F_all <- scale(F_all, center = FALSE, scale = TRUE)
    F_igt <- scale(F_igt, center = FALSE, scale = TRUE)
  }
  cos_mat <- compute_cosine_sim_matrix(F_all, F_igt)
  hungarian_match_with_dummy(cos_mat, dummy_cos = dummy_cos)
}

build_cosine_score_matrix_qs_hungarian <- function(F_all, data_path, subfolder = "igt_specific",
                                                    file_pattern = "^fit_IGT[0-9]+_summary\\.qs$",
                                                    scale_cols = TRUE, align_genes = TRUE, dummy_cos = 0) {
  stopifnot(is.matrix(F_all))
  if (is.null(rownames(F_all))) stop("F_all must have rownames (gene IDs).")
  if (is.null(colnames(F_all))) colnames(F_all) <- paste0("K", seq_len(ncol(F_all)))

  igt_folder <- file.path(data_path, subfolder)
  files <- list.files(igt_folder, pattern = file_pattern, full.names = TRUE)
  if (length(files) == 0) stop("No IGT files found in: ", igt_folder)

  igt_ids <- vapply(files, extract_igt_id_from_name, integer(1))
  o <- order(igt_ids)
  files <- files[o]
  igt_ids <- igt_ids[o]

  K_all <- ncol(F_all)
  score_mat <- matrix(0, nrow = K_all, ncol = length(files), dimnames = list(colnames(F_all), paste0("IGT", igt_ids)))
  match_mat <- matrix(NA_integer_, nrow = K_all, ncol = length(files), dimnames = list(colnames(F_all), paste0("IGT", igt_ids)))
  unmatched_mat <- matrix(FALSE, nrow = K_all, ncol = length(files), dimnames = list(colnames(F_all), paste0("IGT", igt_ids)))

  for (j in seq_along(files)) {
    obj <- qs::qread(files[j])
    F_igt <- extract_Fpm(obj)
    hm <- match_factors_hungarian_allow_unmatched(F_all = F_all, F_igt = F_igt, scale_cols = scale_cols, align_genes = align_genes, dummy_cos = dummy_cos)
    score_mat[, j] <- hm$matched_cos
    match_mat[, j] <- hm$matched_col
    unmatched_mat[, j] <- hm$unmatched
  }
  list(score_mat = score_mat, match_mat = match_mat, unmatched_mat = unmatched_mat)
}

F_all <- flashier_snmf_summary$F_pm
colnames(F_all) <- paste0("K", seq_len(ncol(F_all)))

res <- build_cosine_score_matrix_qs_hungarian(F_all = F_all, data_path = data_path, subfolder = "igt_specific", dummy_cos = 0)
score_mat <- res$score_mat

write.csv(score_mat, file.path(data_path, "igt_specific_cosine_scores.csv"), row.names = TRUE)

threshold <- 0.5
validated <- score_mat >= threshold
write.csv(validated, file.path(data_path, "igt_specific_validated_matrix.csv"), row.names = TRUE)

message(sprintf("GPs validated by >= 1 IGT: %d", sum(apply(validated, 1, any))))
message(sprintf("GPs validated by >= 2 IGTs: %d", sum(apply(validated, 1, function(x) sum(x) >= 2))))
