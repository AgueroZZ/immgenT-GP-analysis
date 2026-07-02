compute_cosine_sim_matrix <- function(L1, L2){
  norms1 <- apply(L1, 2, function(x){sqrt(sum(x^2))})
  norms1[norms1 == 0] <- Inf

  norms2 <- apply(L2, 2, function(x){sqrt(sum(x^2))})
  norms2[norms2 == 0] <- Inf

  L1_normalized <- t(t(L1)/norms1)
  L2_normalized <- t(t(L2)/norms2)

  #compute matrix of cosine similarities
  cosine_sim_matrix <- crossprod(L1_normalized, L2_normalized)

  return(cosine_sim_matrix)
}

match_factors_max_cosine <- function(
    F_all,
    F_1,
    F_2,
    scale_cols = TRUE,
    align_genes = TRUE,
    return_full = FALSE
) {
  # ---- checks ----
  stopifnot(is.matrix(F_all) || is.data.frame(F_all))
  stopifnot(is.matrix(F_1)   || is.data.frame(F_1))
  stopifnot(is.matrix(F_2)   || is.data.frame(F_2))
  F_all <- as.matrix(F_all)
  F_1   <- as.matrix(F_1)
  F_2   <- as.matrix(F_2)

  if (is.null(rownames(F_all)) || is.null(rownames(F_1)) || is.null(rownames(F_2))) {
    stop("All inputs must have rownames (gene IDs).")
  }

  # ---- align genes ----
  if (align_genes) {
    common_genes <- Reduce(intersect, list(rownames(F_all), rownames(F_1), rownames(F_2)))
    if (length(common_genes) == 0) stop("No common genes across F_all/F_1/F_2.")
    F_all <- F_all[common_genes, , drop = FALSE]
    F_1   <- F_1[common_genes, , drop = FALSE]
    F_2   <- F_2[common_genes, , drop = FALSE]
  }

  # ---- optional scaling (column-wise) ----
  if (scale_cols) {
    F_all <- scale(F_all, center = FALSE, scale = TRUE)
    F_1   <- scale(F_1,   center = FALSE, scale = TRUE)
    F_2   <- scale(F_2,   center = FALSE, scale = TRUE)
  }

  # ---- cosine similarity + Hungarian matching ----
  cos_1 <- compute_cosine_sim_matrix(F_all, F_1)
  cos_2 <- compute_cosine_sim_matrix(F_all, F_2)

  assignment_1 <- RcppHungarian::HungarianSolver(-1 * cos_1)
  assignment_2 <- RcppHungarian::HungarianSolver(-1 * cos_2)

  pairs_1 <- assignment_1$pairs
  pairs_2 <- assignment_2$pairs

  match_results <- data.frame(
    Factor_All     = pairs_1[, 1],
    Best_Match_F1  = pairs_1[, 2],
    Best_Match_F2  = pairs_2[, 2],
    cosine_F1      = cos_1[pairs_1],
    cosine_F2      = cos_2[pairs_2]
  )

  match_results$max_cosine_similarity <- pmax(match_results$cosine_F1, match_results$cosine_F2)

  score_df <- match_results |>
    dplyr::select(Factor_All, max_cosine_similarity) |>
    dplyr::arrange(dplyr::desc(max_cosine_similarity))

  if (!return_full) return(score_df)

  list(
    score_df = score_df,
    match_results = match_results,
    cos_1 = cos_1,
    cos_2 = cos_2,
    common_genes = if (align_genes) rownames(F_all) else NULL,
    assignment_1 = assignment_1,
    assignment_2 = assignment_2
  )
}

match_factors_min_cosine <- function(
    F_all,
    F_1,
    F_2,
    scale_cols = TRUE,
    align_genes = TRUE,
    return_full = FALSE
) {
  # ---- checks ----
  stopifnot(is.matrix(F_all) || is.data.frame(F_all))
  stopifnot(is.matrix(F_1)   || is.data.frame(F_1))
  stopifnot(is.matrix(F_2)   || is.data.frame(F_2))
  F_all <- as.matrix(F_all)
  F_1   <- as.matrix(F_1)
  F_2   <- as.matrix(F_2)

  if (is.null(rownames(F_all)) || is.null(rownames(F_1)) || is.null(rownames(F_2))) {
    stop("All inputs must have rownames (gene IDs).")
  }

  # ---- align genes ----
  if (align_genes) {
    common_genes <- Reduce(intersect, list(rownames(F_all), rownames(F_1), rownames(F_2)))
    if (length(common_genes) == 0) stop("No common genes across F_all/F_1/F_2.")
    F_all <- F_all[common_genes, , drop = FALSE]
    F_1   <- F_1[common_genes, , drop = FALSE]
    F_2   <- F_2[common_genes, , drop = FALSE]
  }

  # ---- optional scaling (column-wise) ----
  if (scale_cols) {
    F_all <- scale(F_all, center = FALSE, scale = TRUE)
    F_1   <- scale(F_1,   center = FALSE, scale = TRUE)
    F_2   <- scale(F_2,   center = FALSE, scale = TRUE)
  }

  # ---- cosine similarity + Hungarian matching ----
  cos_1 <- compute_cosine_sim_matrix(F_all, F_1)
  cos_2 <- compute_cosine_sim_matrix(F_all, F_2)

  assignment_1 <- RcppHungarian::HungarianSolver(-1 * cos_1)
  assignment_2 <- RcppHungarian::HungarianSolver(-1 * cos_2)

  pairs_1 <- assignment_1$pairs
  pairs_2 <- assignment_2$pairs

  match_results <- data.frame(
    Factor_All     = pairs_1[, 1],
    Best_Match_F1  = pairs_1[, 2],
    Best_Match_F2  = pairs_2[, 2],
    cosine_F1      = cos_1[pairs_1],
    cosine_F2      = cos_2[pairs_2]
  )

  match_results$min_cosine_similarity <- pmin(match_results$cosine_F1, match_results$cosine_F2)

  score_df <- match_results |>
    dplyr::select(Factor_All, min_cosine_similarity) |>
    dplyr::arrange(dplyr::desc(min_cosine_similarity))

  if (!return_full) return(score_df)

  list(
    score_df = score_df,
    match_results = match_results,
    cos_1 = cos_1,
    cos_2 = cos_2,
    common_genes = if (align_genes) rownames(F_all) else NULL,
    assignment_1 = assignment_1,
    assignment_2 = assignment_2
  )
}



match_factors_cosine_F1_F2 <- function(
    F_all,
    F_1,
    F_2,
    scale_cols = TRUE,
    align_genes = TRUE,
    use_abs = FALSE
) {
  stopifnot(is.matrix(F_all) || is.data.frame(F_all))
  stopifnot(is.matrix(F_1)   || is.data.frame(F_1))
  stopifnot(is.matrix(F_2)   || is.data.frame(F_2))
  F_all <- as.matrix(F_all)
  F_1   <- as.matrix(F_1)
  F_2   <- as.matrix(F_2)

  if (is.null(rownames(F_all)) || is.null(rownames(F_1)) || is.null(rownames(F_2))) {
    stop("All inputs must have rownames (gene IDs).")
  }

  if (align_genes) {
    common_genes <- Reduce(intersect, list(rownames(F_all), rownames(F_1), rownames(F_2)))
    if (length(common_genes) == 0) stop("No common genes across F_all/F_1/F_2.")
    F_all <- F_all[common_genes, , drop = FALSE]
    F_1   <- F_1[common_genes, , drop = FALSE]
    F_2   <- F_2[common_genes, , drop = FALSE]
  }

  if (scale_cols) {
    F_all <- scale(F_all, center = FALSE, scale = TRUE)
    F_1   <- scale(F_1,   center = FALSE, scale = TRUE)
    F_2   <- scale(F_2,   center = FALSE, scale = TRUE)
  }

  # match F_all -> F_1 and F_all -> F_2
  cos_1 <- compute_cosine_sim_matrix(F_all, F_1)
  cos_2 <- compute_cosine_sim_matrix(F_all, F_2)
  if (use_abs) {
    cos_1 <- abs(cos_1)
    cos_2 <- abs(cos_2)
  }

  assignment_1 <- RcppHungarian::HungarianSolver(-1 * cos_1)
  assignment_2 <- RcppHungarian::HungarianSolver(-1 * cos_2)

  pairs_1 <- assignment_1$pairs
  pairs_2 <- assignment_2$pairs

  match_results <- data.frame(
    Factor_All    = pairs_1[, 1],
    Best_Match_F1 = pairs_1[, 2],
    Best_Match_F2 = pairs_2[, 2]
  )

  # cosine between the matched factors (F1 vs F2), per Factor_All
  F1_matched <- F_1[, match_results$Best_Match_F1, drop = FALSE]
  F2_matched <- F_2[, match_results$Best_Match_F2, drop = FALSE]

  cos_F1_F2_mat <- compute_cosine_sim_matrix(F1_matched, F2_matched)
  if (use_abs) cos_F1_F2_mat <- abs(cos_F1_F2_mat)

  data.frame(
    Factor_All = match_results$Factor_All,
    cosine_F1_F2 = diag(cos_F1_F2_mat)
  ) %>%
    dplyr::arrange(dplyr::desc(cosine_F1_F2))
}
