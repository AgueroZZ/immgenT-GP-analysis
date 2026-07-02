library(qs)
library(Matrix)
library(flashier)
library(parallel)

# data_path <- "/project2/mstephens/immgent"
data_path <- "data/"
code_path <- "code/"
figure_path <- "figures/Figure_Saturation/"
# seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))
out_dir <- "igt_specific"
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))

## Build a matrix of max-cosine scores: rows = GP (K1..K200), cols = IGT (IGT1..)
## Files: data_path/igt_specific/fit_IGT<number>_summary.qs
## Requires: qs

# ---- helper: parse IGT id from filename like "fit_IGT12_summary.qs" ----
extract_igt_id_from_name <- function(path) {
  nm <- basename(path)
  m <- regexec("IGT([0-9]+)", nm)
  reg <- regmatches(nm, m)[[1]]
  if (length(reg) < 2) {
    return(NA_integer_)
  }
  as.integer(reg[2])
}

# ---- helper: extract F_pm matrix from loaded object ----
extract_Fpm <- function(obj) {
  if (is.matrix(obj)) {
    return(obj)
  }
  if (is.list(obj) && !is.null(obj$F_pm)) {
    return(obj$F_pm)
  }
  if (is.list(obj) && !is.null(obj$flashier_snmf_summary$F_pm)) {
    return(obj$flashier_snmf_summary$F_pm)
  }
  stop("Cannot find F_pm in this object. Please edit extract_Fpm().")
}

# ---- helper: compute max cosine per column of A vs all columns of B ----
# use Hungarian algorithm to find best matching, then take max cosine for each column of A
# if the IGT F has fewer columns than F_all, we will still do Hungarian matching, but allow some columns of F_all to be unmatched (score = 0)

library(qs)
library(dplyr)

# ---- helper: parse IGT id from filename like "fit_IGT12_summary.qs" ----
extract_igt_id_from_name <- function(path) {
  nm <- basename(path)
  m <- regexec("IGT([0-9]+)", nm)
  reg <- regmatches(nm, m)[[1]]
  if (length(reg) < 2) {
    return(NA_integer_)
  }
  as.integer(reg[2])
}

# ---- helper: extract F_pm matrix from loaded object ----
extract_Fpm <- function(obj) {
  if (is.matrix(obj)) {
    return(obj)
  }
  if (is.list(obj) && !is.null(obj$F_pm)) {
    return(obj$F_pm)
  }
  if (is.list(obj) && !is.null(obj$flashier_snmf_summary$F_pm)) {
    return(obj$flashier_snmf_summary$F_pm)
  }
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

  if (K_igt > K_all) {
    stop("K_igt > K_all: handle separately if needed.")
  }

  if (K_igt == K_all) {
    assignment <- RcppHungarian::HungarianSolver(-1 * cos_mat)
    pairs <- assignment$pairs
    matched_col <- pairs[, 2]
    matched_cos <- cos_mat[pairs]
    unmatched <- rep(FALSE, K_all)
    return(list(
      matched_col = matched_col,
      matched_cos = matched_cos,
      unmatched = unmatched
    ))
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

  list(
    matched_col = matched_col_out,
    matched_cos = matched_cos,
    unmatched = unmatched
  )
}

match_factors_hungarian_allow_unmatched <- function(
  F_all,
  F_igt,
  scale_cols = TRUE,
  align_genes = TRUE,
  dummy_cos = 0
) {
  F_all <- as.matrix(F_all)
  F_igt <- as.matrix(F_igt)

  if (is.null(rownames(F_all)) || is.null(rownames(F_igt))) {
    stop("All inputs must have rownames (gene IDs).")
  }

  if (align_genes) {
    common_genes <- intersect(rownames(F_all), rownames(F_igt))
    if (length(common_genes) == 0) {
      stop("No common genes across F_all/F_igt.")
    }
    F_all <- F_all[common_genes, , drop = FALSE]
    F_igt <- F_igt[common_genes, , drop = FALSE]
  }

  if (scale_cols) {
    F_all <- scale(F_all, center = FALSE, scale = TRUE)
    F_igt <- scale(F_igt, center = FALSE, scale = TRUE)
  }

  cos_mat <- compute_cosine_sim_matrix(F_all, F_igt) # K_all x K_igt
  hungarian_match_with_dummy(cos_mat, dummy_cos = dummy_cos)
}

# ---- build score matrix across all IGTs (Hungarian + dummy) ----
build_cosine_score_matrix_qs_hungarian <- function(
  F_all,
  data_path,
  subfolder = "igt_specific",
  file_pattern = "^fit_IGT[0-9]+_summary\\.qs$",
  scale_cols = TRUE,
  align_genes = TRUE,
  dummy_cos = 0
) {
  stopifnot(is.matrix(F_all))
  if (is.null(rownames(F_all))) {
    stop("F_all must have rownames (gene IDs).")
  }
  if (is.null(colnames(F_all))) {
    colnames(F_all) <- paste0("K", seq_len(ncol(F_all)))
  }

  igt_folder <- file.path(data_path, subfolder)
  files <- list.files(igt_folder, pattern = file_pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No IGT files found in: ", igt_folder)
  }

  igt_ids <- vapply(files, extract_igt_id_from_name, integer(1))
  o <- order(igt_ids)
  files <- files[o]
  igt_ids <- igt_ids[o]

  K_all <- ncol(F_all)
  score_mat <- matrix(
    0,
    nrow = K_all,
    ncol = length(files),
    dimnames = list(colnames(F_all), paste0("IGT", igt_ids))
  )
  match_mat <- matrix(
    NA_integer_,
    nrow = K_all,
    ncol = length(files),
    dimnames = list(colnames(F_all), paste0("IGT", igt_ids))
  )
  unmatched_mat <- matrix(
    FALSE,
    nrow = K_all,
    ncol = length(files),
    dimnames = list(colnames(F_all), paste0("IGT", igt_ids))
  )

  for (j in seq_along(files)) {
    obj <- qs::qread(files[j])
    F_igt <- extract_Fpm(obj)

    hm <- match_factors_hungarian_allow_unmatched(
      F_all = F_all,
      F_igt = F_igt,
      scale_cols = scale_cols,
      align_genes = align_genes,
      dummy_cos = dummy_cos
    )

    score_mat[, j] <- hm$matched_cos
    match_mat[, j] <- hm$matched_col
    unmatched_mat[, j] <- hm$unmatched
  }

  list(
    score_mat = score_mat,
    match_mat = match_mat,
    unmatched_mat = unmatched_mat
  )
}

# ---- run ----
F_all <- flashier_snmf_summary$F_pm
colnames(F_all) <- paste0("K", seq_len(ncol(F_all)))

res <- build_cosine_score_matrix_qs_hungarian(
  F_all = F_all,
  data_path = data_path,
  subfolder = "igt_specific",
  dummy_cos = 0
)

score_mat <- res$score_mat
match_mat <- res$match_mat
unmatched_mat <- res$unmatched_mat

dim(score_mat)


# save score_mat into a csv file
write.csv(
  score_mat,
  file.path(data_path, "igt_specific_cosine_scores.csv"),
  row.names = TRUE
)


# use a threshold, to see what GPs have been validated in each IGT
threshold <- 0.5
validated <- score_mat >= threshold

write.csv(
  validated,
  file.path(data_path, "igt_specific_validated_matrix.csv"),
  row.names = TRUE
)

# how many GPs are validated by at least one IGT?
sum(apply(validated, 1, any))

# how many GPs are validated by at least two IGT?
sum(apply(validated, 1, function(x) sum(x) >= 2))

library(dplyr)
library(tidyr)
library(ggplot2)

thresholds <- seq(0.2, 0.8, by = 0.1)
X_grid <- 1:50

plot_df <- tidyr::crossing(threshold = thresholds, X = X_grid) %>%
  mutate(
    n_GP = purrr::map2_int(threshold, X, \(t, x) {
      rowSums(score_mat >= t, na.rm = TRUE) |> (\(v) sum(v >= x))()
    })
  )

ggplot(plot_df, aes(x = X, y = n_GP, group = factor(threshold))) +
  geom_line() +
  geom_point(size = 1) +
  # log scale for y-axis if needed
  scale_y_log10() +
  scale_x_log10() +
  labs(
    x = "X (validated by at least X IGTs)",
    y = "Number of GPs",
    color = "Threshold"
  ) +
  aes(color = factor(threshold)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# save as pdf
ggsave(paste0(figure_path, "igt_specific_validated_by_X_IGTs.pdf"), width = 6, height = 4)


thresholds <- seq(0.3, 0.7, by = 0.2)
X_grid <- 0:50
# non-cumulative plot
plot_df_non_cum <- tidyr::crossing(threshold = thresholds, X = X_grid) %>%
  mutate(
    n_GP = purrr::map2_int(threshold, X, \(t, x) {
      counts <- rowSums(score_mat >= t, na.rm = TRUE) # per GP: #IGTs passing threshold t
      sum(counts == x, na.rm = TRUE) # non-cumulative: exactly x
    })
  )

ggplot(
  plot_df_non_cum,
  aes(x = X, y = n_GP, color = factor(threshold), group = factor(threshold))
) +
  geom_line() +
  geom_point(size = 1) +
  # log scale for y-axis if needed
  scale_y_log10() +
  # scale_x_log10() +
  labs(
    x = "X (validated by exactly X IGTs)",
    y = "Number of GPs",
    color = "Threshold"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

ggsave(
  paste0(figure_path, "igt_specific_validated_by_exactly_X_IGTs.pdf"),
  width = 6,
  height = 4
)


# ---- make sure columns are ordered by IGT index ----
igt_idx <- as.integer(gsub("^IGT", "", colnames(score_mat)))
o <- order(igt_idx)
score_mat_ord <- score_mat[, o, drop = FALSE]
igt_idx_ord <- igt_idx[o]

# ---- thresholds you want to compare ----
thresholds <- seq(0.2, 0.8, by = 0.1)

# ---- cumulative validated count function ----
cum_validated_counts <- function(score_mat_ord, threshold) {
  validated <- score_mat_ord >= threshold
  # per GP (row), cumulative "ever validated so far" as we include more IGTs
  ever_validated <- t(apply(validated, 1, cummax)) # 200 x nIGT logical
  colSums(ever_validated) # length nIGT
}

# ---- build plotting data ----
plot_df <- lapply(thresholds, function(t) {
  y <- cum_validated_counts(score_mat_ord, t)
  data.frame(
    n_IGTs_included = seq_along(y),
    validated_GPs = y,
    threshold = factor(t)
  )
}) %>%
  bind_rows()

# ---- plot ----
ggplot(
  plot_df,
  aes(x = n_IGTs_included, y = validated_GPs, color = threshold)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  labs(
    x = "Number of IGTs included (in IGT index order)",
    y = "Number of validated GPs (cumulative union)",
    color = "Threshold"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0, ncol(score_mat_ord), by = 5)) +
  scale_y_continuous(breaks = seq(0, max(plot_df$validated_GPs), by = 20))

# save as pdf
ggsave(
  paste0(figure_path, "igt_specific_cumulative_validated.pdf"),
  width = 6,
  height = 4
)
