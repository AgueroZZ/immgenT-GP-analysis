compute_auc_by_group <- function(
    loading,
    group_info,
    groups = NULL,           # specify which groups to evaluate; defaults to all
    sort_groups = TRUE,      # whether to sort group names alphabetically
    remove_na = TRUE,        # whether to remove missing values
    ci = FALSE,              # whether to compute confidence intervals for AUC
    ci_level = 0.95,         # confidence level for AUC CI
    ci_boot_n = 2000,        # number of bootstrap replicates for CI
    quiet = TRUE,            # suppress messages from pROC
    direction = "auto",      # ROC direction ("auto" lets pROC decide)
    ...
) {
  # Check dependency
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required. Please install it with install.packages('pROC').")
  }

  # Sanity check
  if (length(loading) != length(group_info)) {
    stop("`loading` and `group_info` must have the same length.")
  }

  # Handle missing values
  if (remove_na) {
    keep <- !(is.na(loading) | is.na(group_info))
    loading <- loading[keep]
    group_info <- group_info[keep]
  }

  # Determine the set of groups to evaluate
  if (is.null(groups)) {
    unique_groups <- unique(group_info)
  } else {
    unique_groups <- intersect(groups, unique(group_info))
  }
  if (sort_groups) unique_groups <- sort(unique_groups)

  # Compute AUC for each group (one-vs-rest)
  res_list <- lapply(unique_groups, function(g) {
    labels <- as.numeric(group_info == g)  # 1 = positive class, 0 = negative class

    # Skip if there's only one class
    if (length(unique(labels)) < 2L) {
      return(data.frame(
        group = as.character(g),
        auc = NA_real_,
        auc_lower = NA_real_,
        auc_upper = NA_real_,
        n_pos = sum(labels == 1L),
        n_neg = sum(labels == 0L),
        stringsAsFactors = FALSE
      ))
    }

    # Compute ROC and AUC
    roc_obj <- pROC::roc(
      response = labels,
      predictor = loading,
      quiet = quiet,
      direction = direction,
      levels = c(0, 1),
      ...
    )

    auc_val <- as.numeric(pROC::auc(roc_obj))

    # Optionally compute bootstrap confidence interval
    if (ci) {
      ci_vec <- suppressMessages(
        pROC::ci.auc(roc_obj, conf.level = ci_level, boot.n = ci_boot_n)
      )
      auc_lower <- as.numeric(ci_vec[1])
      auc_upper <- as.numeric(ci_vec[3])
    } else {
      auc_lower <- NA_real_
      auc_upper <- NA_real_
    }

    data.frame(
      group = as.character(g),
      auc = auc_val,
      auc_lower = auc_lower,
      auc_upper = auc_upper,
      n_pos = sum(labels == 1L),
      n_neg = sum(labels == 0L),
      stringsAsFactors = FALSE
    )
  })

  # Combine results into a single data.frame
  do.call(rbind, res_list)
}

compute_auc_matrix <- function(loading_mat, group_info, sort_groups = TRUE, remove_na = TRUE, quiet = TRUE) {
  # Check input
  if (!is.matrix(loading_mat) && !is.data.frame(loading_mat)) {
    stop("`loading_mat` must be a matrix or data.frame.")
  }

  # Set up column names
  n_loadings <- ncol(loading_mat)
  if (is.null(colnames(loading_mat)) || any(colnames(loading_mat) == "")) {
    colnames(loading_mat) <- paste0("K", seq_len(n_loadings))
  }

  # Compute AUCs per column
  auc_list <- lapply(seq_len(n_loadings), function(j) {
    res <- compute_auc_by_group(
      loading = loading_mat[, j],
      group_info = group_info,
      sort_groups = sort_groups,
      remove_na = remove_na,
      quiet = quiet
    )
    res
  })

  # Combine AUC columns
  auc_mat <- do.call(cbind, lapply(auc_list, function(x) x$auc))

  # Add row and column names
  rownames(auc_mat) <- names(auc_list[[1]]$auc)
  colnames(auc_mat) <- colnames(loading_mat)

  return(auc_mat)
}

# Helper: per-loading, per-group AUC + threshold
compute_auc_threshold_by_group <- function(
    loading,
    group_info,
    groups      = NULL,
    sort_groups = TRUE,
    remove_na   = TRUE,
    quiet       = TRUE,
    direction   = "auto",
    best_method = c("closest.topleft", "youden"),
    ...
) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required. Please install it with install.packages('pROC').")
  }

  best_method <- match.arg(best_method)

  # Sanity check
  if (length(loading) != length(group_info)) {
    stop("`loading` and `group_info` must have the same length.")
  }

  # Handle missing values
  if (remove_na) {
    keep <- !(is.na(loading) | is.na(group_info))
    loading    <- loading[keep]
    group_info <- group_info[keep]
  }

  # Determine groups
  if (is.null(groups)) {
    unique_groups <- unique(group_info)
  } else {
    unique_groups <- intersect(groups, unique(group_info))
  }
  if (sort_groups) unique_groups <- sort(unique_groups)

  # One-vs-rest ROC per group
  res_list <- lapply(unique_groups, function(g) {
    labels <- as.numeric(group_info == g)  # 1 = positive, 0 = negative

    if (length(unique(labels)) < 2L) {
      return(data.frame(
        group     = as.character(g),
        auc       = NA_real_,
        threshold = NA_real_,
        n_pos     = sum(labels == 1L),
        n_neg     = sum(labels == 0L),
        stringsAsFactors = FALSE
      ))
    }

    roc_obj <- pROC::roc(
      response  = labels,
      predictor = loading,
      quiet     = quiet,
      direction = direction,
      levels    = c(0, 1),
      ...
    )

    auc_val <- as.numeric(pROC::auc(roc_obj))

    thr_val <- as.numeric(
      pROC::coords(
        roc_obj,
        x           = "best",
        best.method = best_method,
        ret         = "threshold",
        transpose   = FALSE
      )
    )

    data.frame(
      group     = as.character(g),
      auc       = auc_val,
      threshold = thr_val,
      n_pos     = sum(labels == 1L),
      n_neg     = sum(labels == 0L),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, res_list)
}

# Main: return AUC matrix + threshold matrix
compute_auc_threshold_matrix <- function(
    loading_mat,
    group_info,
    sort_groups = TRUE,
    remove_na   = TRUE,
    quiet       = TRUE,
    direction   = "auto",
    best_method = c("closest.topleft", "youden"),
    ...
) {
  # Check input
  if (!is.matrix(loading_mat) && !is.data.frame(loading_mat)) {
    stop("`loading_mat` must be a matrix or data.frame.")
  }

  loading_mat <- as.matrix(loading_mat)
  n_loadings  <- ncol(loading_mat)

  # Ensure column names
  if (is.null(colnames(loading_mat)) || any(colnames(loading_mat) == "")) {
    colnames(loading_mat) <- paste0("K", seq_len(n_loadings))
  }

  # Compute per-column AUC + threshold
  res_list <- lapply(seq_len(n_loadings), function(j) {
    compute_auc_threshold_by_group(
      loading    = loading_mat[, j],
      group_info = group_info,
      sort_groups = sort_groups,
      remove_na   = remove_na,
      quiet       = quiet,
      direction   = direction,
      best_method = best_method,
      ...
    )
  })

  # Build matrices
  auc_mat <- do.call(cbind, lapply(res_list, function(x) x$auc))
  thr_mat <- do.call(cbind, lapply(res_list, function(x) x$threshold))

  rownames(auc_mat) <- res_list[[1]]$group
  rownames(thr_mat) <- res_list[[1]]$group
  colnames(auc_mat) <- colnames(loading_mat)
  colnames(thr_mat) <- colnames(loading_mat)

  list(
    auc       = auc_mat,
    threshold = thr_mat
  )
}

