### Re-run AUC analysis
library(ggplot2)
library(fastTopics)
library(dplyr)
library(ggrepel)
library(tidyr)
library(tibble)
library(cowplot)
library(patchwork)

data_path <- "/project2/mstephens/immgent/"
# code_path <- "code/"
code_path <- ""
source(paste0(code_path, "ROC.R"))
seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))

cells_flashier <- rownames(flashier_snmf_summary$L_pm)
cells_seurat <- seurat_meta$cellID
L_pm <- flashier_snmf_summary$L_pm[cells_flashier %in% cells_seurat, ]
colnames(L_pm) <- paste0("K", 1:ncol(L_pm))

# For the AUC analysis, let's do not consider thymocytes
cells_thymocyte <- which(seurat_meta$annotation_level1 == "thymocyte")
L_pm_no_thymocytes <- L_pm[-cells_thymocyte, ]
seurat_meta_no_thymocytes <- seurat_meta[-cells_thymocyte, ]


# Compute the level 2 AUC within each level 1 cell type
compute_auc_threshold_by_group_within_strata <- function(
    loading,
    group_info,   # level-2 annotation
    strata_info,  # level-1 annotation
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

  if (length(loading) != length(group_info) || length(loading) != length(strata_info)) {
    stop("`loading`, `group_info`, and `strata_info` must have the same length.")
  }

  if (remove_na) {
    keep <- !(is.na(loading) | is.na(group_info) | is.na(strata_info))
    loading     <- loading[keep]
    group_info  <- group_info[keep]
    strata_info <- strata_info[keep]
  }

  if (is.null(groups)) {
    unique_groups <- unique(group_info)
  } else {
    unique_groups <- intersect(groups, unique(group_info))
  }
  if (sort_groups) unique_groups <- sort(unique_groups)

  res_list <- lapply(unique_groups, function(g) {
    parent_level1 <- unique(strata_info[group_info == g])
    parent_level1 <- parent_level1[!is.na(parent_level1)]

    if (length(parent_level1) != 1L) {
      warning(sprintf(
        "Group %s has %d parent strata (after removing NAs). Skipping.",
        as.character(g), length(parent_level1)
      ))
      return(data.frame(
        group      = as.character(g),
        stratum    = NA_character_,
        auc        = NA_real_,
        threshold  = NA_real_,
        n_pos      = sum(group_info == g),
        n_neg      = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }

    idx_stratum <- which(strata_info == parent_level1)

    labels <- as.numeric(group_info[idx_stratum] == g)  # 1 = positive, 0 = negative

    if (length(unique(labels)) < 2L) {
      return(data.frame(
        group      = as.character(g),
        stratum    = as.character(parent_level1),
        auc        = NA_real_,
        threshold  = NA_real_,
        n_pos      = sum(labels == 1L),
        n_neg      = sum(labels == 0L),
        stringsAsFactors = FALSE
      ))
    }

    roc_obj <- pROC::roc(
      response  = labels,
      predictor = loading[idx_stratum],
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
      group      = as.character(g),
      stratum    = as.character(parent_level1),
      auc        = auc_val,
      threshold  = thr_val,
      n_pos      = sum(labels == 1L),
      n_neg      = sum(labels == 0L),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, res_list)
}



compute_auc_threshold_matrix_within_strata <- function(
    loading_mat,
    group_info,   # level-2
    strata_info,  # level-1
    sort_groups = TRUE,
    remove_na   = TRUE,
    quiet       = TRUE,
    direction   = "auto",
    best_method = c("closest.topleft", "youden"),
    ...
) {
  if (!is.matrix(loading_mat) && !is.data.frame(loading_mat)) {
    stop("`loading_mat` must be a matrix or data.frame.")
  }

  loading_mat <- as.matrix(loading_mat)
  n_loadings  <- ncol(loading_mat)

  if (length(group_info) != nrow(loading_mat) ||
      length(strata_info) != nrow(loading_mat)) {
    stop("`group_info` and `strata_info` must have length nrow(loading_mat).")
  }

  if (is.null(colnames(loading_mat)) || any(colnames(loading_mat) == "")) {
    colnames(loading_mat) <- paste0("K", seq_len(n_loadings))
  }

  res_list <- lapply(seq_len(n_loadings), function(j) {
    compute_auc_threshold_by_group_within_strata(
      loading     = loading_mat[, j],
      group_info  = group_info,
      strata_info = strata_info,
      sort_groups = sort_groups,
      remove_na   = remove_na,
      quiet       = quiet,
      direction   = direction,
      best_method = best_method,
      ...
    )
  })

  auc_mat <- do.call(cbind, lapply(res_list, function(x) x$auc))
  thr_mat <- do.call(cbind, lapply(res_list, function(x) x$threshold))

  rownames(auc_mat) <- res_list[[1]]$group
  rownames(thr_mat) <- res_list[[1]]$group
  colnames(auc_mat) <- colnames(loading_mat)
  colnames(thr_mat) <- colnames(loading_mat)

  list(
    auc       = auc_mat,
    threshold = thr_mat,
    meta      = res_list[[1]][, c("group", "stratum", "n_pos", "n_neg")]
  )
}

level_2_within_level1_AUC_list <- compute_auc_threshold_matrix_within_strata(
  loading_mat = L_pm_no_thymocytes,
  group_info  = seurat_meta_no_thymocytes$annotation_level2,
  strata_info = seurat_meta_no_thymocytes$annotation_level1,
  best_method = "closest.topleft"
)

saveRDS(level_2_within_level1_AUC_list, file = paste0(data_path, "level_2_within_level_1_AUC_list_figure.rds"))

