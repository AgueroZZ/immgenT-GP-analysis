library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(fastTopics)
library(dplyr)
library(pheatmap)
library(scales)
library(scattermore)
library(lme4)

#####################################################
#####################################################
####################################################
##### Defining directory and loading functions
####################################################
#####################################################
# data_path <- "../data/"
# code_path <- "../code/"
data_path <- "data/"
code_path <- "code/"
figure_path <- "figures/Figure_Cosmos/"
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)
seurat_meta <- readRDS(paste0(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data


# create a dataframe for sample, with column being the proportion of each level 2 subtype, and row being the sample (unique(seurat_meta$IGTHT))
seurat_meta_no_thymocytes <- seurat_meta %>%
  filter(
    !is.na(IGTHT),
    !is.na(annotation_level1),
    !is.na(annotation_level2),
    annotation_level1 != "thymocyte"
  )

sample_ids <- sort(unique(as.character(seurat_meta_no_thymocytes$IGTHT)))
level1_types <- sort(unique(as.character(
  seurat_meta_no_thymocytes$annotation_level1
)))
level2_types <- sort(unique(as.character(
  seurat_meta_no_thymocytes$annotation_level2
)))

count_mat <- table(
  IGTHT = factor(
    as.character(seurat_meta_no_thymocytes$IGTHT),
    levels = sample_ids
  ),
  annotation_level1 = factor(
    as.character(seurat_meta_no_thymocytes$annotation_level1),
    levels = level1_types
  ),
  annotation_level2 = factor(
    as.character(seurat_meta_no_thymocytes$annotation_level2),
    levels = level2_types
  )
)

sample_level2_proportion_array <- prop.table(count_mat, margin = c(1, 2))
sample_level2_proportion_array[!is.finite(sample_level2_proportion_array)] <- 0

level2_to_level1 <- seurat_meta_no_thymocytes %>%
  distinct(annotation_level2, annotation_level1) %>%
  mutate(
    annotation_level2 = as.character(annotation_level2),
    annotation_level1 = as.character(annotation_level1)
  )

sample_level2_proportion_mat <- matrix(
  0,
  nrow = length(sample_ids),
  ncol = length(level2_types),
  dimnames = list(sample_ids, level2_types)
)

for (i in seq_len(nrow(level2_to_level1))) {
  level2_type <- level2_to_level1$annotation_level2[i]
  level1_type <- level2_to_level1$annotation_level1[i]
  sample_level2_proportion_mat[, level2_type] <-
    sample_level2_proportion_array[, level1_type, level2_type]
}

sample_level2_proportion_mat[!is.finite(sample_level2_proportion_mat)] <- 0

sample_level2_proportion_df <- as.data.frame.matrix(
  sample_level2_proportion_mat
)
sample_level2_proportion_df <- data.frame(
  IGTHT = rownames(sample_level2_proportion_df),
  sample_level2_proportion_df,
  check.names = FALSE,
  row.names = NULL
)

# add the other information for each sample, such as sex, age_weeks, condition_detailed_simplified, organ_simplified
first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA)
  }
  x[1]
}

sample_metadata_df <- seurat_meta_no_thymocytes %>%
  group_by(IGTHT) %>%
  summarise(
    sex = first_non_missing(sex),
    age_weeks = first_non_missing(age_weeks),
    condition_detailed_simplified = first_non_missing(
      condition_detailed_simplified
    ),
    organ_simplified = first_non_missing(organ_simplified),
    .groups = "drop"
  )

sample_level2_proportion_df <- sample_level2_proportion_df %>%
  left_join(sample_metadata_df, by = "IGTHT") %>%
  relocate(
    sex,
    age_weeks,
    condition_detailed_simplified,
    organ_simplified,
    .after = IGTHT
  )


boxplot(
  sample_level2_proportion_df$CD4.A ~ sample_level2_proportion_df$condition_detailed_simplified
)


# a function that takes a level2 and a condition,
# it computes the variance of the level2 proportion across samples
# and variance of the mean level2 proportion across samples in each condition,
# and the proportion of variance explained by the condition (by running regression and do ANOVA)
compute_condition_variance_explained <- function(
  level2,
  condition,
  sample_df = sample_level2_proportion_df,
  model = c("fixed", "random")
) {
  model <- match.arg(model)

  if (!level2 %in% colnames(sample_df)) {
    stop("level2 is not a column in sample_df: ", level2)
  }
  if (!condition %in% colnames(sample_df)) {
    stop("condition is not a column in sample_df: ", condition)
  }

  analysis_df <- data.frame(
    level2_proportion = sample_df[[level2]],
    condition = sample_df[[condition]]
  ) %>%
    filter(!is.na(level2_proportion), !is.na(condition))

  analysis_df$condition <- factor(analysis_df$condition)

  if (nrow(analysis_df) < 2 || nlevels(analysis_df$condition) < 2) {
    return(data.frame(
      proportion_variance_explained = NA_real_,
      variance_ratio = NA_real_
    ))
  }

  if (model == "random") {
    fit_random <- tryCatch(
      suppressWarnings(
        suppressMessages(
          lmer(level2_proportion ~ 1 + (1 | condition), data = analysis_df)
        )
      ),
      error = function(e) NULL
    )

    if (is.null(fit_random)) {
      return(data.frame(
        proportion_variance_explained = NA_real_,
        variance_ratio = NA_real_
      ))
    }

    if (isSingular(fit_random)) {
      return(data.frame(
        proportion_variance_explained = 0,
        variance_ratio = 0
      ))
    }

    variance_components <- as.data.frame(VarCorr(fit_random))$vcov
    random_effect_variance <- variance_components[1]
    total_model_variance <- sum(variance_components)
    random_proportion_variance_explained <- ifelse(
      total_model_variance > 0,
      random_effect_variance / total_model_variance,
      NA_real_
    )

    return(data.frame(
      proportion_variance_explained = random_proportion_variance_explained,
      variance_ratio = random_proportion_variance_explained
    ))
  }

  condition_means <- analysis_df %>%
    group_by(condition) %>%
    summarise(
      mean_level2_proportion = mean(level2_proportion),
      .groups = "drop"
    )

  fit <- lm(level2_proportion ~ condition, data = analysis_df)
  anova_fit <- anova(fit)
  condition_sum_sq <- anova_fit["condition", "Sum Sq"]
  residual_sum_sq <- anova_fit["Residuals", "Sum Sq"]
  total_sum_sq <- condition_sum_sq + residual_sum_sq
  total_variance <- var(analysis_df$level2_proportion)
  condition_mean_variance <- var(condition_means$mean_level2_proportion)

  data.frame(
    proportion_variance_explained = condition_sum_sq / total_sum_sq,
    variance_ratio = condition_mean_variance / total_variance
  )
}

# repeatedly use this function for each level 2 and each possible condition
condition_vars <- c(
  "age_weeks",
  "sex",
  "condition_detailed_simplified",
  "organ_simplified"
)

compute_condition_variance_grid <- function(
  level2_types,
  condition_vars,
  sample_df = sample_level2_proportion_df,
  model = c("fixed", "random")
) {
  model <- match.arg(model)

  proportion_variance_explained_mat <- matrix(
    NA_real_,
    nrow = length(level2_types),
    ncol = length(condition_vars),
    dimnames = list(level2_types, condition_vars)
  )

  variance_ratio_mat <- matrix(
    NA_real_,
    nrow = length(level2_types),
    ncol = length(condition_vars),
    dimnames = list(level2_types, condition_vars)
  )

  for (level2 in level2_types) {
    for (condition in condition_vars) {
      variance_result <- compute_condition_variance_explained(
        level2 = level2,
        condition = condition,
        sample_df = sample_df,
        model = model
      )

      proportion_variance_explained_mat[level2, condition] <-
        variance_result$proportion_variance_explained
      variance_ratio_mat[level2, condition] <-
        variance_result$variance_ratio
    }
  }

  list(
    proportion_variance_explained_mat = proportion_variance_explained_mat,
    variance_ratio_mat = variance_ratio_mat,
    proportion_variance_explained_df = data.frame(
      annotation_level2 = rownames(proportion_variance_explained_mat),
      as.data.frame.matrix(proportion_variance_explained_mat),
      check.names = FALSE,
      row.names = NULL
    ),
    variance_ratio_df = data.frame(
      annotation_level2 = rownames(variance_ratio_mat),
      as.data.frame.matrix(variance_ratio_mat),
      check.names = FALSE,
      row.names = NULL
    )
  )
}

plot_condition_variance_heatmap <- function(
  proportion_variance_explained_df,
  output_file,
  plot_title
) {
  heatmap_df <- proportion_variance_explained_df %>%
    left_join(level2_to_level1, by = "annotation_level2") %>%
    arrange(factor(annotation_level1, levels = level1_types), annotation_level2)

  heatmap_mat <- as.matrix(heatmap_df[, condition_vars, drop = FALSE])
  rownames(heatmap_mat) <- heatmap_df$annotation_level2

  row_annotation_df <- data.frame(
    annotation_level1 = factor(
      heatmap_df$annotation_level1,
      levels = level1_types
    ),
    row.names = heatmap_df$annotation_level2
  )

  pheatmap(
    heatmap_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_annotation_df,
    color = colorRampPalette(c("white", "firebrick3"))(100),
    main = plot_title,
    filename = output_file,
    width = 5,
    height = 15
  )

  invisible(heatmap_mat)
}

fixed_effect_variance_results <- compute_condition_variance_grid(
  level2_types = level2_types,
  condition_vars = condition_vars,
  model = "fixed"
)

random_effect_variance_results <- compute_condition_variance_grid(
  level2_types = level2_types,
  condition_vars = condition_vars,
  model = "random"
)

proportion_variance_explained_mat <-
  fixed_effect_variance_results$proportion_variance_explained_mat
variance_ratio_mat <- fixed_effect_variance_results$variance_ratio_mat
proportion_variance_explained_df <-
  fixed_effect_variance_results$proportion_variance_explained_df
variance_ratio_df <- fixed_effect_variance_results$variance_ratio_df

proportion_variance_explained_random_mat <-
  random_effect_variance_results$proportion_variance_explained_mat
variance_ratio_random_mat <- random_effect_variance_results$variance_ratio_mat
proportion_variance_explained_random_df <-
  random_effect_variance_results$proportion_variance_explained_df
variance_ratio_random_df <- random_effect_variance_results$variance_ratio_df

proportion_variance_explained_difference_df <- data.frame(
  annotation_level2 = level2_types,
  as.data.frame.matrix(
    proportion_variance_explained_random_mat -
      proportion_variance_explained_mat
  ),
  check.names = FALSE,
  row.names = NULL
)

proportion_variance_explained_heatmap_mat <- plot_condition_variance_heatmap(
  proportion_variance_explained_df = proportion_variance_explained_df,
  output_file = paste0(
    figure_path,
    "proportion_variance_explained_by_condition_heatmap.pdf"
  ),
  plot_title = "Fixed-effect proportion of variance explained"
)

proportion_variance_explained_random_heatmap_mat <- plot_condition_variance_heatmap(
  proportion_variance_explained_df = proportion_variance_explained_random_df,
  output_file = paste0(
    figure_path,
    "proportion_variance_explained_random_by_condition_heatmap.pdf"
  ),
  plot_title = "Random-effect proportion of variance explained"
)
