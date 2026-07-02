# Small, generic helpers shared across multiple Figure*.R scripts.
# (Figure/topic-specific plotting functions live in their own
# code/R/<topic>_plots.R file instead of here.)

# Multiply each column j of A by b[j].
scale_cols <- function(A, b) t(t(A) * b)

# Iteratively drop cells whose total GP membership (rowSum of L) exceeds
# max_val, renormalizing columns to max 1 after each pass. Used once
# upstream (see code/pipeline/01_extract_data.R) to produce the
# cached L_pm_filtered.rds / F_pm_filtered.rds that every figure script reads.
filter_cells_by_total_membership <- function(L, max_val = 10, numiter = 10) {
  n <- nrow(L)
  rows <- 1:n
  for (iter in 1:numiter) {
    x <- rowSums(L)
    cat(sprintf("%d. Filtered out %d cells.\n", iter, sum(x > max_val)))
    i <- which(x <= max_val)
    L <- L[i, ]
    rows <- rows[i]
    d <- apply(L, 2, max)
    L <- scale_cols(L, 1 / d)
  }
  return(rows)
}

# Identify Tukey outliers per group (matches geom_boxplot's outlier definition).
# Returns rows where value falls outside [Q1 - 1.5*IQR, Q3 + 1.5*IQR] within each
# group defined by `group_cols`. Used to plot outliers as a separate rasterized
# layer (smaller PDF), with `geom_boxplot(outlier.shape = NA)` for the boxes.
tukey_outliers <- function(df, value_col, group_cols) {
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(
      .q1 = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
      .q3 = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
      .iqr = .q3 - .q1
    ) %>%
    dplyr::filter(
      .data[[value_col]] < .q1 - 1.5 * .iqr |
        .data[[value_col]] > .q3 + 1.5 * .iqr
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.q1, -.q3, -.iqr)
}

# Lineage color palette (falls back to a fixed default if ZemmourLib isn't installed).
lineage_colors <- function() {
  tryCatch(
    ZemmourLib::immgent_colors$level1,
    error = function(e) {
      c(
        CD4 = "blue",
        CD8 = "darkorange2",
        Treg = "deeppink",
        gdT = "chartreuse3",
        CD8aa = "darkorchid",
        Tz = "darkgoldenrod1",
        DN = "deepskyblue",
        DP = "red"
      )
    }
  )
}
