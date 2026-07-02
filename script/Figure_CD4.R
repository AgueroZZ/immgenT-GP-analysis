library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(dplyr)
library(scales)
library(scattermore)

#####################################################
#####################################################
#####################################################
##### Defining directory and loading functions
#####################################################
#####################################################
data_path <- "data/"
code_path <- "code/"
figure_path <- "figures/Figure_CD4/"
dir.create(figure_path, showWarnings = FALSE, recursive = TRUE)

#####################################################
#####################################################
#####################################################
### Loading Data
#####################################################
#####################################################
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
colnames(L_pm_filtered) <- paste0("GP", seq_len(ncol(L_pm_filtered)))
seurat_meta <- readRDS(paste0(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]


#####################################################
#####################################################
#####################################################
### CD4 helper GPs
#####################################################
#####################################################
# Map of GPs of interest in CD4 to their helper-cell annotation
cd4_gp_map <- c(
  "GP56" = "Th2",
  "GP36" = "Th17",
  "GP12" = "Tfh",
  "GP10" = "Th1"
)

gp_label <- function(gp) {
  paste0(gp, " (", cd4_gp_map[[gp]], ")")
}


#####################################################
#####################################################
#####################################################
### Scatter plot of two GP loadings within CD4,
### with per-level2 summary dots on top.
#####################################################
#####################################################
# `summary_fun` controls the per-level2 summary point: pass `mean`, `median`,
# or any other vector-summarizing function. Change it to switch between them.
plot_cd4_gp_scatter <- function(
  gp_x,
  gp_y,
  L,
  meta,
  summary_fun = mean,
  transform = c("raw", "log2fc"),
  threshold = NULL,
  pc = 1e-10,
  cap = NULL,
  bg_size = 0.6,
  bg_alpha = 0.5,
  bg_color = "darkorange2",
  point_size = 4,
  point_stroke = 0.8,
  label_size = 3.5,
  pct_size = 3.2,
  seed = 1
) {
  stopifnot(gp_x %in% colnames(L), gp_y %in% colnames(L))
  transform <- match.arg(transform)
  if (is.null(threshold)) {
    threshold <- if (transform == "log2fc") 0 else 0.1
  }

  cd4_cells <- meta$cellID[meta$annotation_level1 == "CD4"]
  cd4_cells <- intersect(cd4_cells, rownames(L))

  level2 <- as.character(meta$annotation_level2[match(cd4_cells, meta$cellID)])
  # Strip the "CD4." or "CD4_" prefix (both separators appear in the meta).
  level2_label <- sub("^CD4[._]", "", level2)

  # Drop the "P" cluster and any "w..." clusters (wM, wW, etc.).
  # Case-insensitive on the stripped label; also matches a raw "CD4.wM"
  # or "CD4_wT" defensively in case the strip above missed a separator.
  excluded <- !is.na(level2_label) &
    (
      level2_label == "P" |
        grepl("^w", level2_label, ignore.case = TRUE) |
        grepl("[._]w", level2, ignore.case = TRUE)
    )
  cd4_cells <- cd4_cells[!excluded]
  level2 <- level2[!excluded]
  level2_label <- level2_label[!excluded]

  x_raw <- L[cd4_cells, gp_x]
  y_raw <- L[cd4_cells, gp_y]

  if (transform == "log2fc") {
    # log2 fold change of each cell's loading vs. the CD4-overall mean
    # for that GP. Pseudocount `pc` avoids log(0) on sparse loadings.
    mu_x <- mean(x_raw, na.rm = TRUE)
    mu_y <- mean(y_raw, na.rm = TRUE)
    x_val <- log2((x_raw + pc) / (mu_x + pc))
    y_val <- log2((y_raw + pc) / (mu_y + pc))
    if (!is.null(cap)) {
      x_val <- pmax(pmin(x_val, cap), -cap)
      y_val <- pmax(pmin(y_val, cap), -cap)
    }
    axis_suffix <- " log2FC vs CD4 mean"
  } else {
    x_val <- x_raw
    y_val <- y_raw
    axis_suffix <- " loading"
  }

  df_cells <- data.frame(
    x = x_val,
    y = y_val,
    level2 = level2_label,
    level2_full = level2,
    stringsAsFactors = FALSE
  )

  # Per-level2 summary uses only cells with a valid level 2 annotation;
  # background scatter shows all CD4 cells regardless.
  df_summary <- df_cells %>%
    filter(!is.na(level2) & !tolower(level2) %in% c("", "na", "nan")) %>%
    group_by(level2) %>%
    summarise(
      x = summary_fun(x, na.rm = TRUE),
      y = summary_fun(y, na.rm = TRUE),
      n = dplyr::n(),
      level2_full = dplyr::first(level2_full),
      .groups = "drop"
    )

  # Canonical level 2 palette. Try both the full "CD4.xxx" key and the
  # stripped label, since the palette is sometimes keyed by either.
  # Single lookup (by stripped label) is reused for both the per-cell
  # background dots and the per-level2 summary dots.
  level2_palette <- tryCatch(
    ZemmourLib::immgent_colors$level2,
    error = function(e) NULL
  )
  unique_l2 <- unique(df_cells$level2)
  unique_l2_full <- df_cells$level2_full[match(unique_l2, df_cells$level2)]
  resolve_palette <- function(full, stripped) {
    if (is.null(level2_palette)) return(NA_character_)
    if (full %in% names(level2_palette)) return(level2_palette[[full]])
    if (stripped %in% names(level2_palette)) return(level2_palette[[stripped]])
    NA_character_
  }
  level2_color_lookup <- setNames(
    vapply(
      seq_along(unique_l2),
      function(i) resolve_palette(unique_l2_full[i], unique_l2[i]),
      character(1)
    ),
    unique_l2
  )
  # Neutral fallback for any level 2 not in the palette (not bg_color, so
  # missing-from-palette is visually distinguishable from a real assignment).
  missing_in_palette <- is.na(level2_color_lookup)
  if (any(missing_in_palette)) {
    level2_color_lookup[missing_in_palette] <- "grey60"
    warning(sprintf(
      "Level 2 clusters not in ZemmourLib::immgent_colors$level2 (colored grey60): %s",
      paste(unique_l2[missing_in_palette], collapse = ", ")
    ))
  }

  # % of CD4 cells in each of the four quadrants formed by the
  # x = threshold and y = threshold dashed lines.
  pct_bl <- mean(
    df_cells$x < threshold & df_cells$y < threshold,
    na.rm = TRUE
  ) *
    100
  pct_br <- mean(
    df_cells$x >= threshold & df_cells$y < threshold,
    na.rm = TRUE
  ) *
    100
  pct_tl <- mean(
    df_cells$x < threshold & df_cells$y >= threshold,
    na.rm = TRUE
  ) *
    100
  pct_tr <- mean(
    df_cells$x >= threshold & df_cells$y >= threshold,
    na.rm = TRUE
  ) *
    100
  fmt_pct <- function(p) sprintf("%.1f%%", p)

  summary_name <- deparse(substitute(summary_fun))

  set.seed(seed)
  ggplot() +
    scattermore::geom_scattermore(
      data = df_cells,
      mapping = aes(x = x, y = y, color = level2),
      pointsize = bg_size,
      alpha = bg_alpha
    ) +
    scale_color_manual(values = level2_color_lookup, guide = "none") +
    geom_vline(
      xintercept = threshold,
      linetype = "dashed",
      color = "grey40"
    ) +
    geom_hline(
      yintercept = threshold,
      linetype = "dashed",
      color = "grey40"
    ) +
    # % of CD4 cells in each quadrant, placed at the plot corners.
    annotate(
      "text",
      x = -Inf,
      y = -Inf,
      label = fmt_pct(pct_bl),
      hjust = -0.2,
      vjust = -0.6,
      size = pct_size,
      color = "grey25"
    ) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = fmt_pct(pct_br),
      hjust = 1.2,
      vjust = -0.6,
      size = pct_size,
      color = "grey25"
    ) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = fmt_pct(pct_tl),
      hjust = -0.2,
      vjust = 1.5,
      size = pct_size,
      color = "grey25"
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = fmt_pct(pct_tr),
      hjust = 1.2,
      vjust = 1.5,
      size = pct_size,
      color = "grey25"
    ) +
    geom_point(
      data = df_summary,
      mapping = aes(x = x, y = y, fill = level2),
      shape = 21,
      color = "black",
      size = point_size,
      stroke = point_stroke
    ) +
    scale_fill_manual(values = level2_color_lookup, guide = "none") +
    ggrepel::geom_text_repel(
      data = df_summary,
      mapping = aes(x = x, y = y, label = level2),
      size = label_size,
      box.padding = 0.4,
      point.padding = 0.3,
      min.segment.length = 0,
      segment.color = "grey50",
      max.overlaps = Inf
    ) +
    labs(
      title = paste0(
        gp_label(gp_y),
        " vs ",
        gp_label(gp_x),
        " in CD4"
      ),
      subtitle = paste0(
        "Per-level2 ",
        summary_name,
        if (transform == "log2fc") " log2FC" else " loading",
        " (n level2 = ",
        nrow(df_summary),
        ", n CD4 cells = ",
        nrow(df_cells),
        ")"
      ),
      x = paste0(gp_label(gp_x), axis_suffix),
      y = paste0(gp_label(gp_y), axis_suffix),
      fill = "Level 2"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )
}


#####################################################
#####################################################
#####################################################
### Generate all pairwise scatter plots
#####################################################
#####################################################
# Set `summary_fun` to `mean` or `median` (or any vector summarizer).
# Set `transform` to "raw" for raw loadings or "log2fc" for log2 fold change
# vs. the CD4-overall mean of each GP (compresses the near-zero blob and
# emphasises per-level2 separation).
summary_fun <- mean
transform <- "log2fc"

gps <- names(cd4_gp_map)
pair_grid <- t(combn(gps, 2))

plots <- lapply(seq_len(nrow(pair_grid)), function(i) {
  plot_cd4_gp_scatter(
    gp_x = pair_grid[i, 1],
    gp_y = pair_grid[i, 2],
    L = L_pm_filtered,
    meta = seurat_meta_filtered,
    summary_fun = summary_fun,
    transform = transform
  )
})

summary_name <- deparse(substitute(summary_fun))
out_tag <- paste0(summary_name, "_", transform)
pdf(
  paste0(figure_path, "CD4_GP_scatter_pairs_", out_tag, ".pdf"),
  width = 6,
  height = 6
)
for (p in plots) {
  print(p)
}
dev.off()

# Combined grid (one page) for quick inspection
combined <- wrap_plots(plots, ncol = 3)
ggsave(
  filename = paste0(
    figure_path,
    "CD4_GP_scatter_pairs_",
    out_tag,
    "_grid.pdf"
  ),
  plot = combined,
  width = 16,
  height = 10
)
