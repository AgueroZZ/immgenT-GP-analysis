## Reproduce the OLD Figure 1C: giant heatmap of all 200 GP loadings across a
## stratified sample of cells (rows = cells ordered by lineage x organ, columns =
## GPs clustered by similarity). Ported verbatim from the pre-network Figure1.R.
##
## ADDITION: the same GPs highlighted in the GP-gene network (Fig 1C) are marked
## here too -- each highlighted GP column gets its network color via (1) a top
## annotation bar, (2) a colored/bold column label, and (3) a box drawn around the
## column ("circling" the GP). GPs are COLUMNS in this heatmap, so a GP is marked
## by its column, not a row.
##
## Run from repo root:
##   Rscript experiments/gene_correlation_network/heatmap_loading.R

suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(dplyr); library(grid)
})
data_path <- "data/"; out_dir <- "experiments/gene_correlation_network/"

# same GP -> color map as the GP-gene network (gp_highlight_selected.R)
GP_HIGHLIGHTS <- c(
  GP68  = "pink2",  GP58  = "orange2", GP35  = "purple", GP171 = "blue",
  GP1   = "cyan2",  GP56  = "red2",    GP161 = "brown",  GP6   = "green2",
  GP7   = "green3", GP196 = "yellow3")

# ── data (mirrors script/Figure1.R setup, restricted to non-thymocytes) ────────
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
seurat_meta <- readRDS(paste0(data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
non_thymo <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 != "thymocyte"]
L_pm_filtered <- L_pm_filtered[non_thymo, ]
seurat_meta_filtered <- seurat_meta_filtered[non_thymo, ]

# ── old 1C heatmap construction (verbatim) ─────────────────────────────────────
set.seed(6173)
MIN_CELLS <- 20; N_SAMPLE <- 80; TOP_ORGANS <- 5; K_ANCHOR <- 5
level1_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN")

organs_all <- as.character(unique(seurat_meta_filtered$organ_simplified))
ln_match <- organs_all[grepl("^LN$|lymph", organs_all, ignore.case = TRUE)]
spleen_match <- organs_all[grepl("spleen", organs_all, ignore.case = TRUE)]
organ_order <- c(spleen_match, ln_match, sort(setdiff(organs_all, c(ln_match, spleen_match))))

top_organ_combos <- seurat_meta_filtered |>
  dplyr::filter(annotation_level1 %in% level1_order) |>
  dplyr::count(annotation_level1, organ_simplified) |>
  dplyr::group_by(annotation_level1) |>
  dplyr::slice_max(n, n = TOP_ORGANS, with_ties = FALSE) |>
  dplyr::ungroup() |> dplyr::select(annotation_level1, organ_simplified)

sampled_random <- seurat_meta_filtered |>
  dplyr::filter(annotation_level1 %in% level1_order) |>
  dplyr::inner_join(top_organ_combos, by = c("annotation_level1", "organ_simplified")) |>
  dplyr::group_by(annotation_level1, organ_simplified) |>
  dplyr::filter(dplyr::n() >= MIN_CELLS) |>
  dplyr::slice_sample(n = N_SAMPLE) |> dplyr::ungroup()

anchor_cellids <- apply(L_pm_filtered, 2, function(x)
  rownames(L_pm_filtered)[order(x, decreasing = TRUE)[seq_len(K_ANCHOR)]]) |>
  as.vector() |> unique()
anchor_meta <- seurat_meta_filtered |>
  dplyr::filter(cellID %in% anchor_cellids, annotation_level1 %in% level1_order) |>
  dplyr::inner_join(top_organ_combos, by = c("annotation_level1", "organ_simplified"))

all_meta <- dplyr::bind_rows(sampled_random, anchor_meta) |>
  dplyr::distinct(cellID, .keep_all = TRUE) |>
  dplyr::arrange(factor(annotation_level1, levels = level1_order),
                 factor(organ_simplified, levels = organ_order))

L_sampled <- L_pm_filtered[all_meta$cellID, ]
clip_val <- quantile(L_sampled, 0.99)
L_display <- pmin(L_sampled, clip_val)
colnames(L_display) <- gsub("^K", "GP", colnames(L_display))

col_fun <- colorRamp2(c(0, clip_val / 2, clip_val), c("white", "#4393c3", "#08306b"))
level1_colors <- ZemmourLib::immgent_colors$level1
organ_colors <- ZemmourLib::immgent_colors$organ_simplified
level1_present <- intersect(level1_order, as.character(unique(all_meta$annotation_level1)))
organ_present <- intersect(organ_order, as.character(unique(all_meta$organ_simplified)))

row_ann <- rowAnnotation(
  Cell_Type = factor(as.character(all_meta$annotation_level1), levels = level1_present),
  Organ = factor(as.character(all_meta$organ_simplified), levels = organ_present),
  col = list(Cell_Type = level1_colors[level1_present], Organ = organ_colors[organ_present]),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(Cell_Type = list(title = "Cell Type"), Organ = list(title = "Organ")))

# ── GP highlighting (columns): top annotation + colored bold labels ────────────
gpn      <- colnames(L_display)
hl_val   <- ifelse(gpn %in% names(GP_HIGHLIGHTS), gpn, NA_character_)
lab_col  <- ifelse(gpn %in% names(GP_HIGHLIGHTS), GP_HIGHLIGHTS[gpn], "grey55")
lab_fs   <- ifelse(gpn %in% names(GP_HIGHLIGHTS), 9, 4)
lab_face <- ifelse(gpn %in% names(GP_HIGHLIGHTS), 2, 1)

top_ann <- HeatmapAnnotation(
  `Highlighted GP` = hl_val,
  col = list(`Highlighted GP` = GP_HIGHLIGHTS),
  na_col = "white", simple_anno_size = unit(4, "mm"),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(`Highlighted GP` = list(title = "Highlighted GP",
    at = names(GP_HIGHLIGHTS), labels = names(GP_HIGHLIGHTS))))

ht <- Heatmap(
  L_display, name = "Loading", col = col_fun, left_annotation = row_ann,
  top_annotation = top_ann,
  cluster_rows = FALSE, cluster_columns = TRUE,
  clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  column_names_gp = gpar(col = lab_col, fontsize = lab_fs, fontface = lab_face),
  column_title = "Gene Programs (GPs)", column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  use_raster = TRUE, raster_quality = 3, border = FALSE,
  heatmap_legend_param = list(title = "Loading", direction = "vertical"))

# draw + box the highlighted GP columns ("circle" each GP)
render <- function() {
  ht_drawn <- draw(ht, merge_legend = TRUE)
  co <- column_order(ht_drawn)
  disp <- colnames(L_display)[co]; n <- length(co)
  decorate_heatmap_body("Loading", {
    for (gp in names(GP_HIGHLIGHTS)) {
      j <- which(disp == gp)
      if (length(j)) grid.rect(x = (j - 0.5) / n, y = 0.5, width = 1.4 / n, height = 1,
                               gp = gpar(col = GP_HIGHLIGHTS[gp], fill = NA, lwd = 2.5))
    }
  })
}
pdf(paste0(out_dir, "heatmap_loading.pdf"), width = 15, height = 20, useDingbats = FALSE)
render(); invisible(dev.off())
png(paste0(out_dir, "heatmap_loading.png"), width = 15, height = 20, units = "in", res = 130)
render(); invisible(dev.off())
cat("highlighted GP columns:", paste(intersect(names(GP_HIGHLIGHTS), gpn), collapse = ", "), "\n")
cat("saved heatmap_loading.{pdf,png}\n")
