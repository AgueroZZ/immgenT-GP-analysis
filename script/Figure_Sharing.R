library(pheatmap)

#####################################################
### Paths and data
#####################################################
data_path <- "data/"
figure_path <- "figures/Figure4_Sharing/"
path_cosine <- paste0(figure_path, "cosine_similarity/")
dir.create(path_cosine, recursive = TRUE, showWarnings = FALSE)

L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
colnames(L_pm_filtered) <- paste0("GP", seq_len(ncol(L_pm_filtered)))

seurat_meta <- readRDS(paste0(
  data_path,
  "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"
))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
rm(seurat_meta)


#####################################################
### Figure 3 GPs and lineages
###
### Color groups match Figure_Activation.R:
###   CD4 only  -> blue
###   CD8 only  -> darkorange2
###   both up   -> darkred
###   both down -> darkgreen
#####################################################
gp_groups <- list(
  "CD4 only" = c(
    "GP56",
    "GP162",
    "GP36",
    "GP152",
    "GP161",
    "GP177",
    "GP35",
    "GP72",
    "GP79",
    "GP12"
  ),
  "CD8 only" = c("GP43", "GP10", "GP58"),
  "both up" = c("GP25", "GP26"),
  "both down" = c("GP9", "GP171")
)
group_colors <- c(
  "CD4 only" = "blue",
  "CD8 only" = "darkorange2",
  "both up" = "darkred",
  "both down" = "darkgreen"
)
GPs_of_interest <- unlist(gp_groups, use.names = FALSE)
lineages <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN")


#####################################################
### Per-lineage cosine similarity heatmaps (Figure 3 GPs only)
###
### A single GP order is derived by clustering the overall cosine matrix
### (all cells, 17 GPs) so all per-lineage heatmaps share axes and can be
### compared directly.
#####################################################
L_subset <- L_pm_filtered[, GPs_of_interest, drop = FALSE]

cosine_sim <- function(mat) {
  norms <- sqrt(colSums(mat^2))
  norms[norms == 0] <- 1
  crossprod(mat) / outer(norms, norms)
}

cos_overall <- cosine_sim(L_subset)
hc_overall <- hclust(as.dist(1 - cos_overall), method = "complete")
gp_order <- hc_overall$labels[hc_overall$order]

# Per-GP group annotation (color bar matches Figure_Activation.R)
gp_to_group <- setNames(
  rep(names(gp_groups), lengths(gp_groups)),
  unlist(gp_groups, use.names = FALSE)
)
annotation_df <- data.frame(
  Group = factor(gp_to_group[gp_order], levels = names(gp_groups)),
  row.names = gp_order
)
annotation_colors <- list(Group = group_colors)

breaks_cosine <- seq(0, 1, length.out = 201)
cols_cosine <- colorRampPalette(c("white", "red"))(length(breaks_cosine) - 1)

for (lin in lineages) {
  cell_idx <- which(seurat_meta_filtered$annotation_level1 == lin)
  if (length(cell_idx) < 10) {
    next
  }
  cos_lin <- cosine_sim(L_subset[cell_idx, , drop = FALSE])[gp_order, gp_order]
  pheatmap(
    cos_lin,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = cols_cosine,
    breaks = breaks_cosine,
    annotation_row = annotation_df,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    main = paste0("GP Cosine Similarity — ", lin, " cells"),
    filename = paste0(path_cosine, "Cosine_Sim_Figure3_GPs_", lin, ".pdf"),
    width = 7.5,
    height = 7
  )
}

#####################################################
### Heatmap showing average GP loadings for selected GPs
###
### Rows: the 17 Figure 3 GPs in semantic group order.
### Cols: annotation_level2 sub-types (parent in `lineages`, >= 50 cells),
###       ordered by parent lineage.
### Color groups (row + column annotation bars) match Figure_Activation.R.
#####################################################
keep_cells <- seurat_meta_filtered$annotation_level1 %in% lineages
meta_sub   <- seurat_meta_filtered[keep_cells, ]
L_keep     <- L_subset[keep_cells, , drop = FALSE]

l2_counts <- table(meta_sub$annotation_level2)
l2_keep   <- names(l2_counts)[l2_counts >= 50]

# Drop the "P" cluster and any "w..." clusters (wM, wW, etc.) across all
# lineages — matches the filter used in Figure_CD4.R. Strips an
# arbitrary lineage prefix ("CD4.", "CD8_", "Treg.", ...) before testing.
l2_stripped <- sub("^[^._]+[._]", "", l2_keep)
exclude_l2  <- l2_stripped == "P" |
  grepl("^w", l2_stripped, ignore.case = TRUE) |
  grepl("[._]w", l2_keep, ignore.case = TRUE)
l2_keep <- l2_keep[!exclude_l2]

mean_mat <- vapply(l2_keep, function(l2) {
  colMeans(L_keep[meta_sub$annotation_level2 == l2, , drop = FALSE])
}, numeric(length(GPs_of_interest)))

# Map each level2 column to its parent level1, then order columns by lineage
l2_to_l1  <- vapply(l2_keep, function(l2) {
  as.character(meta_sub$annotation_level1[meta_sub$annotation_level2 == l2][1])
}, character(1))
col_order <- order(match(l2_to_l1, lineages), l2_keep)
mean_mat  <- mean_mat[GPs_of_interest, col_order]
l2_to_l1  <- l2_to_l1[col_order]

immgen_cols <- ZemmourLib::immgent_colors

# Column annotation bar: parent level1 (uses ZemmourLib$level1)
col_anno <- data.frame(
  Lineage   = factor(l2_to_l1, levels = lineages),
  row.names = colnames(mean_mat)
)
anno_colors_mean <- list(Lineage = immgen_cols$level1[lineages])

# Per-label text colors
#   rows: GP group (matches Figure_Activation.R)
#   cols: level2 colors from ZemmourLib (fall back to black if missing)
row_label_cols <- group_colors[gp_to_group[GPs_of_interest]]
col_label_cols <- immgen_cols$level2[colnames(mean_mat)]
col_label_cols[is.na(col_label_cols)] <- "black"

ph <- pheatmap(
  mean_mat,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  color             = colorRampPalette(c("white", "red"))(200),
  annotation_col    = col_anno,
  annotation_colors = anno_colors_mean,
  gaps_row          = head(cumsum(lengths(gp_groups)), -1),
  gaps_col          = head(cumsum(rle(l2_to_l1)$lengths), -1),
  main              = "Average loading of Figure 3 GPs per Level-2 sub-lineage",
  silent            = TRUE
)

# Recolor row + column label text (pheatmap doesn't expose per-label colors)
row_idx <- which(ph$gtable$layout$name == "row_names")
col_idx <- which(ph$gtable$layout$name == "col_names")
ph$gtable$grobs[[row_idx]]$gp$col <- row_label_cols
ph$gtable$grobs[[col_idx]]$gp$col <- col_label_cols

pdf(
  paste0(figure_path, "MeanLoading_Figure3_GPs_by_Level2.pdf"),
  width = 11, height = 5.5
)
grid::grid.draw(ph$gtable)
dev.off()
