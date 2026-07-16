# Figure 1. Overview.
#
# Panels produced (final files figures/final-selected/bits/Figure 1/1A.pdf .. 1I.pdf;
# no separate caption file was found for Figure 1 during this refactor):
#   1A  Global MDE colored by major lineage (subsampled per lineage).
#   1B  Hand-finished schematic (Adobe Illustrator) -- NOT code-generated,
#       no source to port. Not produced by this script.
#   1C  Gene-program network: 200 GPs linked by shared top signature genes.
#   1D  Giant loading heatmap (200 GP loadings x a stratified cell sample), with
#       the GPs highlighted in the 1C network marked by a top color bar, colored
#       column labels, and a box around each highlighted GP column.
# The former Figure 1 panels 1D-1I (GP1 signature volcano, active-cell/active-gene
# histograms, active-GP-count boxplots, and the CD44 scatter) have moved to
# Figure 2 (script/Figure2.R, panels 2A-2F).
#
# Source: ported from Figure_Overview.R (panels C, E-I) and
# Figure_Lineage.R's "MDE by Lineage" section (panel A, confirmed by
# an exact byte match between figures/final-selected/bits/Figure 1/1A.pdf
# and figures/Figure2_Lineage/UMAP_level1_group.pdf during this refactor).
#
# Required inputs (data/) -- see code/README.md's "Data provenance" table
# for the full picture:
#   flashier_snmf_summary.rds                [code/pipeline/01_extract_data.R]
#   L_pm_filtered.rds, F_pm_filtered.rds      [code/pipeline/01b_filter_cells.R]
#   igt1_96_..._ADTonly.Rds                   [primary input Seurat object]
#   protein_mat_normalized_lognorm.rds        [code/other/prepare_citeseq_protein_matrices_20260206.R]
#   umap_result.rds                           [gap, no producer script here]
#   mean_shifted_log_expr.rds                 [gap, no producer script here]
#   flashier_snmf_fitted_prior.rda            [gap, no producer script here]

library(ggplot2)
library(dplyr)
library(scattermore)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(Matrix) # protein_mat_normalized_lognorm is a dgCMatrix; must be
                # attached (not just loaded) for `[` subsetting to dispatch

data_path <- "data/"
figure_path <- "figures/generated/Figure 1/"
source("code/R/volcano_helpers.R") # plot_gp_signature_volcano() for panel 1D

# ============================================================
# Load data
# ============================================================
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
protein_mat_normalized_lognorm <- readRDS(paste0(data_path, "protein_mat_normalized_lognorm.rds"))
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[rownames(L_pm_filtered), "CD44"]
mde_result <- readRDS(paste0(data_path, "umap_result.rds"))
colnames(mde_result) <- c("MDE_1", "MDE_2")
mde_result <- mde_result[rownames(L_pm_filtered), ]
df_mde <- as.data.frame(mde_result)

# ============================================================
# 1A: Global MDE colored by major lineage, subsampled per lineage
# ============================================================
set.seed(1)
df_mde_a <- df_mde %>% tibble::rownames_to_column("cellID")
plot_df <- df_mde_a %>%
  inner_join(seurat_meta_filtered %>% select(cellID, annotation_level1), by = "cellID") %>%
  filter(annotation_level1 != "thymocyte")
max_total <- 1000000
min_per_group <- 300
cap_per_group <- 20000
group_sizes <- plot_df %>% count(annotation_level1, name = "n")
G <- nrow(group_sizes)
base_per_group <- ceiling(max_total / max(G, 1))
sample_plan <- group_sizes %>% mutate(n_take = pmin(n, pmax(min_per_group, pmin(cap_per_group, base_per_group))))
plot_df_sub <- plot_df %>%
  group_by(annotation_level1) %>%
  group_modify(~ dplyr::slice_sample(.x, n = sample_plan$n_take[sample_plan$annotation_level1 == .y$annotation_level1])) %>%
  ungroup()

p_1A <- ggplot(plot_df_sub, aes(x = MDE_1, y = MDE_2)) +
  scattermore::geom_scattermore(aes(color = annotation_level1), pointsize = 1.2) +
  scale_color_manual(values = ZemmourLib::immgent_colors$level1) +
  coord_equal() +
  theme_classic() +
  labs(title = "MDE: Annotation Level 1", x = "MDE 1", y = "MDE 2", color = "Cell Type") +
  theme(legend.text = element_text(size = 10), legend.key.size = unit(1.5, "lines")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggsave(filename = paste0(figure_path, "1A.pdf"), plot = p_1A, width = 5, height = 5)

# 1B: hand-finished Illustrator schematic -- not code-generated, no output here.

# Drop thymocytes from all downstream cell-level visualizations (matches
# Figure_Overview.R).
non_thymo_cells <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 != "thymocyte"]
L_pm_filtered <- L_pm_filtered[non_thymo_cells, ]
seurat_meta_filtered <- seurat_meta_filtered[non_thymo_cells, ]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[non_thymo_cells]

# ============================================================
# 1C: gene-program network -- the 200 GPs linked by their shared top signature
# genes, with a selected set of GPs highlighted. Each highlighted GP's color runs
# along its edges to its top signature genes (which are labeled); non-highlighted
# GPs are grey. Formal version: no legend / no GP-index labels -- the color->GP
# mapping and interpretation are in the Fig 1C caption (analysis/Figure1.Rmd).
# ============================================================
suppressPackageStartupMessages({
  library(igraph); library(tidygraph); library(ggraph)
})
N_TOP   <- 5                                      # top up-genes taken per GP
SIG_THR <- 0.1                                    # min per-GP-normalized loading
GP_HIGHLIGHTS <- c(                               # GP -> highlight color
  GP68  = "pink2",  GP58  = "orange2", GP35  = "purple", GP171 = "blue",
  GP1   = "cyan2",  GP56  = "red2",    GP161 = "brown",  GP6   = "green2",
  GP7   = "green3", GP196 = "yellow3")
GP_COL   <- "darkgrey"                            # non-highlighted GP nodes
GENE_COL <- "#C7A76C"                             # gene nodes (tan)

Fn <- F_pm_filtered
colnames(Fn) <- paste0("GP", seq_len(ncol(Fn)))
Fn <- sweep(Fn, 2, apply(abs(Fn), 2, max), "/")   # per-GP (column) max-abs norm
GPs <- colnames(Fn)

# each GP -> its top up-regulated signature genes (normalized loading >= SIG_THR)
top_up <- lapply(GPs, function(j) { x <- Fn[, j]
  u <- names(sort(x, decreasing = TRUE))[1:N_TOP]; u[x[u] >= SIG_THR] })
names(top_up) <- GPs
edges <- do.call(rbind, lapply(GPs, function(g)
  if (length(top_up[[g]])) data.frame(GP = g, Gene = top_up[[g]]) else NULL))

# bipartite GP<->gene graph (every GP + all its top genes); plain FR layout
gi <- graph_from_data_frame(edges, directed = FALSE)
g  <- as_tbl_graph(gi, directed = FALSE)
gp_set <- unique(edges$GP)
set.seed(1)
lay <- layout_with_fr(gi)
colnames(lay) <- c("x", "y"); rownames(lay) <- V(gi)$name
nm <- g %>% activate(nodes) %>% pull(name)

# highlight selected GPs: color their node + edges, label the genes they connect to
hl <- intersect(names(GP_HIGHLIGHTS), gp_set)
hl_genes <- setdiff(unique(unlist(
  lapply(hl, function(gp) neighbors(gi, gp)$name))), gp_set)
g <- g %>% activate(nodes) %>% mutate(
  is_gp = name %in% gp_set,
  gp_fill = ifelse(is_gp & name %in% hl, unname(GP_HIGHLIGHTS[name]), GP_COL),
  label_gene = !is_gp & name %in% hl_genes,
  gene_lab = ifelse(label_gene, name, ""),
  gp_size = ifelse(is_gp, 3, NA_real_))
g <- g %>% activate(edges) %>% mutate(
  gp_end = ifelse(.N()$is_gp[from], .N()$name[from], .N()$name[to]),
  gp_edge_highlight = gp_end %in% hl,
  gp_edge_col = ifelse(gp_edge_highlight, unname(GP_HIGHLIGHTS[gp_end]), NA_character_))

p_1C <- ggraph(g, layout = "manual", x = lay[nm, "x"], y = lay[nm, "y"]) +
  geom_edge_link(aes(filter = !gp_edge_highlight), color = "black", alpha = 0.18, width = 0.32) +
  geom_edge_link(aes(filter = gp_edge_highlight, edge_colour = gp_edge_col), alpha = 0.95, width = 1.25) +
  geom_node_point(aes(filter = !is_gp), shape = 16, size = 1, color = GENE_COL, alpha = 0.75) +
  geom_node_point(aes(filter = is_gp, size = gp_size, fill = gp_fill),
                  shape = 21, color = "white", stroke = 0.5) +
  geom_node_text(aes(filter = label_gene, label = gene_lab), repel = TRUE,
                 color = "black", size = 5, fontface = "italic", max.overlaps = Inf) +
  scale_fill_identity() + scale_edge_colour_identity() +
  scale_size_identity() + scale_edge_width_identity() +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  scale_y_continuous(expand = expansion(mult = 0.08)) +
  theme_void(base_size = 12) +
  theme(plot.margin = margin(10, 10, 10, 10), legend.position = "none")
ggsave(filename = paste0(figure_path, "1C.pdf"), plot = p_1C, width = 20, height = 20)

# ============================================================
# 1D: giant loading heatmap (200 GP loadings x a stratified cell sample; rows =
# cells by lineage x organ, columns = GPs clustered). The same GPs highlighted
# in the 1C network are marked here by a top color bar, a colored/bold column
# label, and a box around each GP column. See the Fig 1D caption (Figure1.Rmd).
# ============================================================
suppressPackageStartupMessages({ library(grid) })
GP_HIGHLIGHTS <- c(
  GP68 = "pink2",  GP58 = "orange2", GP35 = "purple", GP171 = "blue",
  GP1  = "cyan2",  GP56 = "red2",    GP161 = "brown", GP6  = "green2",
  GP7  = "green3", GP196 = "yellow3")

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
  rownames(L_pm_filtered)[order(x, decreasing = TRUE)[seq_len(K_ANCHOR)]]) |> as.vector() |> unique()
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
level1_present <- intersect(level1_order, as.character(unique(all_meta$annotation_level1)))
organ_present <- intersect(organ_order, as.character(unique(all_meta$organ_simplified)))
row_ann <- rowAnnotation(
  Cell_Type = factor(as.character(all_meta$annotation_level1), levels = level1_present),
  Organ = factor(as.character(all_meta$organ_simplified), levels = organ_present),
  col = list(Cell_Type = ZemmourLib::immgent_colors$level1[level1_present],
             Organ = ZemmourLib::immgent_colors$organ_simplified[organ_present]),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(Cell_Type = list(title = "Cell Type"), Organ = list(title = "Organ")))

gpn      <- colnames(L_display)
hl_val   <- ifelse(gpn %in% names(GP_HIGHLIGHTS), gpn, NA_character_)
lab_col  <- ifelse(gpn %in% names(GP_HIGHLIGHTS), GP_HIGHLIGHTS[gpn], "grey55")
lab_fs   <- ifelse(gpn %in% names(GP_HIGHLIGHTS), 9, 4)
lab_face <- ifelse(gpn %in% names(GP_HIGHLIGHTS), 2, 1)
# Publication version: only the highlighted GP columns keep an index label
# (background GP indices dropped), and the highlighted-GP legend is hidden -- the
# colour->GP mapping is given in the Fig 1D caption. The full-label, legended
# version for collaborators lives in experiments/gene_correlation_network/
# (heatmap_loading.*).
top_ann <- HeatmapAnnotation(
  `Highlighted GP` = hl_val, col = list(`Highlighted GP` = GP_HIGHLIGHTS),
  na_col = "white", simple_anno_size = unit(4, "mm"), annotation_name_gp = gpar(fontsize = 8),
  show_legend = FALSE)
ht_1D <- Heatmap(
  L_display, name = "Loading", col = col_fun, left_annotation = row_ann, top_annotation = top_ann,
  cluster_rows = FALSE, cluster_columns = TRUE,
  clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  column_labels = ifelse(gpn %in% names(GP_HIGHLIGHTS), gpn, ""),
  column_names_gp = gpar(col = lab_col, fontsize = lab_fs, fontface = lab_face),
  column_title = "Gene Programs (GPs)", column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  use_raster = TRUE, raster_quality = 3, border = FALSE,
  heatmap_legend_param = list(title = "Loading", direction = "vertical"))
pdf(paste0(figure_path, "1D.pdf"), width = 15, height = 20, useDingbats = FALSE)
ht_drawn <- draw(ht_1D, merge_legend = TRUE)
co <- column_order(ht_drawn); disp <- colnames(L_display)[co]; n_col_ht <- length(co)
decorate_heatmap_body("Loading", {
  for (gp in names(GP_HIGHLIGHTS)) { j <- which(disp == gp)
    if (length(j)) grid.rect(x = (j - 0.5) / n_col_ht, y = 0.5, width = 1.4 / n_col_ht, height = 1,
                             gp = gpar(col = GP_HIGHLIGHTS[gp], fill = NA, lwd = 2.5)) } })
dev.off()
