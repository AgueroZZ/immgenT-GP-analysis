# Figure 1. Overview.
#
# Panels produced (final files figures/Previous/bits/Figure 1/1A.pdf .. 1I.pdf;
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
# an exact byte match between figures/Previous/bits/Figure 1/1A.pdf
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

# --- doc:setup ---
library(ggplot2)
library(dplyr)
library(scattermore)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(Matrix) # protein_mat_normalized_lognorm is a dgCMatrix; must be
                # attached (not just loaded) for `[` subsetting to dispatch

data_path <- "data/"
figure_path <- "figures/final-selected/Figure 1/"
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

# --- doc:setup2 ---
# Drop thymocytes from all downstream cell-level visualizations.
non_thymo_cells <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 != "thymocyte"]
L_pm_filtered <- L_pm_filtered[non_thymo_cells, ]
seurat_meta_filtered <- seurat_meta_filtered[non_thymo_cells, ]
protein_mat_normalized_lognorm <- protein_mat_normalized_lognorm[non_thymo_cells]

# ============================================================
# 1C: gene-program network -- the 200 GPs linked by their shared top signature
# genes, with a selected set of GPs highlighted. Each highlighted GP's color runs
# along its edges to its top signature genes (which are labeled); non-highlighted
# GPs are grey. No legend and no GP-index labels.
# --- internal ---
# The colour -> GP mapping used to be spelled out in the Fig. 1c caption on
# analysis/Figure1.Rmd. That page now carries the published caption verbatim,
# which does not list it; GP_HIGHLIGHTS in the 1D block below is the only
# record of which GP takes which colour.
# --- end internal ---
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
# 1D: giant loading heatmap. Rows are the 200 GPs (clustered); columns are a
# stratified cell sample blocked lineage > annotation_level2_group > level2.
# The GPs highlighted in the 1C network are marked here by a left colour bar,
# a coloured/bold row label, and a box around each highlighted GP row.
#
# Two things to know about the cell set:
#   * Cells whose annotation_level2_group is "miniverse" are excluded. Those are
#     exactly the seven ".wM" level2 clusters (CD8.wM, CD4.wM, Treg.wM, gdT.wM,
#     CD8aa.wM, Tz.wM, DN.wM), so 95 level2 clusters become 88. The exclusion is
#     applied BEFORE sampling and before the anchor-cell search, so every GP
#     still contributes K_ANCHOR top-loading cells that survive into the plot.
#   * Cells are stratified by level2: every cluster with at least MIN_CELLS
#     cells is sampled at N_SAMPLE.
# --- internal ---
# Rebuilt 2026-09-09. The previous version of this panel put the GPs on columns
# and stratified cells by organ_simplified (top 5 organs per lineage); it had no
# level2_group bar and did not exclude miniverse. The internal experiment behind
# the change -- colour-scale variants, and a reconciliation of every level2
# cluster in the atlas against what this panel draws -- is in gitignored
# experiments/fig1d_gp_rows_by_level2/.
#
# The Fig. 1d caption on analysis/Figure1.Rmd is an interim stand-in: the
# published legend says "from different organs and lineages", which this panel
# no longer shows, and Ziang is sending the replacement wording separately.
# --- end internal ---
#
# Column order: within each lineage the groups are sorted by the first
# (alphabetically lowest) level2 cluster they contain, and level2 orders cells
# inside each group. level2 is therefore alphabetical within a group but not
# across a whole lineage -- each lineage's ".P" (proliferating) cluster is the
# single exception, landing at the end of its lineage rather than mid-alphabet.
# LEVEL2_GROUP_ORDER below fixes the legend order only, not the block order.
# ============================================================
suppressPackageStartupMessages({ library(grid) })
GP_HIGHLIGHTS <- c(
  GP68 = "pink2",  GP58 = "orange2", GP35 = "purple", GP171 = "blue",
  GP1  = "cyan2",  GP56 = "red2",    GP161 = "brown", GP6  = "green2",
  GP7  = "green3", GP196 = "yellow3")

# No annotation_level2_group order or palette exists in ZemmourLib, so both are
# defined here. The colours are a deliberately different family from the level1
# primaries so the two annotation bars cannot be confused, and none is
# near-white (against the white heatmap body a pale category would read as
# missing data rather than as a level).
LEVEL2_GROUP_ORDER  <- c("resting", "activated", "proliferating", "miniverse", "other")
LEVEL2_GROUP_COLORS <- c(resting       = "#80cdc1", activated = "#b2182b",
                         proliferating = "#542788", miniverse = "#8c510a",
                         other         = "#bdbdbd")
EXCLUDE_LEVEL2_GROUPS <- c("miniverse")

set.seed(6173)
MIN_CELLS <- 20; N_SAMPLE <- 80; K_ANCHOR <- 5
level1_order <- c("CD8", "CD4", "Treg", "gdT", "CD8aa", "Tz", "DN")

meta_1d <- seurat_meta_filtered |>
  dplyr::filter(annotation_level1 %in% level1_order,
                !is.na(annotation_level2),
                !annotation_level2_group %in% EXCLUDE_LEVEL2_GROUPS)
stopifnot(nrow(meta_1d) > 0)
L_1d <- L_pm_filtered[meta_1d$cellID, ]

# level2 blocks follow the level1 order, alphabetised within each lineage.
level2_order <- meta_1d |>
  dplyr::distinct(annotation_level1, annotation_level2) |>
  dplyr::arrange(factor(annotation_level1, levels = level1_order),
                 as.character(annotation_level2)) |>
  dplyr::pull(annotation_level2) |> as.character()

sampled_random <- meta_1d |>
  dplyr::group_by(annotation_level2) |>
  dplyr::filter(dplyr::n() >= MIN_CELLS) |>
  dplyr::slice_sample(n = N_SAMPLE) |> dplyr::ungroup()
anchor_cellids <- apply(L_1d, 2, function(x)
  rownames(L_1d)[order(x, decreasing = TRUE)[seq_len(K_ANCHOR)]]) |>
  as.vector() |> unique()
anchor_meta <- meta_1d |> dplyr::filter(cellID %in% anchor_cellids)

all_meta <- dplyr::bind_rows(sampled_random, anchor_meta) |>
  dplyr::distinct(cellID, .keep_all = TRUE) |>
  dplyr::filter(as.character(annotation_level2) %in% level2_order)

# Sort keys: lineage, then group by its first level2 cluster, then level2.
l2g     <- as.character(all_meta$annotation_level2_group)
l1_rank <- match(as.character(all_meta$annotation_level1), level1_order)
l2_rank <- match(as.character(all_meta$annotation_level2), level2_order)
stopifnot(!anyNA(l1_rank), !anyNA(l2_rank),
          all(l2g %in% LEVEL2_GROUP_ORDER))
blk      <- paste(l1_rank, l2g, sep = "|")
grp_rank <- tapply(l2_rank, blk, min)[blk]
all_meta <- all_meta[order(l1_rank, grp_rank, l2_rank), ]
l2g <- as.character(all_meta$annotation_level2_group)

level1_present       <- intersect(level1_order, unique(as.character(all_meta$annotation_level1)))
level2_present       <- intersect(level2_order, unique(as.character(all_meta$annotation_level2)))
level2_group_present <- intersect(LEVEL2_GROUP_ORDER, unique(l2g))
# Each lineage x group pair must be one contiguous block on the axis.
stopifnot(!any(duplicated(rle(paste(l1_rank[order(l1_rank, grp_rank, l2_rank)], l2g))$values)))

L_sampled <- L_1d[all_meta$cellID, ]
clip_val <- quantile(L_sampled, 0.99)
L_display <- pmin(L_sampled, clip_val)
colnames(L_display) <- gsub("^K", "GP", colnames(L_display))
mat_1D <- t(L_display)   # rows = GPs, columns = cells
col_fun <- colorRamp2(c(0, clip_val / 2, clip_val), c("white", "#4393c3", "#08306b"))

col_ann <- HeatmapAnnotation(
  Cell_Type    = factor(as.character(all_meta$annotation_level1), levels = level1_present),
  Level2_Group = factor(l2g, levels = level2_group_present),
  col = list(Cell_Type    = ZemmourLib::immgent_colors$level1[level1_present],
             Level2_Group = LEVEL2_GROUP_COLORS[level2_group_present]),
  simple_anno_size = unit(4, "mm"), annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(Cell_Type    = list(title = "Cell Type"),
                                 Level2_Group = list(title = "Level2 Group")))

gpn      <- rownames(mat_1D)
hl_val   <- ifelse(gpn %in% names(GP_HIGHLIGHTS), gpn, NA_character_)
lab_col  <- ifelse(gpn %in% names(GP_HIGHLIGHTS), GP_HIGHLIGHTS[gpn], "grey55")
lab_fs   <- ifelse(gpn %in% names(GP_HIGHLIGHTS), 9, 4)
lab_face <- ifelse(gpn %in% names(GP_HIGHLIGHTS), 2, 1)
# Only the highlighted GP rows keep an index label (background GP indices
# dropped), and the highlighted-GP legend is hidden.
row_ann <- rowAnnotation(
  `Highlighted GP` = hl_val, col = list(`Highlighted GP` = GP_HIGHLIGHTS),
  na_col = "white", simple_anno_size = unit(4, "mm"),
  annotation_name_gp = gpar(fontsize = 8), show_legend = FALSE)
ht_1D <- Heatmap(
  mat_1D, name = "Loading", col = col_fun, left_annotation = row_ann, top_annotation = col_ann,
  cluster_rows = TRUE, cluster_columns = FALSE,
  clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
  show_column_names = FALSE,
  row_labels = ifelse(gpn %in% names(GP_HIGHLIGHTS), gpn, ""),
  row_names_gp = gpar(col = lab_col, fontsize = lab_fs, fontface = lab_face),
  row_title = "Gene Programs (GPs)", row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  use_raster = TRUE, raster_quality = 3, border = FALSE,
  heatmap_legend_param = list(title = "Loading", direction = "vertical"))
pdf(paste0(figure_path, "1D.pdf"), width = 20, height = 15, useDingbats = FALSE)
ht_drawn <- draw(ht_1D, merge_legend = TRUE)
ro <- row_order(ht_drawn); disp <- rownames(mat_1D)[ro]; n_row_ht <- length(ro)
decorate_heatmap_body("Loading", {
  for (gp in names(GP_HIGHLIGHTS)) { i <- which(disp == gp)
    if (length(i)) grid.rect(x = 0.5, y = 1 - (i - 0.5) / n_row_ht, width = 1, height = 1.4 / n_row_ht,
                             gp = gpar(col = GP_HIGHLIGHTS[gp], fill = NA, lwd = 2.5)) } })
dev.off()
cat(sprintf("1D: %d cells, %d level2 clusters, %d groups, clip_val=%.4f\n",
            nrow(all_meta), length(level2_present), length(level2_group_present), clip_val))
