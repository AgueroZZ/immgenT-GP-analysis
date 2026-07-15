## GP-highlight GP-gene network.
##
## Highlight a chosen set of GPs on the GP-gene bipartite network: each GP's color
## flows to its node and its edges, and the genes each highlighted GP connects to
## are labeled. Produces two versions:
##   * gp_highlight_formal.{png,pdf}   -- no legend, no GP-index labels (for figures)
##   * gp_highlight_internal.{png,pdf} -- legend + GP-index labels (for the team)
##
## To highlight a DIFFERENT set of GPs, edit GP_HIGHLIGHTS below and rerun -- the
## graph, layout, colors and everything else stay the same.
##
## Run from repo root:
##   Rscript experiments/gene_correlation_network/gp_highlight_selected.R

suppressPackageStartupMessages({
  library(tidygraph); library(ggraph); library(igraph); library(ggplot2)
})

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG  -- edit GP_HIGHLIGHTS to change which GPs are highlighted (that's all)
# ══════════════════════════════════════════════════════════════════════════════
GP_HIGHLIGHTS <- c(                       # GP -> highlight color
  GP68  = "pink2",  GP58  = "orange2", GP35  = "purple", GP171 = "blue",
  GP1   = "cyan2",  GP56  = "red2",    GP161 = "brown",  GP6   = "green2",
  GP7   = "green3", GP196 = "yellow3"
)
N_TOP    <- 5            # top up-regulated genes per GP (|loading| >= SIG_THR)
SIG_THR  <- 0.1
SEED     <- 1            # layout seed -- fixes the arrangement
GP_COL   <- "darkgrey"   # non-highlighted GP nodes
GENE_COL <- "#C7A76C"    # gene nodes (tan)
CANVAS_W <- 20; CANVAS_H <- 20
data_path <- "data/"
out_dir   <- "experiments/gene_correlation_network/"

# ══════════════════════════════════════════════════════════════════════════════
# 1. build the bipartite GP<->gene graph + fixed FR layout
# ══════════════════════════════════════════════════════════════════════════════
build_gp_gene_graph <- function(data_path, n_top, sig_thr, seed) {
  F <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
  colnames(F) <- paste0("GP", seq_len(ncol(F)))
  F <- sweep(F, 2, apply(abs(F), 2, max), "/")          # per-GP max-abs norm
  GPs <- colnames(F)
  top_up <- lapply(GPs, function(j) { x <- F[, j]
    u <- names(sort(x, decreasing = TRUE))[1:n_top]; u[x[u] >= sig_thr] })
  names(top_up) <- GPs
  edges <- do.call(rbind, lapply(GPs, function(gp)
    if (length(top_up[[gp]])) data.frame(GP = gp, Gene = top_up[[gp]]) else NULL))
  gi <- graph_from_data_frame(edges, directed = FALSE)
  set.seed(seed)
  lay <- layout_with_fr(gi)                              # plain Fruchterman-Reingold
  colnames(lay) <- c("x", "y"); rownames(lay) <- V(gi)$name
  g <- as_tbl_graph(gi, directed = FALSE)
  list(g = g, gi = gi, gp_set = unique(edges$GP), lay = lay,
       nm = g %>% activate(nodes) %>% pull(name))
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. annotate nodes/edges for a given highlight set
# ══════════════════════════════════════════════════════════════════════════════
add_gp_highlights <- function(gr, gp_highlights, gp_col, gene_col) {
  hl <- intersect(names(gp_highlights), gr$gp_set)
  miss <- setdiff(names(gp_highlights), gr$gp_set)
  if (length(miss)) warning("Highlighted GPs absent from graph: ", paste(miss, collapse = ", "))
  hl_genes <- setdiff(unique(unlist(
    lapply(hl, function(gp) neighbors(gr$gi, gp)$name))), gr$gp_set)
  g <- gr$g %>% activate(nodes) %>% mutate(
    is_gp = name %in% gr$gp_set,
    gp_lab = ifelse(is_gp, sub("^GP", "", name), ""),
    gp_highlight = is_gp & name %in% hl,
    gp_fill = ifelse(gp_highlight, unname(gp_highlights[name]), gp_col),
    label_gene = !is_gp & name %in% hl_genes,
    gene_lab = ifelse(label_gene, name, ""),
    gp_size = ifelse(is_gp, 3, NA_real_))
  g <- g %>% activate(edges) %>% mutate(
    gp_end = ifelse(.N()$is_gp[from], .N()$name[from], .N()$name[to]),
    gp_edge_highlight = gp_end %in% hl,
    gp_edge_col = ifelse(gp_edge_highlight, unname(gp_highlights[gp_end]), NA_character_))
  list(g = g, highlighted_gps = hl)
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. plot; toggle legend + GP-index labels
# ══════════════════════════════════════════════════════════════════════════════
make_gp_plot <- function(g, lay, nm, gp_highlights, highlighted_gps, gene_col,
                         show_legend = TRUE, show_gp_labels = TRUE) {
  p <- ggraph(g, layout = "manual", x = lay[nm, "x"], y = lay[nm, "y"]) +
    geom_edge_link(aes(filter = !gp_edge_highlight), color = "black",
                   alpha = 0.18, width = 0.32) +
    geom_edge_link(aes(filter = gp_edge_highlight, edge_colour = gp_edge_col),
                   alpha = 0.95, width = 1.25) +
    geom_node_point(aes(filter = !is_gp), shape = 16, size = 1,
                    color = gene_col, alpha = 0.75) +
    geom_node_point(aes(filter = is_gp, size = gp_size, fill = gp_fill),
                    shape = 21, color = "white", stroke = 0.5) +
    geom_node_text(aes(filter = label_gene, label = gene_lab), repel = TRUE,
                   color = "black", size = 5, fontface = "italic", max.overlaps = Inf) +
    scale_fill_identity(name = "Highlighted GP",
                        breaks = unname(gp_highlights[highlighted_gps]),
                        labels = highlighted_gps,
                        guide = if (show_legend) "legend" else "none") +
    scale_edge_colour_identity() + scale_size_identity() + scale_edge_width_identity() +
    scale_x_continuous(expand = expansion(mult = 0.08)) +
    scale_y_continuous(expand = expansion(mult = 0.08)) +
    theme_void(base_size = 12) +
    theme(plot.margin = margin(10, 10, 10, 10),
          legend.position = if (show_legend) "right" else "none",
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 10)) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 5)))
  if (show_gp_labels)
    p <- p + geom_node_text(aes(filter = is_gp, label = gp_lab), color = "black", size = 1.9)
  p
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. run -- formal (no legend/labels) + internal (legend + GP indices)
# ══════════════════════════════════════════════════════════════════════════════
gr  <- build_gp_gene_graph(data_path, N_TOP, SIG_THR, SEED)
hgl <- add_gp_highlights(gr, GP_HIGHLIGHTS, GP_COL, GENE_COL)

save_plot <- function(p, name)
  for (ext in c("png", "pdf"))
    ggsave(sprintf("%s%s.%s", out_dir, name, ext), p,
           width = CANVAS_W, height = CANVAS_H,
           dpi = 150, bg = "white", limitsize = FALSE)

save_plot(make_gp_plot(hgl$g, gr$lay, gr$nm, GP_HIGHLIGHTS, hgl$highlighted_gps,
                       GENE_COL, show_legend = FALSE, show_gp_labels = FALSE),
          "gp_highlight_formal")
save_plot(make_gp_plot(hgl$g, gr$lay, gr$nm, GP_HIGHLIGHTS, hgl$highlighted_gps,
                       GENE_COL, show_legend = TRUE,  show_gp_labels = TRUE),
          "gp_highlight_internal")

cat("highlighted GPs:", paste(hgl$highlighted_gps, collapse = ", "), "\n")
cat("saved gp_highlight_formal.{png,pdf} and gp_highlight_internal.{png,pdf}\n")
