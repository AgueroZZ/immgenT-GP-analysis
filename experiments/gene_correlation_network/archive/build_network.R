## Gene-gene correlation network from the GP factor matrix.
##
## Motivation: if genes naturally organize into programs, then correlating each
## gene's GP loading profile and drawing edges between highly-correlated genes
## should reveal coherent clusters -- and those clusters should line up with the
## GPs themselves. This is a graph-based counterpart to the gene t-SNE atlas.
##
## Pipeline (per user spec):
##   1. F_pm_filtered.rds, per-GP (per-column) max-abs normalization
##   2. each gene = its 200-dim GP loading vector; Pearson corr between genes
##   3. edge between genes with corr >= 0.5  (=> ~1.37M edges, near-complete
##      hairball: 99.7% of genes connected, median degree 108)
##   4. big-picture layout: DrL force-directed on the FULL corr>=0.5 graph to
##      declutter clusters; nodes colored by dominant GP; edges NOT drawn
##      (1.37M edges render as a black smear -- the layout already encodes them).
##
## Run from repo root:
##   Rscript experiments/gene_correlation_network/build_network.R

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(ggrastr)
})

data_path <- "data/"
out_dir   <- "experiments/gene_correlation_network/"
corr_thr  <- 0.5
set.seed(1)

# ── standardized gene matrix (reuse exploration output if present) ─────────────
expl <- paste0(out_dir, "explore_scale.rds")
F_pm <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F_pm) <- paste0("GP", seq_len(ncol(F_pm)))
F_pm <- sweep(F_pm, 2, apply(abs(F_pm), 2, max), "/")     # per-GP max-abs
n_gp <- ncol(F_pm)

if (file.exists(expl)) {
  Z <- readRDS(expl)$Z
} else {
  rmean <- rowMeans(F_pm); rsd <- apply(F_pm, 1, sd)
  Z <- (F_pm[rsd > 0, ] - rmean[rsd > 0]) / rsd[rsd > 0]
}
G     <- nrow(Z)
genes <- rownames(Z)
cat("standardized matrix:", G, "genes x", n_gp, "GPs\n")

# ── dominant GP per gene (argmax normalized loading, original GP index) ─────────
dom_gp_num <- max.col(F_pm[genes, , drop = FALSE], ties.method = "first")

# ── build edge list at corr >= thr (chunked upper triangle) ────────────────────
cat("extracting edges at corr >=", corr_thr, "...\n")
chunk <- 1000
from_l <- integer(0); to_l <- integer(0)
frag_from <- vector("list"); frag_to <- vector("list"); k <- 0
for (start in seq(1, G, by = chunk)) {
  idx <- start:min(start + chunk - 1, G)
  C   <- tcrossprod(Z[idx, , drop = FALSE], Z) / (n_gp - 1)   # |idx| x G
  # mask to strict upper triangle (global j > global i)
  for (r in seq_along(idx)) C[r, 1:idx[r]] <- NA
  hit <- which(C >= corr_thr, arr.ind = TRUE)
  if (nrow(hit)) {
    k <- k + 1
    frag_from[[k]] <- idx[hit[, 1]]
    frag_to[[k]]   <- hit[, 2]
  }
}
from_l <- unlist(frag_from); to_l <- unlist(frag_to)
cat("edges:", length(from_l), "\n")

g <- graph_from_data_frame(
  data.frame(from = genes[from_l], to = genes[to_l]),
  directed = FALSE,
  vertices = data.frame(name = genes, dom_gp = dom_gp_num))
cat("graph:", vcount(g), "nodes,", ecount(g), "edges\n")
cat("components:", components(g)$no,
    "| largest:", max(components(g)$csize), "\n")

# ── layout: DrL on full graph (declutters dense hairballs) ─────────────────────
cat("running DrL layout (this can take a few minutes)...\n")
t0 <- Sys.time()
lay <- layout_with_drl(g, options = drl_defaults$default)
cat("layout done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

df <- data.frame(x = lay[, 1], y = lay[, 2],
                 gene = V(g)$name, dom_gp = V(g)$dom_gp)
saveRDS(list(layout = df, graph_summary =
               list(n = vcount(g), e = ecount(g))),
        paste0(out_dir, "network_layout.rds"))

# ── plot: nodes colored by dominant GP, edges NOT drawn, rasterized ────────────
p <- ggplot(df, aes(x, y, color = dom_gp)) +
  rasterise(geom_point(size = 0.35, alpha = 0.8), dpi = 150) +
  scale_color_gradientn(colors = grDevices::rainbow(200, end = 0.9),
                        limits = c(1, 200), breaks = c(1, 50, 100, 150, 200),
                        name = "dominant GP") +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), legend.key.height = unit(1.1, "cm")) +
  labs(title = "Gene-gene correlation network (corr >= 0.5 on GP loadings)",
       subtitle = paste0(vcount(g), " genes, ", format(ecount(g), big.mark = ","),
                         " edges | DrL layout | nodes colored by dominant GP | edges not drawn"),
       x = NULL, y = NULL)
ggsave(paste0(out_dir, "gene_network_domGP.pdf"), p, width = 8.5, height = 7)
ggsave(paste0(out_dir, "gene_network_domGP.png"), p, width = 8.5, height = 7, dpi = 150)
cat("saved gene_network_domGP.{pdf,png}\n")
