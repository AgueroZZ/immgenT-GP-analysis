## Signature-gene network via kNN BACKBONE (instead of a global corr threshold).
## Rationale: a global "edge if corr >= 0.5" graph is a dense hairball whose
## communities overlap too much for a force layout to separate, and it glues
## small specific programs (Treg) onto the broadly-correlated housekeeping mass.
## A kNN backbone keeps, per gene, only its top-k strongest correlates -- sparse,
## adaptive, and the standard basis for single-cell graph layouts. Small programs
## keep their own edges and separate cleanly.
##
## Pipeline: signature genes -> per-GP max-abs norm -> gene-gene Pearson corr ->
## per gene keep top-k neighbors with corr >= floor -> symmetrize -> Louvain ->
## FR layout -> color by dominant GP, annotate Treg landmarks.
## Run from repo root:
##   Rscript experiments/gene_correlation_network/knn_network.R [k] [floor]

suppressPackageStartupMessages({
  library(igraph); library(ggplot2); library(ggrastr)
})
data_path <- "data/"; out_dir <- "experiments/gene_correlation_network/"
args <- commandArgs(trailingOnly = TRUE)
k     <- if (length(args) >= 1) as.integer(args[1]) else 10
floor <- if (length(args) >= 2) as.numeric(args[2]) else 0.3
sig_thr <- 0.1; set.seed(1)
cat("kNN backbone: k =", k, "| corr floor =", floor, "\n")

F_pm <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F_pm) <- paste0("GP", seq_len(ncol(F_pm)))
F_pm <- sweep(F_pm, 2, apply(abs(F_pm), 2, max), "/")
sig  <- rownames(F_pm)[apply(abs(F_pm), 1, max) > sig_thr]
F_sig <- F_pm[sig, ]; n_gp <- ncol(F_sig)
dom_gp <- max.col(F_sig, ties.method = "first"); names(dom_gp) <- sig

rm_ <- rowMeans(F_sig); rs_ <- apply(F_sig, 1, sd)
keep <- rs_ > 0; Z <- (F_sig[keep, ] - rm_[keep]) / rs_[keep]
genes <- rownames(Z); G <- nrow(Z)
C <- tcrossprod(Z) / (n_gp - 1); diag(C) <- -Inf
cat("signature genes:", G, "\n")

# ── per-gene top-k neighbors (corr >= floor) ──────────────────────────────────
from <- integer(0); to <- integer(0)
for (i in seq_len(G)) {
  ord <- order(C[i, ], decreasing = TRUE)[1:k]
  ord <- ord[C[i, ord] >= floor]
  if (length(ord)) { from <- c(from, rep(i, length(ord))); to <- c(to, ord) }
}
g <- simplify(graph_from_data_frame(
  data.frame(from = genes[from], to = genes[to]), directed = FALSE,
  vertices = data.frame(name = genes, dom_gp = dom_gp[genes])))
cat("kNN graph:", vcount(g), "nodes,", ecount(g), "edges | components:",
    components(g)$no, "| largest:", max(components(g)$csize), "\n")

# ── Louvain communities ───────────────────────────────────────────────────────
cl <- cluster_louvain(g); memb <- membership(cl)
cat("Louvain:", length(cl), "communities | modularity:", round(modularity(cl), 3), "\n")

# Foxp3 community check
treg <- c("Foxp3","Il2ra","Ctla4","Ikzf2","Ikzf4","Tnfrsf18","Il10","Nrp1",
          "Tnfrsf4","Lrrc32","Ikzf1","Gpr83","Izumo1r","Cish","Lta")
if ("Foxp3" %in% genes) {
  mem <- names(memb)[memb == memb["Foxp3"]]
  cat("Foxp3 community size:", length(mem),
      "| Treg markers:", paste(intersect(treg, mem), collapse = ", "), "\n")
  dtab <- sort(table(dom_gp[mem]), decreasing = TRUE)
  cat("  top GPs:", paste(names(head(dtab, 4)), "=", head(dtab, 4), collapse = ", "), "\n")
}

# ── FR layout ─────────────────────────────────────────────────────────────────
cat("FR layout...\n"); t0 <- Sys.time()
lay <- layout_with_fr(g, niter = 2000)
cat("done in", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")
rownames(lay) <- V(g)$name
df <- data.frame(x = lay[, 1], y = lay[, 2], gene = V(g)$name,
                 dom_gp = V(g)$dom_gp, comm = as.integer(memb[V(g)$name]))
rownames(df) <- df$gene
saveRDS(list(layout = df, graph = g, membership = memb, dom_gp = dom_gp),
        paste0(out_dir, "knn_layout.rds"))

# ── edges for rendering ───────────────────────────────────────────────────────
el <- as_edgelist(g, names = TRUE)
edf <- data.frame(x = df[el[, 1], "x"], y = df[el[, 1], "y"],
                  xend = df[el[, 2], "x"], yend = df[el[, 2], "y"])
treg_df <- df[intersect(treg, df$gene), ]

mk <- function(color_by, title) {
  aes_col <- if (color_by == "gp") aes(x, y, color = dom_gp) else
                                   aes(x, y, color = factor(comm))
  p <- ggplot() +
    rasterise(geom_segment(data = edf, aes(x, y, xend = xend, yend = yend),
                color = "grey75", linewidth = 0.1, alpha = 0.25), dpi = 150) +
    rasterise(geom_point(data = df, aes_col, size = 0.6, alpha = 0.9), dpi = 150) +
    geom_point(data = treg_df, aes(x, y), shape = 21, size = 2.4,
               color = "black", fill = NA, stroke = 0.7) +
    ggrepel::geom_text_repel(data = treg_df, aes(x, y, label = gene),
               size = 3, fontface = "bold", max.overlaps = 30,
               min.segment.length = 0, segment.size = 0.3) +
    coord_equal() + theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank()) +
    labs(title = title, x = NULL, y = NULL,
         subtitle = paste0(vcount(g), " signature genes | kNN(k=", k,
                           ") backbone | FR layout | Treg landmarks labeled"))
  if (color_by == "gp")
    p + scale_color_gradientn(colors = grDevices::rainbow(200, end = 0.9),
          limits = c(1, 200), breaks = c(1, 50, 100, 150, 200), name = "dominant GP") +
        theme(legend.key.height = unit(1.1, "cm"))
  else p + scale_color_manual(values = rep(grDevices::rainbow(20), 20), guide = "none")
}

ggsave(paste0(out_dir, "knn_network_domGP.png"), mk("gp",
       "Signature-gene kNN correlation network -- colored by dominant GP"),
       width = 9, height = 7.5, dpi = 150)
ggsave(paste0(out_dir, "knn_network_community.png"), mk("comm",
       "Signature-gene kNN correlation network -- colored by Louvain community"),
       width = 9, height = 7.5, dpi = 150)
ggsave(paste0(out_dir, "knn_network_domGP.pdf"), mk("gp",
       "Signature-gene kNN correlation network -- colored by dominant GP"),
       width = 9, height = 7.5)
cat("saved knn_network_{domGP,community}.png + pdf\n")
