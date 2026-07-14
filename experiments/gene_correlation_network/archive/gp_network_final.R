## Broad-picture GP-GP network (polished deliverable).
## 200 GPs as nodes, linked by specificity-weighted shared top genes; Louvain
## families as colored convex-hull regions; node size ~ PVE. This is the
## 200-GP "broad picture" analog of the Fig-3d bipartite network (genes are
## folded into the edges rather than drawn as separate nodes).
##
## Edge rule (see gp_network_specific.R): each GP -> top 5 up + 5 down genes
## (per-GP max-abs norm, |loading|>=0.1); GP-GP edge weight = sum 1/df(shared
## gene); genes shared by > df_max GPs are dropped as non-specific glue.
##
## Run from repo root:
##   Rscript experiments/gene_correlation_network/gp_network_final.R [df_max] [min_w]

suppressPackageStartupMessages({ library(igraph); library(ggplot2); library(ggrepel) })
data_path <- "data/"; out_dir <- "experiments/gene_correlation_network/"
args <- commandArgs(trailingOnly = TRUE)
df_max <- if (length(args) >= 1) as.integer(args[1]) else 10
min_w  <- if (length(args) >= 2) as.numeric(args[2]) else 0.3
set.seed(1)

F_pm <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F_pm) <- paste0("GP", seq_len(ncol(F_pm)))
F_pm <- sweep(F_pm, 2, apply(abs(F_pm), 2, max), "/")
nGP <- ncol(F_pm)

sets <- lapply(seq_len(nGP), function(j) {
  v  <- F_pm[, j]
  up <- names(sort(v, decreasing = TRUE))[1:5]; up <- up[v[up] >=  0.1]
  dn <- names(sort(v))[1:5];                    dn <- dn[v[dn] <= -0.1]
  unique(c(up, dn))
})
names(sets) <- paste0("GP", seq_len(nGP))
top_gene <- sapply(seq_len(nGP), function(j) names(which.max(F_pm[, j])))

# PVE for node size
pve <- setNames(rep(NA_real_, nGP), paste0("GP", seq_len(nGP)))
pve_file <- "data/factors_comments_updated.xlsx - factor_pve.csv"
if (file.exists(pve_file)) {
  pt <- read.csv(pve_file, check.names = FALSE)
  pve[paste0("GP", as.integer(gsub("F", "", pt$factor)))] <- pt$pve
}

df_gene <- table(unlist(sets))
drop <- names(df_gene)[df_gene > df_max]
w_gene <- 1 / as.numeric(df_gene); names(w_gene) <- names(df_gene)

ed <- data.frame(from = character(), to = character(), w = numeric())
for (i in 1:(nGP - 1)) for (j in (i + 1):nGP) {
  sh <- setdiff(intersect(sets[[i]], sets[[j]]), drop)
  if (length(sh)) { w <- sum(w_gene[sh])
    if (w >= min_w) ed <- rbind(ed, data.frame(from = paste0("GP", i), to = paste0("GP", j), w = w)) }
}
g <- graph_from_data_frame(ed, directed = FALSE,
       vertices = data.frame(name = paste0("GP", seq_len(nGP)), top = top_gene, pve = pve))
g <- delete_vertices(g, degree(g) == 0)
memb <- membership(cluster_louvain(g, weights = E(g)$w))
cat("GPs:", vcount(g), "| edges:", ecount(g), "| families:", length(unique(memb)),
    "| modularity:", round(modularity(cluster_louvain(g, weights = E(g)$w)), 3), "\n")

lay <- layout_with_fr(g, weights = E(g)$w, niter = 3000)
df <- data.frame(x = lay[, 1], y = lay[, 2], gp = V(g)$name,
                 comm = factor(memb[V(g)$name]), top = V(g)$top, pve = V(g)$pve)
df$size <- if (all(is.na(df$pve))) 3.2 else 2.4 + 6 * sqrt(pmax(df$pve, 0) / max(df$pve, na.rm = TRUE))

el <- as_edgelist(g, names = TRUE)
edf <- data.frame(x = df$x[match(el[, 1], df$gp)], y = df$y[match(el[, 1], df$gp)],
                  xend = df$x[match(el[, 2], df$gp)], yend = df$y[match(el[, 2], df$gp)], w = E(g)$w)

# convex hull per family (>=3 GPs) for the shaded region
hull <- do.call(rbind, lapply(split(df, df$comm), function(d)
  if (nrow(d) >= 3) d[chull(d$x, d$y), ] else NULL))

ncomm <- length(levels(df$comm))
pal <- setNames(rep(grDevices::rainbow(12, end = 0.88), length.out = ncomm), levels(df$comm))

p <- ggplot() +
  geom_polygon(data = hull, aes(x, y, group = comm, fill = comm),
               alpha = 0.12, color = NA) +
  geom_segment(data = edf, aes(x, y, xend = xend, yend = yend, linewidth = w),
               color = "grey65", alpha = 0.45) +
  geom_point(data = df, aes(x, y, color = comm, size = size)) +
  geom_text_repel(data = df, aes(x, y, label = paste0(gp, "\n", top), color = comm),
                  size = 2.3, lineheight = 0.8, max.overlaps = 60,
                  segment.size = 0.15, show.legend = FALSE) +
  scale_color_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_linewidth(range = c(0.2, 1.8), guide = "none") +
  scale_size_identity() +
  coord_equal() + theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.margin = margin(12, 12, 12, 12)) +
  labs(title = "Gene-program network -- broad picture of 200 GPs",
       subtitle = paste0(vcount(g), " GPs | edges = specificity-weighted shared top genes | ",
         length(unique(memb)), " Louvain families (shaded) | node size ~ PVE"))

ggsave(paste0(out_dir, "gp_network_final.png"), p, width = 12, height = 10, dpi = 150, bg = "white")
ggsave(paste0(out_dir, "gp_network_final.pdf"), p, width = 12, height = 10, bg = "white")
cat("saved gp_network_final.{png,pdf}\n")
