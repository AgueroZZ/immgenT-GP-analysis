## GP-GP network v2: edges weighted by SHARED-GENE SPECIFICITY (TF-IDF style).
##
## v1 (gp_network.R) glued all immune GPs into one central blob because broadly-
## expressed genes (Tmsb4x/Tmsb10, Ccl5, Gzma/Gzmb, ribosomal) sit in dozens of
## GP top-sets and link everything. Here a shared gene contributes 1/df(gene),
## where df = # of GP top-sets containing it, and ultra-promiscuous genes
## (df > df_max) are dropped entirely. Rare shared genes (a marker in just 2-3
## GPs) dominate the weight, so specific families resolve.
##
## Run from repo root:
##   Rscript experiments/gene_correlation_network/gp_network_specific.R [df_max] [min_w]

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

# document frequency of each top-gene across GP sets
df_gene <- table(unlist(sets))
cat("most promiscuous shared genes (df):\n"); print(head(sort(df_gene, decreasing = TRUE), 12))
drop <- names(df_gene)[df_gene > df_max]
cat("dropping", length(drop), "genes with df >", df_max, "as glue\n")
w_gene <- 1 / as.numeric(df_gene); names(w_gene) <- names(df_gene)

# specificity-weighted GP-GP edges
ed <- data.frame(from = character(), to = character(), w = numeric())
for (i in 1:(nGP - 1)) for (j in (i + 1):nGP) {
  sh <- setdiff(intersect(sets[[i]], sets[[j]]), drop)
  if (!length(sh)) next
  w <- sum(w_gene[sh])
  if (w >= min_w) ed <- rbind(ed, data.frame(
    from = paste0("GP", i), to = paste0("GP", j), w = w))
}
cat("edges (specificity weight >=", min_w, "):", nrow(ed), "\n")

g <- graph_from_data_frame(ed, directed = FALSE,
       vertices = data.frame(name = paste0("GP", seq_len(nGP)), top = top_gene))
g <- delete_vertices(g, degree(g) == 0)
cat("graph:", vcount(g), "GP nodes,", ecount(g), "edges | components:",
    components(g)$no, "| largest:", max(components(g)$csize), "\n")

cl <- cluster_louvain(g, weights = E(g)$w); memb <- membership(cl)
cat("Louvain:", length(cl), "families | modularity:", round(modularity(cl), 3), "\n")

cat("\n== GP families (specificity-weighted) ==\n")
for (c in sort(unique(memb))) {
  gps <- names(memb)[memb == c]
  gtab <- sort(table(unlist(lapply(sets[gps], setdiff, drop))), decreasing = TRUE)
  cat(sprintf("family %d (%d GPs): %s\n   hallmark genes: %s\n", c, length(gps),
              paste(gps, collapse = ", "), paste(names(head(gtab, 12)), collapse = ", ")))
}
foxp3_gp <- paste0("GP", which.max(F_pm["Foxp3", ]))
if (foxp3_gp %in% names(memb)) {
  fam <- names(memb)[memb == memb[foxp3_gp]]
  cat("\nFoxp3 GP =", foxp3_gp, "| family:", paste(fam, collapse = ", "), "\n")
  cat("  neighbors:", paste(names(neighbors(g, foxp3_gp)), collapse = ", "), "\n")
}

lay <- layout_with_fr(g, weights = E(g)$w, niter = 3000)
df <- data.frame(x = lay[, 1], y = lay[, 2], gp = V(g)$name,
                 comm = factor(memb[V(g)$name]), top = V(g)$top)
el <- as_edgelist(g, names = TRUE)
edf <- data.frame(x = df$x[match(el[, 1], df$gp)], y = df$y[match(el[, 1], df$gp)],
                  xend = df$x[match(el[, 2], df$gp)], yend = df$y[match(el[, 2], df$gp)],
                  w = E(g)$w)
pal <- rep(grDevices::rainbow(12, end = 0.85), 10)[seq_along(levels(df$comm))]
p <- ggplot() +
  geom_segment(data = edf, aes(x, y, xend = xend, yend = yend, linewidth = w),
               color = "grey70", alpha = 0.5) +
  geom_point(data = df, aes(x, y, color = comm), size = 3.2) +
  geom_text_repel(data = df, aes(x, y, label = paste0(gp, "\n", top), color = comm),
                  size = 2.4, lineheight = 0.8, max.overlaps = 50,
                  segment.size = 0.2, show.legend = FALSE) +
  scale_color_manual(values = pal, guide = "none") +
  scale_linewidth(range = c(0.2, 1.8), guide = "none") +
  coord_equal() + theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
  labs(title = "GP-GP network (specificity-weighted shared genes)",
       subtitle = paste0(vcount(g), " GPs | edge weight = sum 1/df(shared gene), df<=",
         df_max, " | color = Louvain family"), x = NULL, y = NULL)
ggsave(paste0(out_dir, "gp_network_specific.png"), p, width = 11, height = 9, dpi = 150)
ggsave(paste0(out_dir, "gp_network_specific.pdf"), p, width = 11, height = 9)
saveRDS(list(graph = g, membership = memb, layout = df),
        paste0(out_dir, "gp_network_specific.rds"))
cat("\nsaved gp_network_specific.{png,pdf}\n")
