## GP-centric network: nodes = GPs, edges = shared top signature genes.
##
## Each GP is summarized by its top 5 up + top 5 down genes (per-GP max-abs
## normalized loading, |loading| >= 0.1). Two GPs are linked if their top-gene
## sets overlap; edge weight = number of shared genes. GPs that share marker
## genes form communities (GP "families"). Only ~200 nodes -> a clean, layout-
## able, structured community plot (unlike the 8k-gene hairball).
##
## Run from repo root:
##   Rscript experiments/gene_correlation_network/gp_network.R [min_shared]

suppressPackageStartupMessages({
  library(igraph); library(ggplot2); library(ggrepel)
})
data_path <- "data/"; out_dir <- "experiments/gene_correlation_network/"
args <- commandArgs(trailingOnly = TRUE)
min_shared <- if (length(args) >= 1) as.integer(args[1]) else 2
set.seed(1)

F_pm <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F_pm) <- paste0("GP", seq_len(ncol(F_pm)))
F_pm <- sweep(F_pm, 2, apply(abs(F_pm), 2, max), "/")
nGP <- ncol(F_pm)

# top 5 up + 5 down per GP (|loading| >= 0.1)
sets <- lapply(seq_len(nGP), function(j) {
  v  <- F_pm[, j]
  up <- names(sort(v, decreasing = TRUE))[1:5]; up <- up[v[up] >=  0.1]
  dn <- names(sort(v))[1:5];                    dn <- dn[v[dn] <= -0.1]
  unique(c(up, dn))
})
names(sets) <- paste0("GP", seq_len(nGP))
top_gene <- sapply(seq_len(nGP), function(j) names(which.max(F_pm[, j])))

# optional PVE for node size
pve <- setNames(rep(NA_real_, nGP), paste0("GP", seq_len(nGP)))
pve_file <- "data/factors_comments_updated.xlsx - factor_pve.csv"
if (file.exists(pve_file)) {
  pt <- read.csv(pve_file, check.names = FALSE)
  pve[paste0("GP", as.integer(gsub("F", "", pt$factor)))] <- pt$pve
}

# GP-GP shared-gene edges
ed <- data.frame(from = character(), to = character(), w = integer())
for (i in 1:(nGP - 1)) for (j in (i + 1):nGP) {
  s <- length(intersect(sets[[i]], sets[[j]]))
  if (s >= min_shared) ed <- rbind(ed, data.frame(
    from = paste0("GP", i), to = paste0("GP", j), w = s))
}
cat("edges (>=", min_shared, "shared genes):", nrow(ed), "\n")

g <- graph_from_data_frame(ed, directed = FALSE,
       vertices = data.frame(name = paste0("GP", seq_len(nGP)),
                             pve = pve, top = top_gene))
g <- delete_vertices(g, degree(g) == 0)
cat("graph:", vcount(g), "GP nodes,", ecount(g), "edges | components:",
    components(g)$no, "| largest:", max(components(g)$csize), "\n")

# weighted Louvain communities
cl <- cluster_louvain(g, weights = E(g)$w)
memb <- membership(cl)
cat("Louvain:", length(cl), "GP-families | modularity:", round(modularity(cl), 3), "\n")

# describe each community: member GPs + most-shared genes among them
cat("\n== GP families ==\n")
for (c in sort(unique(memb))) {
  gps <- names(memb)[memb == c]
  gtab <- sort(table(unlist(sets[gps])), decreasing = TRUE)
  hallmark <- names(gtab)[gtab >= 2]; if (!length(hallmark)) hallmark <- names(head(gtab, 3))
  cat(sprintf("family %d (%d GPs): %s\n   shared genes: %s\n", c, length(gps),
              paste(gps, collapse = ", "), paste(head(hallmark, 12), collapse = ", ")))
}

# where is the Treg GP (Foxp3's dominant GP)?
foxp3_gp <- paste0("GP", which.max(F_pm["Foxp3", ]))
if (foxp3_gp %in% names(memb)) {
  fam <- names(memb)[memb == memb[foxp3_gp]]
  cat("\nFoxp3 dominant GP =", foxp3_gp, "| its family:", paste(fam, collapse = ", "), "\n")
  cat("  neighbors of", foxp3_gp, ":", paste(names(neighbors(g, foxp3_gp)), collapse = ", "), "\n")
}

# ── layout + plot ─────────────────────────────────────────────────────────────
lay <- layout_with_fr(g, weights = E(g)$w, niter = 3000)
df <- data.frame(x = lay[, 1], y = lay[, 2], gp = V(g)$name,
                 comm = factor(memb[V(g)$name]),
                 top = V(g)$top, pve = V(g)$pve)
el <- as_edgelist(g, names = TRUE)
edf <- data.frame(x = df$x[match(el[, 1], df$gp)], y = df$y[match(el[, 1], df$gp)],
                  xend = df$x[match(el[, 2], df$gp)], yend = df$y[match(el[, 2], df$gp)],
                  w = E(g)$w)
df$size <- if (all(is.na(df$pve))) 3 else 2 + 6 * sqrt(pmax(df$pve, 0) / max(df$pve, na.rm = TRUE))

pal <- rep(grDevices::rainbow(12, end = 0.85), 10)[seq_along(levels(df$comm))]
p <- ggplot() +
  geom_segment(data = edf, aes(x, y, xend = xend, yend = yend, linewidth = w),
               color = "grey70", alpha = 0.5) +
  geom_point(data = df, aes(x, y, color = comm, size = size)) +
  geom_text_repel(data = df, aes(x, y, label = paste0(gp, "\n", top), color = comm),
                  size = 2.4, lineheight = 0.8, max.overlaps = 40,
                  segment.size = 0.2, show.legend = FALSE) +
  scale_color_manual(values = pal, guide = "none") +
  scale_linewidth(range = c(0.2, 1.6), guide = "none") +
  scale_size_identity() +
  coord_equal() + theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(title = "GP-GP network: GPs linked by shared top signature genes",
       subtitle = paste0(vcount(g), " GPs | edge = >=", min_shared,
         " shared top-genes (5 up/5 down) | node size ~ PVE | color = Louvain family"),
       x = NULL, y = NULL)
ggsave(paste0(out_dir, "gp_network.png"), p, width = 11, height = 9, dpi = 150)
ggsave(paste0(out_dir, "gp_network.pdf"), p, width = 11, height = 9)
saveRDS(list(graph = g, membership = memb, sets = sets, layout = df),
        paste0(out_dir, "gp_network.rds"))
cat("\nsaved gp_network.{png,pdf}\n")
