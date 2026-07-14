## Signature-gene correlation network -- LAYOUT + figure.
## 8053 signature genes (max per-GP-normalized |loading| > 0.1), edges at
## corr >= 0.5 on GP-loading profiles. FR force-directed layout; nodes colored
## by dominant GP; edges drawn faintly (subsampled for rendering). Landmark
## Treg genes annotated to show a biologically coherent cluster.
##
## Saves fr_layout.rds so the figure can be re-styled without recomputing.
## Run from repo root:
##   Rscript experiments/gene_correlation_network/signature_layout.R

suppressPackageStartupMessages({
  library(igraph); library(ggplot2); library(ggrastr)
})
data_path <- "data/"; out_dir <- "experiments/gene_correlation_network/"
set.seed(1)

x <- readRDS(paste0(out_dir, "signature_communities.rds"))
g <- x$graph; dom_gp <- x$dom_gp
cat("graph:", vcount(g), "nodes,", ecount(g), "edges\n")

# ── layout (reuse if cached) ──────────────────────────────────────────────────
lay_file <- paste0(out_dir, "fr_layout.rds")
if (file.exists(lay_file)) {
  lay <- readRDS(lay_file)
} else {
  cat("running FR layout...\n"); t0 <- Sys.time()
  lay <- layout_with_fr(g, niter = 500)
  cat("layout done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  rownames(lay) <- V(g)$name
  saveRDS(lay, lay_file)
}

df <- data.frame(x = lay[, 1], y = lay[, 2],
                 gene = V(g)$name, dom_gp = V(g)$dom_gp)
rownames(df) <- df$gene

# ── edges for rendering: subsample to keep the figure legible ─────────────────
el <- as_edgelist(g, names = TRUE)
set.seed(2)
sub <- if (nrow(el) > 60000) sample(nrow(el), 60000) else seq_len(nrow(el))
edf <- data.frame(x = df[el[sub, 1], "x"], y = df[el[sub, 1], "y"],
                  xend = df[el[sub, 2], "x"], yend = df[el[sub, 2], "y"])

# ── Treg landmark genes to annotate ───────────────────────────────────────────
treg <- c("Foxp3","Il2ra","Ctla4","Ikzf2","Nrp1","Lrrc32","Gpr83","Il10")
tdf  <- df[intersect(treg, df$gene), ]

p <- ggplot() +
  rasterise(geom_segment(data = edf,
              aes(x = x, y = y, xend = xend, yend = yend),
              color = "grey80", linewidth = 0.08, alpha = 0.15), dpi = 150) +
  rasterise(geom_point(data = df, aes(x, y, color = dom_gp),
              size = 0.5, alpha = 0.85), dpi = 150) +
  geom_point(data = tdf, aes(x, y), shape = 21, size = 2.2,
             color = "black", fill = NA, stroke = 0.6) +
  ggrepel::geom_text_repel(data = tdf, aes(x, y, label = gene),
             size = 3, fontface = "bold", max.overlaps = 20,
             min.segment.length = 0, segment.size = 0.3) +
  scale_color_gradientn(colors = grDevices::rainbow(200, end = 0.9),
              limits = c(1, 200), breaks = c(1, 50, 100, 150, 200),
              name = "dominant GP") +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), legend.key.height = unit(1.1, "cm")) +
  labs(title = "Signature-gene correlation network (corr >= 0.5 on GP loadings)",
       subtitle = paste0(vcount(g), " genes | FR layout | nodes = dominant GP | Treg landmarks labeled"),
       x = NULL, y = NULL)
ggsave(paste0(out_dir, "signature_network_domGP.pdf"), p, width = 9, height = 7.5)
ggsave(paste0(out_dir, "signature_network_domGP.png"), p, width = 9, height = 7.5, dpi = 150)
cat("saved signature_network_domGP.{pdf,png}\n")
