## Signature-gene network on a UMAP embedding.
## Force layouts (FR/DrL) collapse this small-world graph into a blob. Instead we
## embed genes with UMAP on their GP-loading profiles (cosine metric == Pearson
## corr on the row-standardized matrix), which spreads coherent gene sets into
## separated islands, then overlay the kNN correlation edges on top -- a clean
## "network" that actually shows community structure. Nodes colored by dominant
## GP; Treg landmarks annotated.
## Run from repo root:
##   Rscript experiments/gene_correlation_network/umap_network.R

suppressPackageStartupMessages({
  library(uwot); library(igraph); library(ggplot2); library(ggrastr)
})
data_path <- "data/"; out_dir <- "experiments/gene_correlation_network/"
sig_thr <- 0.1; set.seed(1)

F_pm <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F_pm) <- paste0("GP", seq_len(ncol(F_pm)))
F_pm <- sweep(F_pm, 2, apply(abs(F_pm), 2, max), "/")
sig  <- rownames(F_pm)[apply(abs(F_pm), 1, max) > sig_thr]
F_sig <- F_pm[sig, ]; n_gp <- ncol(F_sig)
dom_gp <- max.col(F_sig, ties.method = "first"); names(dom_gp) <- sig

# row-standardize so cosine == Pearson correlation between genes
rm_ <- rowMeans(F_sig); rs_ <- apply(F_sig, 1, sd)
keep <- rs_ > 0; Z <- (F_sig[keep, ] - rm_[keep]) / rs_[keep]
genes <- rownames(Z)
cat("signature genes:", nrow(Z), "\n")

# ── UMAP embedding (cosine metric) ────────────────────────────────────────────
cat("UMAP...\n"); t0 <- Sys.time()
emb <- umap(Z, n_neighbors = 15, min_dist = 0.3, metric = "cosine",
            n_threads = 4, verbose = FALSE)
cat("done in", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")
df <- data.frame(x = emb[, 1], y = emb[, 2], gene = genes,
                 dom_gp = dom_gp[genes]); rownames(df) <- genes
saveRDS(df, paste0(out_dir, "umap_coords.rds"))

# ── kNN correlation edges to overlay (reuse if present) ───────────────────────
kf <- paste0(out_dir, "knn_layout.rds")
edf <- NULL
if (file.exists(kf)) {
  g <- readRDS(kf)$graph
  el <- as_edgelist(g, names = TRUE)
  el <- el[el[, 1] %in% genes & el[, 2] %in% genes, ]
  edf <- data.frame(x = df[el[, 1], "x"], y = df[el[, 1], "y"],
                    xend = df[el[, 2], "x"], yend = df[el[, 2], "y"])
  cat("overlaying", nrow(edf), "kNN edges\n")
}

treg <- c("Foxp3","Il2ra","Ctla4","Ikzf2","Ikzf4","Tnfrsf18","Il10","Nrp1",
          "Tnfrsf4","Lrrc32","Ikzf1","Gpr83","Izumo1r","Cish","Lta")
treg_df <- df[intersect(treg, df$gene), ]

p <- ggplot()
if (!is.null(edf))
  p <- p + rasterise(geom_segment(data = edf, aes(x, y, xend = xend, yend = yend),
              color = "grey80", linewidth = 0.08, alpha = 0.2), dpi = 150)
p <- p +
  rasterise(geom_point(data = df, aes(x, y, color = dom_gp),
              size = 0.55, alpha = 0.9), dpi = 150) +
  geom_point(data = treg_df, aes(x, y), shape = 21, size = 2.4,
             color = "black", fill = NA, stroke = 0.7) +
  ggrepel::geom_text_repel(data = treg_df, aes(x, y, label = gene),
             size = 3, fontface = "bold", max.overlaps = 30,
             min.segment.length = 0, segment.size = 0.3) +
  scale_color_gradientn(colors = grDevices::rainbow(200, end = 0.9),
              limits = c(1, 200), breaks = c(1, 50, 100, 150, 200), name = "dominant GP") +
  coord_equal() + theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), legend.key.height = unit(1.1, "cm")) +
  labs(title = "Signature-gene correlation network on UMAP (cosine = Pearson corr)",
       subtitle = paste0(nrow(df), " signature genes | UMAP of GP loadings + kNN edges | nodes = dominant GP | Treg labeled"),
       x = NULL, y = NULL)
ggsave(paste0(out_dir, "umap_network_domGP.png"), p, width = 9, height = 7.5, dpi = 150)
ggsave(paste0(out_dir, "umap_network_domGP.pdf"), p, width = 9, height = 7.5)
cat("saved umap_network_domGP.{png,pdf}\n")
