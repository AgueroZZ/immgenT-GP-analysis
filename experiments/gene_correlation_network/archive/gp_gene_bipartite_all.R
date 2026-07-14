## Fig-3d-style GP-gene signature network, but for ALL 200 GPs (broad picture).
##
## Same visual language as Figure 3d: GP nodes (squares) linked to their top-5
## positive (red edges) and top-5 negative (blue edges) signature genes, gene
## nodes (circles), ggraph "stress" layout. Extending to all 200 GPs turns the
## shared genes into hub nodes that pull related GPs together, so GP families
## emerge from the layout itself. GP squares are colored by Louvain community
## (detected on this bipartite graph); the most-shared hub genes are labeled.
##
## Normalization matches code/R/activation_shared_setup.R:
##   F_pm_filtered_norm = per-column (per-GP) max-abs.
##
## Run from repo root:
##   Rscript experiments/gene_correlation_network/gp_gene_bipartite_all.R [drop_df]
##   drop_df (optional): drop genes shared by > drop_df GPs to declutter hubs.

suppressPackageStartupMessages({
  library(tidygraph); library(ggraph); library(igraph); library(ggplot2)
})
data_path <- "data/"; out_dir <- "experiments/gene_correlation_network/"
args <- commandArgs(trailingOnly = TRUE)
drop_df <- if (length(args) >= 1) as.integer(args[1]) else Inf   # Inf = keep all
set.seed(42)

F <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F) <- paste0("GP", seq_len(ncol(F)))
F <- sweep(F, 2, apply(abs(F), 2, max), "/")     # per-GP max-abs (== _norm)
GPs <- colnames(F)

# top 5 UP + top 5 DOWN per GP (|loading| >= 0.1): richer sharing than the
# top-5-by-abs rule (which leaves many GPs as isolated stars) for a 200-GP view
top_pos <- lapply(GPs, function(j) { x <- F[, j]
  u <- names(sort(x, decreasing = TRUE))[1:5]; u[x[u] >=  0.1] })
top_neg <- lapply(GPs, function(j) { x <- F[, j]
  d <- names(sort(x))[1:5];                    d[x[d] <= -0.1] })
names(top_pos) <- names(top_neg) <- GPs

mk <- function(lst, type, col) do.call(rbind, lapply(GPs, function(g)
  if (length(lst[[g]])) data.frame(GP = g, Gene = lst[[g]], Type = type, Color = col)
  else NULL))
edges <- rbind(mk(top_pos, "Positive", "red"), mk(top_neg, "Negative", "blue"))
edges <- edges[edges$Gene != "" & !is.na(edges$Gene), ]

# optional: drop promiscuous hub genes to declutter
df_gene <- table(edges$Gene)
if (is.finite(drop_df)) {
  drop <- names(df_gene)[df_gene > drop_df]
  edges <- edges[!edges$Gene %in% drop, ]
  cat("dropped", length(drop), "genes shared by >", drop_df, "GPs\n")
}
cat("edges:", nrow(edges), "| GP nodes:", length(unique(edges$GP)),
    "| gene nodes:", length(unique(edges$Gene)), "\n")

g <- as_tbl_graph(edges[, c("GP", "Gene", "Type", "Color")], directed = FALSE)
gene_names <- setdiff(unique(edges$Gene), unique(edges$GP))
deg <- igraph::degree(g)
comm <- membership(cluster_louvain(g))
cat("Louvain communities:", length(unique(comm)),
    "| modularity:", round(modularity(cluster_louvain(g)), 3), "\n")

g <- g %>% activate(nodes) %>% mutate(
  is_gp   = name %in% GPs,
  comm    = factor(comm[name]),
  gp_lab  = ifelse(name %in% GPs, sub("GP", "", name), ""),
  # label only hub genes shared by many GPs
  hub_lab = ifelse(!name %in% GPs & df_gene[name] >= 6, name, ""),
  gp_deg  = deg[name])

# color GP squares by community; gene circles stay grey
ncomm <- length(unique(comm))
pal <- rep(grDevices::rainbow(12, end = 0.88), length.out = ncomm)

p <- ggraph(g, layout = "stress") +
  geom_edge_link(aes(color = Color), alpha = 0.25, width = 0.35) +
  geom_node_point(aes(filter = !is_gp), shape = 16, size = 0.7,
                  color = "grey55", alpha = 0.6) +
  geom_node_point(aes(filter = is_gp, color = comm), shape = 15,
                  size = 3.2, alpha = 0.85) +
  geom_node_text(aes(filter = is_gp, label = gp_lab), color = "black",
                 size = 1.7, fontface = "bold") +
  geom_node_text(aes(filter = !is_gp, label = hub_lab), color = "grey20",
                 size = 2.2, repel = TRUE, max.overlaps = 30, segment.size = 0.15) +
  scale_edge_color_identity() +
  scale_color_manual(values = pal, guide = "none") +
  theme_void() +
  labs(title = "GP-gene signature network -- all 200 GPs (broad picture)",
       subtitle = paste0(length(unique(edges$GP)), " GPs linked to top-5 up (red) / top-5 down (blue) genes | stress layout | GP squares colored by community"),
       caption = "Gene hubs (shared by >=6 GPs) labeled") +
  theme(plot.title = element_text(face = "bold"))

ggsave(paste0(out_dir, "gp_gene_bipartite_all.png"), p, width = 13, height = 11, dpi = 150)
ggsave(paste0(out_dir, "gp_gene_bipartite_all.pdf"), p, width = 13, height = 11)
saveRDS(g, paste0(out_dir, "gp_gene_bipartite_all.rds"))
cat("saved gp_gene_bipartite_all.{png,pdf}\n")
