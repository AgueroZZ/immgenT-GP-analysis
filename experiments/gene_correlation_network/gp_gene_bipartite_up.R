## Fig-3d-style GP-gene bipartite network, UP-regulated genes only, with marker
## genes highlighted. Broad picture over all 200 GPs.
##
## Design choices vs the archived all-edges version:
##   * UP genes only  -> half the nodes, no blue edges, far less clutter.
##   * drop promiscuous "glue" genes (shared by > drop_df GPs) -> GP families
##     stop collapsing into one central blob and spread out.
##   * GP squares colored by Louvain family, sized by PVE.
##   * a curated MARKER set (CD8, FR4/Izumo1r, Foxp3, ...) drawn large + labeled;
##     all other gene dots stay small/grey so the anchors pop.
##
## Run from repo root:
##   Rscript experiments/gene_correlation_network/gp_gene_bipartite_up.R [n_top] [drop_df]

suppressPackageStartupMessages({
  library(tidygraph); library(ggraph); library(igraph); library(ggplot2)
})
data_path <- "data/"; out_dir <- "experiments/gene_correlation_network/"
args <- commandArgs(trailingOnly = TRUE)
n_top   <- if (length(args) >= 1) as.integer(args[1]) else 10
drop_df <- if (length(args) >= 2) as.numeric(args[2]) else Inf   # Inf = keep all genes
set.seed(42)

# ── marker genes to highlight (edit freely) ───────────────────────────────────
MARKERS <- c("Foxp3", "Cd8a", "Cd8b1", "Izumo1r", "Il2ra")     # FR4 = Izumo1r

F <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F) <- paste0("GP", seq_len(ncol(F)))
F <- sweep(F, 2, apply(abs(F), 2, max), "/")
GPs <- colnames(F)

# top n_top UP genes per GP (loading >= 0.1)
top_up <- lapply(GPs, function(j) { x <- F[, j]
  u <- names(sort(x, decreasing = TRUE))[1:n_top]; u[x[u] >= 0.1] })
names(top_up) <- GPs

edges <- do.call(rbind, lapply(GPs, function(g)
  if (length(top_up[[g]])) data.frame(GP = g, Gene = top_up[[g]]) else NULL))

# drop promiscuous glue genes (but never drop a marker)
df_gene <- table(edges$Gene)
drop <- setdiff(names(df_gene)[df_gene > drop_df], MARKERS)
edges <- edges[!edges$Gene %in% drop, ]
cat("dropped", length(drop), "glue genes (df >", drop_df, "); edges:", nrow(edges),
    "| GP:", length(unique(edges$GP)), "| genes:", length(unique(edges$Gene)), "\n")
present <- intersect(MARKERS, edges$Gene)
cat("markers present as up-genes:", paste(present, collapse = ", "), "\n")

# PVE for GP node size
pve <- setNames(rep(NA_real_, 200), paste0("GP", 1:200))
pve_file <- "data/factors_comments_updated.xlsx - factor_pve.csv"
if (file.exists(pve_file)) {
  pt <- read.csv(pve_file, check.names = FALSE)
  pve[paste0("GP", as.integer(gsub("F", "", pt$factor)))] <- pt$pve
}

g0 <- graph_from_data_frame(edges, directed = FALSE)
gp_set0 <- unique(edges$GP)
# keep only components that contain >= 2 GPs (drop lonely single-GP stars)
comp <- components(g0)
gp_per_comp <- tapply(names(comp$membership) %in% gp_set0, comp$membership, sum)
keep_comp <- as.integer(names(gp_per_comp)[gp_per_comp >= 2])
keep_v <- names(comp$membership)[comp$membership %in% keep_comp]
n_drop_gp <- length(gp_set0) - sum(keep_v %in% gp_set0)
cat("dropped", n_drop_gp, "single-GP star components;", length(keep_v), "nodes kept\n")
gi <- igraph::induced_subgraph(g0, keep_v)
# hide private (degree-1) gene dots: keep only genes shared by >=2 GPs, plus all
# markers and all GPs -- leaves just the connective backbone linking GP families
gdeg   <- igraph::degree(gi)
vnames <- igraph::V(gi)$name
rm_priv <- vnames[!(vnames %in% gp_set0) & gdeg <= 1 & !(vnames %in% MARKERS)]
gi <- igraph::delete_vertices(gi, rm_priv)
cat("hid", length(rm_priv), "private gene dots; kept shared genes + markers\n")
g <- as_tbl_graph(gi, directed = FALSE)
gp_set <- intersect(gp_set0, igraph::V(gi)$name)
comm <- membership(cluster_louvain(g))
cat("Louvain families:", length(unique(comm)), "| modularity:",
    round(modularity(cluster_louvain(g)), 3), "\n")

# ── one distinct highlight color per marker (shown via a legend, not text). ────
present <- intersect(MARKERS, igraph::V(gi)$name)
distinct_pal <- rep(c(
  "#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#1ba3c4","#f032e6",
  "#469990","#9A6324","#800000","#808000","#000075","#ff4500","#0a9396",
  "#556b2f","#ff1493","#1e90ff","#7f00ff","#00b300","#b30086","#005c31",
  "#a05a2c","#c800c8","#8a6d00"), length.out = max(length(present), 1))
mkcol    <- setNames(distinct_pal[seq_along(present)], present)
GP_COL   <- "#5B8CB8"   # uniform GP color (blue)
GENE_COL <- "#C7A76C"   # shared-gene color (tan) -- distinct from GPs & markers

g <- g %>% activate(nodes) %>% mutate(
  is_gp   = name %in% gp_set,
  gp_lab  = ifelse(name %in% gp_set, sub("GP", "", name), ""),
  marker  = !name %in% gp_set & name %in% present,
  mk_name = ifelse(!name %in% gp_set & name %in% present, name, NA_character_),
  gp_size = ifelse(name %in% gp_set, 5.5, NA))   # uniform GP ball size

# marker-incident edges carry that marker's color; the rest stay light
g <- g %>% activate(edges) %>% mutate(
  gene_end = ifelse(.N()$marker[from], .N()$name[from],
              ifelse(.N()$marker[to], .N()$name[to], NA_character_)),
  is_mk    = !is.na(gene_end),
  edge_col = ifelse(is_mk, mkcol[gene_end], NA_character_))

# ── FR layout with radial compression: pull outer nodes in (shorten long
#    peripheral edges, reclaim whitespace) while the dense core keeps its spread.
set.seed(2)
lay <- igraph::layout_with_fr(gi)
cx <- mean(lay[, 1]); cy <- mean(lay[, 2])
rr <- sqrt((lay[, 1] - cx)^2 + (lay[, 2] - cy)^2)
th <- atan2(lay[, 2] - cy, lay[, 1] - cx)
rr <- rr^0.65                                   # <1 exponent = compress outward
lay <- cbind(x = cx + rr * cos(th), y = cy + rr * sin(th))
# jitter to separate structural twins (GPs with identical neighborhoods that FR
# stacks on the exact same spot)
jit <- 0.012 * max(diff(range(lay[, 1])), diff(range(lay[, 2])))
lay[, 1] <- lay[, 1] + rnorm(nrow(lay), 0, jit)
lay[, 2] <- lay[, 2] + rnorm(nrow(lay), 0, jit)
rownames(lay) <- igraph::V(gi)$name
node_names <- g %>% activate(nodes) %>% pull(name)
LX <- lay[node_names, "x"]; LY <- lay[node_names, "y"]

p <- ggraph(g, layout = "manual", x = LX, y = LY) +
  geom_edge_link(aes(filter = !is_mk), color = "black", alpha = 0.22, width = 0.32) +
  geom_edge_link(aes(filter = is_mk, edge_colour = edge_col), alpha = 0.9, width = 0.7) +
  # shared-gene dots (tan)
  geom_node_point(aes(filter = !is_gp & !marker), shape = 16, size = 1.0,
                  color = GENE_COL, alpha = 0.75) +
  # GP balls: uniform blue, thin white outline
  geom_node_point(aes(filter = is_gp, size = gp_size),
                  shape = 21, fill = GP_COL, color = "white", stroke = 0.4) +
  geom_node_text(aes(filter = is_gp, label = gp_lab), color = "grey97",
                 size = 1.9) +
  # marker genes: colored balls, identity -> legend (no on-plot text)
  geom_node_point(aes(filter = marker, color = mk_name), shape = 16, size = 3.2) +
  scale_color_manual(values = mkcol, breaks = present, name = "Marker gene",
                     na.translate = FALSE) +
  scale_edge_colour_identity() +
  scale_size_identity() + scale_edge_width_identity() +
  scale_x_continuous(expand = expansion(mult = 0.03)) +
  scale_y_continuous(expand = expansion(mult = 0.03)) +
  theme_void(base_size = 12) +
  theme(plot.margin = margin(4, 4, 4, 4),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        legend.text  = element_text(face = "italic", size = 10)) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(paste0(out_dir, "gp_gene_bipartite_up.png"), p, width = 17, height = 14, dpi = 160, bg = "white", limitsize = FALSE)
ggsave(paste0(out_dir, "gp_gene_bipartite_up.pdf"), p, width = 17, height = 14, bg = "white", limitsize = FALSE)
cat("saved gp_gene_bipartite_up.{png,pdf}\n")
