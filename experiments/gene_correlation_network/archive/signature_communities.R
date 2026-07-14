## Signature-gene correlation network + Louvain communities (DIAGNOSTIC, no layout).
## Goal: check whether communities are biologically coherent -- e.g. does Foxp3
## land in a community with other Treg genes, and does that community map to a
## single GP? Fast to run; if the biology checks out we invest in the layout.
##
## Signature genes = max per-GP-normalized |loading| > 0.1 (repo convention).
## Run from repo root:
##   Rscript experiments/gene_correlation_network/signature_communities.R

suppressPackageStartupMessages(library(igraph))
data_path <- "data/"; out_dir <- "experiments/gene_correlation_network/"
corr_thr <- 0.5; sig_thr <- 0.1; set.seed(1)

F_pm <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F_pm) <- paste0("GP", seq_len(ncol(F_pm)))
F_pm <- sweep(F_pm, 2, apply(abs(F_pm), 2, max), "/")
sig  <- rownames(F_pm)[apply(abs(F_pm), 1, max) > sig_thr]
F_sig <- F_pm[sig, ]; n_gp <- ncol(F_sig)
cat("signature genes:", length(sig), "\n")

dom_gp <- max.col(F_sig, ties.method = "first")            # dominant GP per gene
names(dom_gp) <- sig

# standardize rows, edges at corr >= thr
rm_ <- rowMeans(F_sig); rs_ <- apply(F_sig, 1, sd)
keep <- rs_ > 0; Z <- (F_sig[keep, ] - rm_[keep]) / rs_[keep]
genes <- rownames(Z); G <- nrow(Z)

C <- tcrossprod(Z) / (n_gp - 1)
C[lower.tri(C, diag = TRUE)] <- NA
hit <- which(C >= corr_thr, arr.ind = TRUE)
cat("edges at corr >=", corr_thr, ":", nrow(hit), "\n")

g <- graph_from_data_frame(
  data.frame(from = genes[hit[, 1]], to = genes[hit[, 2]]),
  directed = FALSE,
  vertices = data.frame(name = genes, dom_gp = dom_gp[genes]))
g <- simplify(g)
cat("graph:", vcount(g), "nodes,", ecount(g), "edges | components:",
    components(g)$no, "| largest:", max(components(g)$csize), "\n")

# Louvain communities
cl <- cluster_louvain(g)
cat("Louvain communities:", length(cl), "| modularity:",
    round(modularity(cl), 3), "\n")
memb <- membership(cl)
sz <- sort(table(memb), decreasing = TRUE)
cat("\ntop 15 community sizes:\n"); print(head(sz, 15))

# For each of the large communities: dominant GP purity + top-degree genes
deg <- degree(g)
cat("\n== large communities (size>=20): dominant-GP composition ==\n")
big <- as.integer(names(sz))[sz >= 20]
for (c in big[1:min(12, length(big))]) {
  mem <- names(memb)[memb == c]
  gp_tab <- sort(table(dom_gp[mem]), decreasing = TRUE)
  top_gp <- names(gp_tab)[1]
  purity <- round(gp_tab[1] / length(mem), 2)
  hub <- names(sort(deg[mem], decreasing = TRUE))[1:min(8, length(mem))]
  cat(sprintf("comm %2d (n=%3d): top GP=%s (%.0f%% of members) | hubs: %s\n",
              c, length(mem), top_gp, 100 * purity, paste(hub, collapse = ", ")))
}

# Foxp3 community
if ("Foxp3" %in% genes) {
  fc <- memb["Foxp3"]
  mem <- names(memb)[memb == fc]
  cat("\n== Foxp3 community (comm", fc, ", n=", length(mem), ") ==\n")
  treg <- c("Foxp3","Il2ra","Ctla4","Ikzf2","Ikzf4","Tnfrsf18","Il10","Nrp1",
            "Tnfrsf4","Lrrc32","Ikzf1","Tnfrsf9","Il2rb","Gpr83","Cd83","Izumo1r")
  cat("known Treg genes in this community:",
      paste(intersect(treg, mem), collapse = ", "), "\n")
  cat("dominant GP of Foxp3:", dom_gp["Foxp3"], "\n")
  gp_tab <- sort(table(dom_gp[mem]), decreasing = TRUE)
  cat("community dominant-GP breakdown (top 5):\n"); print(head(gp_tab, 5))
  cat("all members (top 40 by degree):\n")
  print(names(sort(deg[mem], decreasing = TRUE))[1:min(40, length(mem))])
} else cat("\nFoxp3 not in signature set / graph\n")

saveRDS(list(graph = g, membership = memb, dom_gp = dom_gp),
        paste0(out_dir, "signature_communities.rds"))
cat("\nsaved signature_communities.rds\n")
