## Exploratory: gauge the scale of the gene-gene correlation network before
## committing to a layout. Correlation is computed on the per-GP (per-column)
## max-abs normalized F matrix; each gene = its 200-dim GP loading vector.
## Reports: nonzero-sd genes, edge counts at several corr thresholds, degree
## distribution, connected-component structure at corr >= 0.5.
## Run from repo root: Rscript experiments/gene_correlation_network/explore_scale.R

data_path <- "data/"

F_pm <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
colnames(F_pm) <- paste0("GP", seq_len(ncol(F_pm)))
cat("F matrix:", nrow(F_pm), "genes x", ncol(F_pm), "GPs\n")

## per-GP (per-column) max-abs normalization
F_pm <- sweep(F_pm, 2, apply(abs(F_pm), 2, max), "/")

## standardize each gene (row) over the GP dimensions: z-score
n_gp <- ncol(F_pm)
rmean <- rowMeans(F_pm)
rsd   <- apply(F_pm, 1, sd)
keep  <- rsd > 0
cat("genes with nonzero sd (can have edges):", sum(keep),
    "| all-zero/constant genes:", sum(!keep), "\n")

Z <- (F_pm[keep, ] - rmean[keep]) / rsd[keep]   # rows standardized
G <- nrow(Z)
cat("standardized matrix:", G, "x", n_gp, "\n")

## chunked upper-triangle edge counting at several thresholds
thr_grid <- c(0.3, 0.4, 0.5, 0.6, 0.7)
edge_count <- setNames(numeric(length(thr_grid)), thr_grid)
deg <- setNames(numeric(G), rownames(Z))          # degree at corr >= 0.5
chunk <- 1000
for (start in seq(1, G, by = chunk)) {
  idx <- start:min(start + chunk - 1, G)
  C   <- tcrossprod(Z[idx, , drop = FALSE], Z) / (n_gp - 1)   # |idx| x G
  # only count j > i to avoid double counting / self
  for (r in seq_along(idx)) {
    i <- idx[r]
    row <- C[r, ]
    row[1:i] <- NA                                 # keep j > i
    for (t in thr_grid) edge_count[as.character(t)] <-
      edge_count[as.character(t)] + sum(row >= t, na.rm = TRUE)
    # degree at 0.5 (count both directions later; here count j>i and add)
  }
  # degree at 0.5 for this chunk vs all (exclude self)
  C05 <- C >= 0.5
  C05[cbind(seq_along(idx), idx)] <- FALSE
  deg[idx] <- deg[idx] + rowSums(C05, na.rm = TRUE)
}
cat("\n== edge counts (undirected, corr >= thr) ==\n")
print(edge_count)

cat("\n== degree distribution at corr >= 0.5 ==\n")
print(summary(deg))
cat("genes with degree >= 1:", sum(deg >= 1), "\n")
cat("top-degree genes:\n")
print(head(sort(deg, decreasing = TRUE), 15))

saveRDS(list(Z = Z, deg = deg, keep = keep, edge_count = edge_count),
        file = "experiments/gene_correlation_network/explore_scale.rds")
cat("\nsaved explore_scale.rds\n")
