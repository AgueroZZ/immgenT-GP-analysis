# Extended Data Table 1: summary of GP characteristics and comments.
#
# One row per GP (GP1..GP200) with:
#   - Top Genes + / Top Genes -: the top 15 up- and top 15 down-regulated genes
#     by factor score, among those with |score| > 0.1 on the max|.|=1-per-GP
#     scaled gene factor matrix (F_pm_filtered). GPs with fewer than 15
#     qualifying genes in a direction list all of them.
#   - Comments: the manual, human-written description of the GP, read from the
#     curated curation/GP_manual_annotations.csv (blank where a GP has no
#     comment).
#   - Prop. active cells: the fraction of non-thymocyte cells whose loading
#     exceeds 0.1 after scaling the GP's loading column by its maximum -- the
#     quantity Figure 2b histograms.
#   - N active genes: the number of genes with |score| > 0.25 on the
#     max|.|=1-per-GP scaled gene factor matrix -- the quantity Figure 2c
#     histograms. Note this is a stricter cutoff than the 0.1 used for the
#     signature-gene columns above.
#
# Writes figures/final-selected/ExtendedDataTable1_GP_summary.csv.
#
# The full signature gene list behind the Top Genes columns (up to 100 per
# direction) is Extended Data Table 2.

data_path <- "data/"
curation_path <- "curation/"
output_path <- "figures/final-selected/"

# ---- Loadings and gene factor matrix ----
# L_pm_filtered's columns are the flashier fit's "K1".."K200"; F_pm_filtered's
# are "F1".."F200". Both are in the same factor order, so position -- not name
# -- is what aligns them, and both are relabeled "GP1".."GP200" for output.
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
stopifnot(ncol(L_pm_filtered) == ncol(F_pm_filtered))
gp_labels <- paste0("GP", seq_len(ncol(L_pm_filtered)))

# ---- Prop. active cells (Figure 2b) ----
# Non-thymocyte cells only, matching Figure 2; each loading column is scaled by
# its own maximum before the 0.1 cut.
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]
rm(seurat_meta); gc()
non_thymo_cells <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 != "thymocyte"]
L_non_thymo <- L_pm_filtered[non_thymo_cells, ]
L_norm_col <- L_non_thymo / matrix(apply(L_non_thymo, 2, max),
                                   nrow = nrow(L_non_thymo), ncol = ncol(L_non_thymo),
                                   byrow = TRUE)
prop_active_cells <- colSums(L_norm_col > 1e-1) / nrow(L_norm_col)
rm(L_pm_filtered, L_non_thymo, L_norm_col); gc()

# ---- N active genes (Figure 2c) and signature genes ----
F_pm_norm <- apply(F_pm_filtered, 2, function(x) x / max(abs(x)))
colnames(F_pm_norm) <- gp_labels
n_active_genes <- colSums(abs(F_pm_norm) > 0.25)

top_signatures <- function(score_mat, gps, cutoff = 0.1, n = 15) {
  pos <- lapply(gps, function(gp) {
    vals <- score_mat[, gp]
    cand <- names(vals)[vals > cutoff]
    paste(head(cand[order(vals[cand], decreasing = TRUE)], n), collapse = "; ")
  })
  neg <- lapply(gps, function(gp) {
    vals <- score_mat[, gp]
    cand <- names(vals)[vals < -cutoff]
    paste(head(cand[order(abs(vals[cand]), decreasing = TRUE)], n), collapse = "; ")
  })
  list(pos = unlist(pos), neg = unlist(neg))
}
gene_sig <- top_signatures(F_pm_norm, gp_labels)

# ---- Manual comments (curated input; blank where a GP has no comment) ----
# Keyed by the GP*n* output labels. Insisting on an exact GP1..GP200 match means
# a renamed or dropped row fails loudly here instead of silently blanking a
# comment in the published table.
manual_comments <- read.csv(paste0(curation_path, "GP_manual_annotations.csv"),
                            stringsAsFactors = FALSE, colClasses = "character")
stopifnot(identical(manual_comments$GP, gp_labels))

supp_table <- data.frame(
  GP = gp_labels,
  `Top Genes +` = gene_sig$pos,
  `Top Genes -` = gene_sig$neg,
  Comments = trimws(manual_comments$Comments),
  `Prop. active cells` = as.numeric(prop_active_cells),
  `N active genes` = as.integer(n_active_genes),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.csv(supp_table,
          file = paste0(output_path, "ExtendedDataTable1_GP_summary.csv"),
          row.names = FALSE)
