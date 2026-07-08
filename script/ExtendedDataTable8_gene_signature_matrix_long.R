# Extended Data Table 8 (long format): comprehensive gene signature matrix.
#
# The same signature genes as the wide-format version
# (script/ExtendedDataTable8_gene_signature_matrix_wide.R), reshaped as a tidy
# long table: one row per gene, with gene_symbol, signature_name (GP +
# direction, underscore-joined, e.g. "GP1_up"), and score. Every gene has
# |score| > 0.1 on the same max|.|=1-per-GP-scaled gene factor matrix Table 1
# uses, ranked by |score| within direction and silently capped at the top 100
# up- and top 100 down-regulated genes per GP (no truncation-count row here,
# unlike the wide version -- see
# code/R/gene_signature_helpers.R::build_gp_gene_signature_blocks()'s
# `annotate_truncation` argument).
#
# Required inputs (data/) -- see code/README.md's "Data provenance" table.

data_path <- "data/"
output_path <- "figures/generated/"

source("code/R/gene_signature_helpers.R")

F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
# Normalize so each GP column has max|score| = 1 (same normalization Extended
# Data Table 1 uses for its gene signatures).
F_pm_filtered <- apply(F_pm_filtered, 2, function(x) x / max(abs(x)))
# F_pm_filtered's raw columns ("F1".."F200") are already in the same factor
# order as L_pm_filtered's "K1".."K200" (same underlying flashier fit --
# 01b_filter_cells.R only filters L's rows/cells, never F's columns), so
# column i is simply GPi; no cross-matrix name matching needed here.
n_gp <- ncol(F_pm_filtered)
gp_labels <- paste0("GP", seq_len(n_gp))

gp_blocks <- build_gp_gene_signature_blocks(F_pm_filtered, cutoff = 0.1, cap = 100, annotate_truncation = FALSE)
long_table <- do.call(rbind, Map(function(gp, block) {
  data.frame(
    gene_symbol = block$Gene,
    signature_name = paste0(gp, "_", block$Direction),
    score = block$Score,
    stringsAsFactors = FALSE
  )
}, gp_labels, gp_blocks))
rownames(long_table) <- NULL

write.csv(
  long_table,
  file = paste0(output_path, "ExtendedDataTable8_gene_signature_matrix_long.csv"),
  row.names = FALSE
)
