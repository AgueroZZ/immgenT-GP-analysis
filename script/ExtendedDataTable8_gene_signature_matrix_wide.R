# Extended Data Table 8 (wide format): comprehensive gene signature matrix.
#
# For every GP, its full signature gene list (not just the top 5 shown in
# Extended Data Table 1): every gene with |score| > 0.1 on the same
# max|.|=1-per-GP-scaled gene factor matrix Table 1 uses, ranked by |score|
# within direction and capped at the top 100 up- and top 100 down-regulated
# genes. When a direction has more than 100 qualifying genes, one extra row
# right after the (capped) gene rows -- e.g. row 101 after 100 up-regulated
# genes -- notes how many more were left out (see
# code/R/gene_signature_helpers.R::build_gp_gene_signature_blocks()), so a GP
# contributes at most 202 rows. Different GPs have different signature-gene
# counts, so this is a ragged "wide" table: one 3-column (Gene, Direction,
# Score) block per GP, side by side, padded with blank cells up to the
# tallest block. A flat CSV can't express that per-GP grouping, so the output
# is an .xlsx with a merged GP header spanning each block.
#
# This is the wide, per-GP-column layout; script/ExtendedDataTable8_gene_signature_matrix_long.R
# produces the same signature genes as a tidy long table (one row per gene)
# instead.
#
# Required inputs (data/) -- see code/README.md's "Data provenance" table.

data_path <- "data/"
output_path <- "figures/generated/"

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package 'openxlsx' is required. Please install it with install.packages('openxlsx').")
}

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

gp_blocks <- build_gp_gene_signature_blocks(F_pm_filtered, cutoff = 0.1, cap = 100)
max_rows <- max(vapply(gp_blocks, nrow, integer(1)))

pad_block <- function(block, n) {
  if (nrow(block) < n) {
    block <- rbind(block, data.frame(
      Gene = rep(NA_character_, n - nrow(block)),
      Direction = rep(NA_character_, n - nrow(block)),
      Score = rep(NA_real_, n - nrow(block)),
      stringsAsFactors = FALSE
    ))
  }
  block
}
wide_df <- do.call(cbind, lapply(gp_blocks, pad_block, n = max_rows))

# Write as .xlsx with a merged GP header spanning each Gene/Direction/Score
# triplet -- a flat CSV can't express this meta-column grouping.
wb <- openxlsx::createWorkbook()
sheet <- "Gene signatures"
openxlsx::addWorksheet(wb, sheet)

for (i in seq_len(n_gp)) {
  col_start <- (i - 1) * 3 + 1
  openxlsx::mergeCells(wb, sheet, cols = col_start:(col_start + 2), rows = 1)
  openxlsx::writeData(wb, sheet, gp_labels[i], startCol = col_start, startRow = 1, colNames = FALSE)
}
openxlsx::writeData(
  wb, sheet,
  matrix(rep(c("Gene", "Direction", "Score"), n_gp), nrow = 1),
  startCol = 1, startRow = 2, colNames = FALSE
)
openxlsx::writeData(wb, sheet, wide_df, startCol = 1, startRow = 3, colNames = FALSE, na.string = "")

bold_center <- openxlsx::createStyle(textDecoration = "bold", halign = "center")
bold <- openxlsx::createStyle(textDecoration = "bold")
openxlsx::addStyle(wb, sheet, bold_center, rows = 1, cols = seq_len(n_gp * 3), gridExpand = TRUE)
openxlsx::addStyle(wb, sheet, bold, rows = 2, cols = seq_len(n_gp * 3), gridExpand = TRUE)
openxlsx::freezePane(wb, sheet, firstActiveRow = 3, firstActiveCol = 1)

openxlsx::saveWorkbook(
  wb,
  file = paste0(output_path, "ExtendedDataTable8_gene_signature_matrix_wide.xlsx"),
  overwrite = TRUE
)
