library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(fastTopics)
library(dplyr)
library(pheatmap)
library(scales)
library(scattermore)

data_path <- "../data/"
code_path <- "../code/"
data_path <- "data/"
code_path <- "code/"

seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
cells_flashier <- rownames(flashier_snmf_summary$L_pm)
cells_seurat <- seurat_meta$cellID
which(!cells_flashier %in% cells_seurat)
L_pm <- flashier_snmf_summary$L_pm[cells_flashier %in% cells_seurat, ]
F_pm <- flashier_snmf_summary$F_pm

# normalize the F_pm matrix by abs max of each column
denom <- apply(F_pm, 2, function(x) max(abs(x), na.rm = TRUE))
F_pm_norm_col <- sweep(F_pm, 2, denom, "/")

# take the top 50 genes in each column in terms of absolute value
top_genes <- apply(F_pm_norm_col, 2, function(x) {
  top_genes <- names(sort(abs(x), decreasing = TRUE))[1:50]
  return(top_genes)
})

top_50_pos <- apply(F_pm_norm_col,2,function(x) {
  idx <- order(abs(x),decreasing=TRUE)[1:50]
  idx <- idx[x[idx] > 0]
  top_genes_pos <- names(x)[idx]
  return(top_genes_pos)
})

top_50_neg <- apply(F_pm_norm_col,2,function(x) {
  idx <- order(abs(x),decreasing=TRUE)[1:50]
  idx <- idx[x[idx] < 0]
  top_genes_neg <- names(x)[idx]
  return(top_genes_neg)
})

# # make sure the total number of pos and neg genes is 50
# # number of pos genes
# num_pos_genes <- sapply(top_50_pos, length)
# # number of neg genes
# num_neg_genes <- sapply(top_50_neg, length)
# num_pos_genes + num_neg_genes


# top_50_pos, top_50_neg: list of length K (e.g., 200), each element is a character vector of genes



stopifnot(is.list(top_50_pos), is.list(top_50_neg))
stopifnot(length(top_50_pos) == length(top_50_neg))
K <- length(top_50_pos)
sep <- ";"

collapse_genes <- function(x) {
  x <- unique(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) "" else paste(x, collapse = sep)
}

df <- data.frame(
  gp = seq_len(K),
  pos_genes = vapply(top_50_pos, collapse_genes, character(1)),
  neg_genes = vapply(top_50_neg, collapse_genes, character(1)),
  stringsAsFactors = FALSE
)

# save top_pos and top_neg as R objects
saveRDS(top_50_pos, file = "GP_pos.rds")
saveRDS(top_50_neg, file = "GP_neg.rds")

# save csv
write.table(
  df,
  file = "top_genes_by_gp_pos_neg.csv",
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = TRUE,
  qmethod = "double"
)


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(knitr)
  library(kableExtra)
})

infile <- "GP_top_genes.csv"   # adjust path if needed
out_tsv <- "GP_top_genes_pretty.tsv"
out_html <- "GP_top_genes_pretty.html"

# ---- helper to format "A;B;C;..." into multi-line text ----
format_gene_list <- function(x, n_per_line = 6, sep = ";") {
  if (is.null(x) || length(x) == 0 || is.na(x) || trimws(x) == "") return("")
  genes <- unlist(strsplit(x, sep, fixed = TRUE))
  genes <- trimws(genes)
  genes <- genes[genes != ""]

  if (length(genes) == 0) return("")

  # group into chunks of n_per_line and put each chunk on its own line
  idx <- ceiling(seq_along(genes) / n_per_line)
  lines <- tapply(genes, idx, function(g) paste(g, collapse = ", "))
  paste(unname(lines), collapse = "\n")
}

# ---- read & format ----
df <- read_csv(infile, show_col_types = FALSE) %>%
  mutate(
    GP = gp,
    Positive_genes = vapply(pos_genes, format_gene_list, character(1)),
    Negative_genes = vapply(neg_genes, format_gene_list, character(1))
  ) %>%
  select(GP, Positive_genes, Negative_genes) %>%
  arrange(GP)

# ---- write TSV (keeps newlines inside cells; best viewed in text editors / some spreadsheet apps) ----
write_tsv(df, out_tsv)

# ---- write HTML table (recommended for viewing) ----
html_tbl <- df %>%
  knitr::kable(
    format = "html",
    escape = TRUE,
    col.names = c("GP", "Positive genes", "Negative genes")
  ) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover")) %>%
  column_spec(2, width = "45em") %>%
  column_spec(3, width = "45em") %>%
  row_spec(0, bold = TRUE)

save_kable(html_tbl, file = out_html)

message("Done!\n- ", out_tsv, "\n- ", out_html)


