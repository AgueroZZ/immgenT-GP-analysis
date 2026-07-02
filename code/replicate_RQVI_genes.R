library(dplyr)
library(stringr)
library(purrr)
library(arrow)
library(tibble)
library(tidyr)
library(readr)

# ---- config ----
dir_gsea <- "~/Desktop/Immgen/immgen-t-factors/data/gsea_RQVI"
pattern  <- "^gsea_seed\\d+_pos\\.parquet$"

read_gsea_one <- function(path) {
  df <- arrow::read_parquet(path) %>% as_tibble()

  stopifnot(all(c("Term", "NES") %in% names(df)))

  df
  # df %>%
  #   group_by(Term) %>%
  #   slice_max(order_by = abs(NES), n = 1, with_ties = FALSE) %>%
  #   ungroup()
}

# ---- find files ----
files <- list.files(dir_gsea, pattern = pattern, full.names = TRUE)
stopifnot(length(files) > 0)

# ---- make nice list names: seed0_pos, seed1_pos, ... ----
list_names <- basename(files) %>%
  str_remove("\\.parquet$") %>%
  str_replace("^gsea_", "")   # "seed0_pos"

# ---- load into a list ----
gsea_list <- set_names(files, list_names) %>%
  map(read_gsea_one)

gsea_long <- imap_dfr(gsea_list, ~{
  m <- str_match(.y, "^seed(\\d+)_(pos|neg)$")
  .x %>%
    mutate(
      seed = as.integer(m[, 2]),
      sign = m[, 3]
    )
})

# keep the result with largest abs NES per Term per Name
gsea_long <- gsea_long %>%
  group_by(Term, Name) %>%
  slice_max(order_by = abs(NES), n = 1, with_ties = FALSE) %>%
  ungroup()

# rename the term just to keep GP number (i.e. GP100_pos to GP100)
gsea_long <- gsea_long %>%
  mutate(Term = str_remove(Term, "_pos$|_neg$")) %>%
  # arrange based on the Term index
  mutate(TermIDX = as.numeric(str_replace(Term, "GP", ""))) %>%
  arrange(TermIDX)


# gsea_long: must contain columns Term, Name, NES
# 1) make a wide table: rows=GPs (Term), cols=Name, values=NES
nes_mat_df <- gsea_long %>%
  mutate(
    Term = as.character(Term),
    Name = as.character(Name),
    gp_id = parse_number(Term),
    name_id = parse_number(Name)
  ) %>%
  group_by(Term, Name, gp_id, name_id) %>%
  slice_max(order_by = abs(NES), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(gp_id, name_id) %>%
  select(Term, Name, NES) %>%
  pivot_wider(names_from = Name, values_from = NES, values_fill = NA_real_) %>%
  arrange(parse_number(Term))

# 2) turn into a numeric matrix with rownames = Term
nes_mat <- nes_mat_df %>%
  column_to_rownames("Term") %>%
  as.matrix()

# optional: reorder columns by numeric order of Name (if Name is numeric-like)
nes_mat <- nes_mat[, order(parse_number(colnames(nes_mat))), drop = FALSE]

abs_nes_mat <- abs(nes_mat)
# produce a heatmap
pheatmap(mat = abs_nes_mat,
         # cluster_rows = F, # row is GP
         # cluster_cols = F, # col is RQVI
         show_rownames = F,
         show_colnames = F,
         color = colorRampPalette(c("white", "red"))(300),
         breaks = seq(1.5, 2.5, length.out = 301),
         main = "Absolute NES values for GPs vs RQVI programs")




###############################################################
###############################################################
###############################################################
# FWER p-val heatmap

gsea_long <- imap_dfr(gsea_list, ~{
  m <- str_match(.y, "^seed(\\d+)_(pos|neg)$")
  .x %>%
    mutate(
      seed = as.integer(m[, 2]),
      sign = m[, 3]
    )
})
# keep the result with smallest pval per Term per Name
gsea_long <- gsea_long %>%
  group_by(Term, Name) %>%
  slice_max(order_by = -(`FWER p-val`), n = 1, with_ties = FALSE) %>%
  ungroup()

# rename the term just to keep GP number (i.e. GP100_pos to GP100)
gsea_long <- gsea_long %>%
  mutate(Term = str_remove(Term, "_pos$|_neg$")) %>%
  # arrange based on the Term index
  mutate(TermIDX = as.numeric(str_replace(Term, "GP", ""))) %>%
  arrange(TermIDX)


pval_mat_df <- gsea_long %>%
  mutate(
    Term = as.character(Term),
    Name = as.character(Name),
    gp_id = parse_number(Term),
    name_id = parse_number(Name)
  ) %>%
  group_by(Term, Name, gp_id, name_id) %>%
  slice_max(order_by = -(`FWER p-val`), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(gp_id, name_id) %>%
  select(Term, Name, `FWER p-val`) %>%
  pivot_wider(names_from = Name, values_from = `FWER p-val`, values_fill = NA_real_) %>%
  arrange(parse_number(Term))

# convert to -log10 p-value matrix (+ 1e-10 for stability)
log10_pval_mat_df <- pval_mat_df %>%
  mutate(Term = as.character(Term)) %>%
  column_to_rownames("Term") %>%
  as.matrix() %>%
  apply(2, function(col) -log10(col + 1e-10))

pheatmap(mat = log10_pval_mat_df,
         cluster_rows = F, # row is GP
         # cluster_cols = F, # col is RQVI
         show_rownames = F,
         show_colnames = F,
         color = colorRampPalette(c("white", "red"))(300),
         breaks = seq(1, 10, length.out = 301),
         main = "-log10 FWER p-values for GPs vs RQVI programs")



###############################################################
###############################################################
###############################################################
# FDR q-val heatmap
gsea_long <- imap_dfr(gsea_list, ~{
  m <- str_match(.y, "^seed(\\d+)_(pos|neg)$")
  .x %>%
    mutate(
      seed = as.integer(m[, 2]),
      sign = m[, 3]
    )
})
# keep the result with smallest qval per Term per Name
gsea_long <- gsea_long %>%
  group_by(Term, Name) %>%
  slice_max(order_by = -(`FDR q-val`), n = 1, with_ties = FALSE) %>%
  ungroup()



# rename the term just to keep GP number (i.e. GP100_pos to GP100)
gsea_long <- gsea_long %>%
  mutate(Term = str_remove(Term, "_pos$|_neg$")) %>%
  # arrange based on the Term index
  mutate(TermIDX = as.numeric(str_replace(Term, "GP", ""))) %>%
  arrange(TermIDX)

qval_mat_df <- gsea_long %>%
  mutate(
    Term = as.character(Term),
    Name = as.character(Name),
    gp_id = parse_number(Term),
    name_id = parse_number(Name)
  ) %>%
  group_by(Term, Name, gp_id, name_id) %>%
  slice_max(order_by = -(`FDR q-val`), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(gp_id, name_id) %>%
  select(Term, Name, `FDR q-val`) %>%
  pivot_wider(names_from = Name, values_from = `FDR q-val`, values_fill = NA_real_) %>%
  arrange(parse_number(Term))

# convert to -log10 q-value matrix (+ 1e-10 for stability)
log10_qval_mat_df <- qval_mat_df %>%
  mutate(Term = as.character(Term)) %>%
  column_to_rownames("Term") %>%
  as.matrix() %>%
  apply(2, function(col) -log10(col + 1e-10))

pheatmap(mat = log10_qval_mat_df,
         cluster_rows = F, # row is GP
         # cluster_cols = F, # col is RQVI
         show_rownames = F,
         show_colnames = F,
         color = colorRampPalette(c("white", "red"))(300),
         breaks = seq(1, 10, length.out = 301),
         main = "-log10 FDR q-values for GPs vs RQVI programs",
         filename = "~/Desktop/Immgen/immgen-t-factors/RQVI_GeneLevel_qval.pdf",
         width = 6, height = 4
         )
