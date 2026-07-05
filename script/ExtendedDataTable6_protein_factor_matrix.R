# Extended Data Table 6: protein factor matrix.
#
# The final processed CITE-seq protein factor matrix used by Figure 6 (panel 6b
# heatmap), exported as a plain table: one row per GP (GP1..GP200), one column
# per protein, holding the protein's factor score in that GP. Built exactly as
# script/Figure6.R (lines 58-65) does from the projected protein EBMF fit --
# drop isotype / THY1.1 / excluded / non-"good" proteins, then scale each GP
# column so its maximum absolute score is 1 (the same max|.|=1 normalization the
# heatmap uses). See code/pipeline/04_protein_projection.R for how the raw
# protein factor matrix is fit, and code/R/citeseq_shared_setup.R:35-40 for the
# protein-selection filters reproduced here.
#
# The r-object orientation is proteins x GPs; we transpose to GP x protein so GP
# is the row, as requested for the table.

data_path <- "data/"
output_path <- "figures/generated/"

Protein_flash_result <- readRDS(paste0(data_path, "protein_flash_selected_summary_lognorm_backfit200.rds"))
Protein_F_pm_raw <- Protein_flash_result$F_pm

# Protein-selection filters (mirror code/R/citeseq_shared_setup.R:35-40)
isotype_proteins <- grep("^Isotype", rownames(Protein_F_pm_raw), value = TRUE)
proteins_quality <- read.csv(paste0(data_path, "TableS4_citeseq_qc_20250513.csv"), header = TRUE, stringsAsFactors = FALSE, skip = 1)
good_proteins <- c(proteins_quality$protein[proteins_quality$classification == "good"], "IL2RA.CD25", "ITB7", "CD69")
exclude_proteins <- c("CD19", "CD34", "CD45.1", "CD45.2", "CD138", "TCRVA2", "TER119")
thy11_proteins <- grep("THY1.1", rownames(Protein_F_pm_raw), value = TRUE)

# Filter + per-GP-column max|.|=1 scaling (mirror script/Figure6.R:58-65)
Protein_F_pm <- Protein_F_pm_raw[!rownames(Protein_F_pm_raw) %in% isotype_proteins, ]
Protein_F_pm <- Protein_F_pm[rownames(Protein_F_pm) %in% good_proteins, ]
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% exclude_proteins, ]
Protein_F_pm <- Protein_F_pm[!rownames(Protein_F_pm) %in% thy11_proteins, ]
D_lognorm <- diag(1 / apply(Protein_F_pm, 2, function(x) max(abs(x))))
Protein_F_pm <- Protein_F_pm %*% D_lognorm
colnames(Protein_F_pm) <- paste0("GP", seq_len(ncol(Protein_F_pm)))
Protein_F_pm[is.na(Protein_F_pm)] <- 0

# Transpose to GP (row) x protein (column) and write out.
protein_factor_gp <- as.data.frame(t(Protein_F_pm), check.names = FALSE, stringsAsFactors = FALSE)
protein_factor_gp <- data.frame(GP = rownames(protein_factor_gp), protein_factor_gp,
                                check.names = FALSE, stringsAsFactors = FALSE)
rownames(protein_factor_gp) <- NULL

write.csv(
  protein_factor_gp,
  file = paste0(output_path, "ExtendedDataTable6_protein_factor_matrix.csv"),
  row.names = FALSE
)
