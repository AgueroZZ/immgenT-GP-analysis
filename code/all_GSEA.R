### Run GSEA on F matrix
library(qs)
library(Matrix)
library(flashier)
library(parallel)
library(singlecelljamboreeR)
library(pathways)
library(openxlsx)
library(ZemmourLib)
library(dplyr)
library(tidyr)

# data_path <- "/project2/mstephens/immgent"
data_path <- "data/"
code_path <- "code/"
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
F_pm <- flashier_snmf_summary$F_pm
colnames(F_pm) <- paste0("K", 1:ncol(F_pm))

# Load gene sets
data(gene_sets_mouse)

# Subset to MSigDB-C7 gene sets
gene_sets     <- gene_sets_mouse$gene_sets
gene_info     <- gene_sets_mouse$gene_info
gene_set_info <- gene_sets_mouse$gene_set_info
j <- which(with(gene_set_info,
                  database == "MSigDB-C7"))
genes <- sort(intersect(rownames(F_pm),gene_info$Symbol))
i1 <- match(genes,gene_info$Symbol)
i2 <- match(genes,rownames(F_pm))
gene_info     <- gene_info[i1,]
gene_set_info <- gene_set_info[j,]
gene_sets     <- gene_sets[i1,j]
F_pm_considered             <- F_pm[i2,]
rownames(gene_sets)     <- gene_info$Symbol
rownames(gene_set_info) <- gene_set_info$id

# factors_to_analyze <- paste0("K", c(58,68))
# gsea <- pathways::perform_gsea(
#   Z = F_pm[,factors_to_analyze],
#   gene_sets = gene_sets,
#   verbose = FALSE)
# sum(qvalue::qvalue(gsea$pval[,1])$qvalues < 0.05)
# sum(qvalue::qvalue(gsea$pval[,2])$qvalues < 0.05)

# run result on all factors
factors_to_analyze <- colnames(F_pm)
gsea <- pathways::perform_gsea(
  Z = F_pm[,factors_to_analyze],
  gene_sets = gene_sets,
  verbose = FALSE)

NES  <- gsea$NES
pval <- gsea$pval

qval <- pval
for (j in seq_len(ncol(pval))) {
  pj <- pval[, j]
  ok <- !is.na(pj)
  qj <- rep(NA_real_, length(pj))
  if (sum(ok) > 0) {
    qj[ok] <- qvalue::qvalue(pj[ok])$qvalues
  }
  qval[, j] <- qj
}

wb <- createWorkbook()
addWorksheet(wb, "NES")
writeData(wb, "NES", data.frame(gene_set = rownames(NES), NES, check.names = FALSE))
addWorksheet(wb, "pval")
writeData(wb, "pval", data.frame(gene_set = rownames(pval), pval, check.names = FALSE))
addWorksheet(wb, "qval")
writeData(wb, "qval", data.frame(gene_set = rownames(qval), qval, check.names = FALSE))
addWorksheet(wb, "gene_set_info")
writeData(wb, "gene_set_info", gene_set_info)
saveWorkbook(wb, file.path(data_path, "gsea_results.xlsx"), overwrite = TRUE)


# edit gsea_results.xlsx
gsea_NES_results <- read_excel("data/gsea_results.xlsx",
                               sheet = "NES")
gsea_qval_results <- read_excel("data/gsea_results.xlsx",
                                sheet = "qval")
gene_set_info <- read_excel("data/gsea_results.xlsx",
                            sheet = "gene_set_info")

# add gene-set-name to each sheet then store gsea_results_named.xlsx
gsea_NES_results_named <- gsea_NES_results %>%
  left_join(gene_set_info %>% select(id, name), by = c("gene_set" = "id")) %>%
  select(gene_set, name, everything())
gsea_qval_results_named <- gsea_qval_results %>%
  left_join(gene_set_info %>% select(id, name), by = c("gene_set" = "id")) %>%
  select(gene_set, name, everything())

wb <- createWorkbook()
addWorksheet(wb, "NES")
writeData(wb, "NES", gsea_NES_results_named)
addWorksheet(wb, "qval")
writeData(wb, "qval", gsea_qval_results_named)
addWorksheet(wb, "gene_set_info")
writeData(wb, "gene_set_info", gene_set_info)
saveWorkbook(wb, file.path(data_path, "gsea_results_named.xlsx"), overwrite = TRUE)


#####################################
#####################################
#####################################
#####################################
#####################################
## Compare active CD4 versus CD8
seurat_meta <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
L_pm <- flashier_snmf_summary$L_pm[match(seurat_meta$cellID, rownames(flashier_snmf_summary$L_pm)),]
colnames(L_pm) <- paste0("K", 1:ncol(L_pm))

# filter out and normalize
filter_cells_by_total_membership <- function (L, max_val = 10, numiter = 10) {
  n <- nrow(L)
  rows <- 1:n
  for (iter in 1:numiter) {
    x <- rowSums(L)
    cat(sprintf("%d. Filtered out %d cells.\n",iter,sum(x > max_val)))
    i <- which(x <= max_val)
    L <- L[i,]
    rows <- rows[i]
    d <- apply(L,2,max)
    L <- scale_cols(L,1/d)
  }
  return(rows)
}
scale_cols <- function (A, b) t(t(A) * b)
colnames(L_pm) <- paste0("K", seq_len(ncol(L_pm)))
colnames(F_pm) <- paste0("F", seq_len(ncol(F_pm)))
D <- diag(1 / apply(L_pm, 2, function(x) max(x)))
L <- L_pm %*% D
cells <- filter_cells_by_total_membership(L,numiter = 12)
seurat_meta_filtered <- seurat_meta[cells,]
L_pm_filtered <- L_pm[cells,]
d <- apply(L_pm_filtered,2,max)
L_pm_filtered <- scale_cols(L_pm_filtered,1/d)
F_pm_filtered <- scale_cols(F_pm,d)
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered),]

activated_cd4 <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level2_group == "activated" & seurat_meta_filtered$annotation_level1 == "CD4"]
activated_cd8 <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level2_group == "activated" & seurat_meta_filtered$annotation_level1 == "CD8"]
all_cd4 <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 == "CD4"]
all_cd8 <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 == "CD8"]

# take all cells within L_pm_filtered
FlashierDGE_corrected <- function (F1, L1, group1, group2, title_plot = ""){
  loadings_group1 = colMeans(L1[group1, ])
  loadings_group2 = colMeans(L1[group2, ])
  loadings_groups = colMeans(L1[c(group1, group2), ])
  mean_genes_group1 = F1 %*% loadings_group1
  mean_genes_group2 = F1 %*% loadings_group2
  mean_genes = F1 %*% colMeans(L1[c(group1, group2), ])
  # fc_loadings = loadings_group1 - loadings_group2
  mean_change_loadings = loadings_group1 - loadings_group2
  # fc_genes = F1 %*% fc_loadings %>% as.data.frame()
  fc_genes = F1 %*% mean_change_loadings %>% as.data.frame()
  fc_loadings = log2((loadings_group1 + 1e-10)/(loadings_group2 + 1e-10))

  # compute median
  median_loadings_group1 = apply(L1[group1, ], 2, median)
  median_loadings_group2 = apply(L1[group2, ], 2, median)
  median_change_loadings = median_loadings_group1 - median_loadings_group2
  log2FC_median_loadings = log2((median_loadings_group1 + 1e-10)/(median_loadings_group2 + 1e-10))

  vplot = data.frame(SYMBOL = names(mean_change_loadings),
                     mean_change = mean_change_loadings,
                     median_change = median_change_loadings,
                     log2FC = fc_loadings,
                     log2FC_median = log2FC_median_loadings,
                     AveExpr = loadings_groups)
  max_mc = ceiling(max(abs(vplot$mean_change)))
  top_genes = vplot %>% dplyr::arrange(dplyr::desc(abs(mean_change))) %>%
    utils::head(50)
  p1 = ggplot2::ggplot(data = vplot) + ggplot2::geom_point(ggplot2::aes(x = mean_change,
                                                                        y = AveExpr), colour = "black", alpha = I(1), size = I(1)) +
    ggplot2::xlim(-max_mc, max_mc) + ggrepel::geom_text_repel(data = top_genes,
                                                              ggplot2::aes(x = mean_change, y = AveExpr, label = SYMBOL),
                                                              size = 3, color = "red", box.padding = 0.35, point.padding = 0.5,
                                                              segment.color = "grey50", max.overlaps = 20) + ggplot2::theme_minimal() +
    ggplot2::labs(x = "Difference in Mean Loading", y = "Average Loading",
                  title = title_plot)

  vplot_genes = data.frame(SYMBOL = rownames(fc_genes), log2FC = fc_genes[,
                                                                          1]/log(2), AveExpr = mean_genes[, 1])
  max_fc_genes = ceiling(max(abs(vplot_genes$log2FC)))
  top_genes_genes <- vplot_genes %>% dplyr::arrange(dplyr::desc(abs(log2FC))) %>%
    utils::head(50)
  p2 <- ggplot2::ggplot(data = vplot_genes) + scattermore::geom_scattermore(ggplot2::aes(x = log2FC,
                                                                                         y = AveExpr), colour = "black", alpha = I(1), size = I(1),
                                                                            pixels = c(512, 512)) + ggplot2::xlim(-max_fc_genes,
                                                                                                                  max_fc_genes) + ggrepel::geom_text_repel(data = top_genes_genes,
                                                                                                                                                           ggplot2::aes(x = log2FC, y = AveExpr, label = SYMBOL),
                                                                                                                                                           size = 3, color = "red", box.padding = 0.35, point.padding = 0.5,
                                                                                                                                                           segment.color = "grey50", max.overlaps = 20) + ggplot2::theme_minimal() +
    ggplot2::labs(x = "Fold Change (log2)", y = "Average Expression",
                  title = title_plot)
  list(p1 = p1, p2 = p2, diff_factors = vplot, diff_genes = vplot_genes)
}

# construct the list of LDF
Dl <- apply(L_pm_filtered,2,max)
Df <- apply(abs(F_pm_filtered),2,max)
LDF_list <- list()
LDF_list$L_pm_filtered_norm <- scale_cols(L_pm_filtered,1/Dl)
LDF_list$F_pm_filtered_norm <- scale_cols(F_pm_filtered,1/Df)
LDF_list$D <- Dl * Df
saveRDS(LDF_list, file.path(data_path, "LDF_filtered_list.rds"))

# plot the relative importance of each GP
plot(LDF_list$D, xlab = "GP", ylab = "LDF", main = "Relative importance of each GP")
lines(smooth.spline(LDF_list$D), col = "red")

plot(flashier_snmf_summary$pve ~ LDF_list$D, log = "xy")
lm_fit <- lm(log10(flashier_snmf_summary$pve) ~ log10(LDF_list$D))
abline(lm_fit, col = "blue")
text(x = min(LDF_list$D),
     col = "blue",
     y = max(flashier_snmf_summary$pve),
     labels = paste0("R2 = ", round(summary(lm_fit)$r.squared, 3)), pos = 4)

# DGE
DGE_result <- FlashierDGE_corrected(F1 = F_pm_filtered,
                          L1 = L_pm_filtered,
                          group1 = activated_cd4,
                          group2 = activated_cd8,
                          title_plot = "Activated CD4 vs CD8")

# top five up regulated programs in CD8 relative to CD4
top_cd8_programs <- DGE_result$diff_factors %>%
  dplyr::arrange(dplyr::desc(-mean_change)) %>%
  utils::head(5) %>%
  pull(SYMBOL)

# top five down regulated programs in CD8 relative to CD4
top_cd4_programs <- DGE_result$diff_factors %>%
  dplyr::arrange(dplyr::desc(mean_change)) %>%
  utils::head(5) %>%
  pull(SYMBOL)

# pull all relevant GPs out
gp_of_interest <- c(top_cd8_programs, top_cd4_programs)
gsea_results_of_interest <- qval[,gp_of_interest, drop = FALSE]



# Greedy, non-overlapping selection:
# For each gp in order, take top_n smallest qvals among remaining gene sets,
# then remove them from the pool.
keep_top_gene_sets_per_gp_greedy <- function(qmat, gp_order, top_n = 5) {
  stopifnot(!is.null(rownames(qmat)), !is.null(colnames(qmat)))
  stopifnot(all(gp_order %in% colnames(qmat)))

  remaining <- rownames(qmat)
  picked <- vector("list", length(gp_order))
  names(picked) <- gp_order

  for (gp in gp_order) {
    v <- qmat[remaining, gp]
    ok <- !is.na(v)
    if (sum(ok) < top_n) {
      stop(paste0("Not enough remaining non-NA gene sets for ", gp,
                  ": need ", top_n, ", have ", sum(ok), "."))
    }
    # take top_n smallest qvals
    ord <- order(v[ok], na.last = NA)
    chosen <- remaining[ok][ord][seq_len(top_n)]

    picked[[gp]] <- chosen
    remaining <- setdiff(remaining, chosen)
  }

  # Build row ids to guarantee 50 rows even if gene_set names repeat (they won't here)
  out_ids <- unlist(mapply(
    function(gp, gs) paste0(gs, "__", gp, "_r", seq_along(gs)),
    gp_order, picked, SIMPLIFY = FALSE
  ))

  out_gene_sets <- unlist(picked, use.names = FALSE)

  list(row_id = out_ids, gene_set = out_gene_sets)
}

sel <- keep_top_gene_sets_per_gp_greedy(
  qmat = gsea_results_of_interest,
  gp_order = gp_of_interest,
  top_n = 5
)

top5_mat <- gsea_results_of_interest[sel$gene_set, gp_of_interest, drop = FALSE]
rownames(top5_mat) <- sel$row_id




