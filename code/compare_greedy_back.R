flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))
flashier_snmf_summary1 <- qs::qread(paste0(data_path, "fit1_summary_noback.qs"))
flashier_snmf_summary2 <- qs::qread(paste0(data_path, "fit2_summary_noback.qs"))

L_all <- flashier_snmf_summary$L_pm
L_1 <- flashier_snmf_summary1$L_pm
L_2 <- flashier_snmf_summary2$L_pm

D_all <- diag(1 / apply(L_all, 2, function(x) max(x)))
L_all <- L_all %*% D_all
seurat_meta_all <- seurat_meta[match(rownames(L_all), seurat_meta$cellID),]

D_1 <- diag(1 / apply(L_1, 2, function(x) max(x)))
L_1 <- L_1 %*% D_1
seurat_meta_1 <- seurat_meta[match(rownames(L_1), seurat_meta$cellID),]

D_2 <- diag(1 / apply(L_2, 2, function(x) max(x)))
L_2 <- L_2 %*% D_2
seurat_meta_2 <- seurat_meta[match(rownames(L_2), seurat_meta$cellID),]

library(ggplot2)

factor_id <- 2

df <- data.frame(
  value = c(L_all[, factor_id], L_1[, factor_id], L_2[, factor_id]),
  group = factor(
    rep(c("All", "First Half", "Second Half"),
        times = c(nrow(L_all), nrow(L_1), nrow(L_2)))
  )
)

ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outliers = TRUE, alpha = 0.7) +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb")) +
  labs(y = "Loading value", x = NULL, title = paste("Factor", factor_id, "Comparison")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")



