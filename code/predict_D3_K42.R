is_D3 <- as.numeric(seurat_meta$condition_detailed == "D3")

roc_obj <- pROC::roc(
  response = is_D3,
  predictor = L_pm[,42],
  levels = c(0, 1),
  ci = TRUE,
  plot = TRUE
)

# plot ROC curve
plot(roc_obj, main = "ROC Curve for L_pm[,42] predicting D3 condition",
     xlim = c(1, 0), ylim = c(0, 1),
     col = "blue", lwd = 2)
