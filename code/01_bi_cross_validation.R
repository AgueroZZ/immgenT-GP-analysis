library(tools)
library(qs)
library(Matrix)
library(fastTopics)
library(flashier)
# data_path <- "data/"
data_path <- "/project2/mstephens/immgent/"

# Read in the meta.information factors_comments_updated.xlsx - factor_pve.csv
factor.info <- read.csv(file = paste0(data_path, "factors_comments_updated.xlsx - factor_pve.csv"), header = TRUE, stringsAsFactors = FALSE)
# sort(factor.info$X[factor.info$pve >= 1e-4])
factors_considered <- sort(factor.info$X[factor.info$strongly.replicated])


# select training cells
flashier_snmf_summary1 <- qs::qread(paste0(data_path, "fit1_summary.qs"))
training_cells <- rownames(flashier_snmf_summary1$L_pm)
shifted_log_counts <- qread(file.path(data_path,"flashier_snmf_matrix.qs"))

# select testing cells
set.seed(1234)
test_indices <- which(!(rownames(shifted_log_counts) %in% training_cells))
second_half <- shifted_log_counts[test_indices,]

# remove genes that are entirely zero
genes_removed_2 <- colSums(second_half) == 0
second_half <- second_half[,!genes_removed_2]

# keep genes that are not removed, and also in the training set
genes_in_training <- rownames(flashier_snmf_summary1$F_pm)
genes_to_keep <- intersect(colnames(second_half), genes_in_training)
second_half <- second_half[, genes_to_keep]

## Model fitting
n2 <- nrow(second_half)
x2 <- rpois(1e7, 1/n2)
s2 <- sd(log(x2 + 1))

# run flashier fixed factors
L_mat_int <- matrix(1, nrow = nrow(second_half), ncol = length(factors_considered))
F_mat <- flashier_snmf_summary1$F_pm[genes_to_keep, factors_considered, drop = FALSE]

flash_fixed_factor <- flash_init(second_half, S = s2, var_type = 2) |>
  flash_set_verbose(1) |>
  flash_factors_init(list(L_mat_int,
                          F_mat), ebnm_fn = c(ebnm_point_exponential, ebnm_point_laplace)) |>
  flash_factors_fix(kset = 1:length(factors_considered),
                    which_dim = "factors")

flash_fixed_factor <- flash_backfit(flash_fixed_factor,extrapolate = FALSE,maxiter = 100,verbose = 2)
flash_fixed_factor <- flash_backfit(flash_fixed_factor,extrapolate = TRUE, maxiter = 40,verbose = 2)

flash_test_B_summary <- list(L_pm = flash_fixed_factor$L_pm,
                             F_pm = flash_fixed_factor$F_pm,
                             pve = flash_fixed_factor$pve,
                             elbo = flash_fixed_factor$elbo,
                             residuals_sd = flash_fixed_factor$residuals_sd)

saveRDS(flash_test_B_summary, file = paste0(data_path, "GEP_flash_selected_summary_testB.rds"))






