setwd("/Users/ziangzhang/Desktop/Immgen/immgenT-GP-analysis")

library(ggplot2)
library(dplyr)

data_path  <- "data/"
fig1_path  <- "figures/generated/Figure 1/"
figs1_path <- "figures/generated/Figure S1/"
png1_path  <- "analysis/assets/Figure1/"
pngs1_path <- "analysis/assets/FigureS1/"

# ── shared data ──────────────────────────────────────────────────────────────
L_pm_filtered <- readRDS(paste0(data_path, "L_pm_filtered.rds"))
F_pm_filtered <- readRDS(paste0(data_path, "F_pm_filtered.rds"))
seurat_meta   <- readRDS(paste0(data_path, "igt1_96_withtotalvi20260206_clean_ADTonly.Rds"))@meta.data
seurat_meta_filtered <- seurat_meta[rownames(L_pm_filtered), ]

# drop thymocytes (matching Figure1.R)
non_thymo <- seurat_meta_filtered$cellID[seurat_meta_filtered$annotation_level1 != "thymocyte"]
L_pm_filtered_nt <- L_pm_filtered[non_thymo, ]
F_pm_filtered_nt <- F_pm_filtered  # F is gene-level, no thymocyte filter

# ── normalise ─────────────────────────────────────────────────────────────────
L_pm_norm_col <- L_pm_filtered_nt / matrix(
  apply(L_pm_filtered_nt, 2, max), nrow = nrow(L_pm_filtered_nt),
  ncol = ncol(L_pm_filtered_nt), byrow = TRUE
)
gp_active_cell_prop <- colSums(L_pm_norm_col > 1e-1) / nrow(L_pm_norm_col)

F_pm_norm_col <- F_pm_filtered / matrix(
  apply(F_pm_filtered, 2, function(x) max(abs(x))), nrow = nrow(F_pm_filtered),
  ncol = ncol(F_pm_filtered), byrow = TRUE
)
gp_active_gene_counts <- colSums(abs(F_pm_norm_col) > 0.25)

# ── 1E: log-scale histogram of active-cell proportions ───────────────────────
p_1E <- ggplot(data.frame(prop = gp_active_cell_prop), aes(x = prop)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  scale_x_log10(labels = scales::label_percent()) +
  annotation_logticks(sides = "b") +
  labs(x = "Proportion of highly active cells per GP (log scale)", y = "Count",
       title = "Histogram of highly active cells per GP (proportion)") +
  theme_minimal(base_size = 13)
ggsave(paste0(fig1_path,  "1E.pdf"), plot = p_1E, width = 6, height = 4, dpi = 300)
ggsave(paste0(png1_path,  "1E.png"), plot = p_1E, width = 6, height = 4, dpi = 150)

# ── 1F: log-scale histogram of active-gene counts ────────────────────────────
p_1F <- ggplot(data.frame(count = gp_active_gene_counts), aes(x = count)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "white") +
  scale_x_log10(labels = scales::label_comma()) +
  annotation_logticks(sides = "b") +
  labs(x = "Number of highly active genes per GP (log scale)", y = "Count",
       title = "Histogram of highly active genes per GP") +
  theme_minimal(base_size = 13)
ggsave(paste0(fig1_path, "1F.pdf"), plot = p_1F, width = 6, height = 4, dpi = 300)
ggsave(paste0(png1_path, "1F.png"), plot = p_1F, width = 6, height = 4, dpi = 150)

# ── S1E: EBMF sparsity scatter (expected active-gene prop x vs cell prop y) ──
# Uses the point-mass prior weights learned by EBMF (p = 1 - pi0).
# Same quantities as Figure 1's scatter_e_active_cells_vs_genes panel.
load(paste0(data_path, "flashier_snmf_fitted_prior.rda"))
n_genes  <- nrow(F_pm_filtered)   # total genes in the factorization

n_gp     <- length(flashier_snmf_fitted_prior$L_ghat)
l_pi_vec <- sapply(seq_len(n_gp), function(i) flashier_snmf_fitted_prior$L_ghat[[i]]$pi[1])
f_pi_vec <- sapply(seq_len(n_gp), function(i) flashier_snmf_fitted_prior$F_ghat[[i]]$pi[1])
p_cells     <- 1 - l_pi_vec
n_genes_act <- (1 - f_pi_vec) * n_genes

scatter_df <- data.frame(n_genes = n_genes_act, prop_cells = p_cells)

pct_breaks <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 1)
p_S1E <- ggplot(scatter_df, aes(x = n_genes, y = prop_cells)) +
  geom_point(size = 2, alpha = 0.7, color = "steelblue") +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(breaks = pct_breaks, labels = function(x) paste0(x * 100, "%")) +
  annotation_logticks(sides = "bl") +
  labs(
    x = "Expected number of active genes per GP (log scale)",
    y = "Expected proportion of active cells per GP (log scale)",
    title = "Expected active genes vs. active-cell proportion per GP\n(EBMF sparsity priors)"
  ) +
  theme_minimal(base_size = 13)
ggsave(paste0(figs1_path, "S1E.pdf"), plot = p_S1E, width = 6, height = 5, dpi = 300)
ggsave(paste0(pngs1_path, "S1E.png"), plot = p_S1E, width = 6, height = 5, dpi = 150)

# ── sync PNGs into docs/assets so the built site picks them up ───────────────
# (wflow_build() on individual Rmds does NOT re-run workflowr's site-wide
# static-asset copy, so docs/assets/ can otherwise go stale.)
dir.create("docs/assets/Figure1",  recursive = TRUE, showWarnings = FALSE)
dir.create("docs/assets/FigureS1", recursive = TRUE, showWarnings = FALSE)
file.copy(paste0(png1_path,  "1E.png"),  "docs/assets/Figure1/1E.png",   overwrite = TRUE)
file.copy(paste0(png1_path,  "1F.png"),  "docs/assets/Figure1/1F.png",   overwrite = TRUE)
file.copy(paste0(pngs1_path, "S1E.png"), "docs/assets/FigureS1/S1E.png", overwrite = TRUE)

message("Done: 1E.pdf/png, 1F.pdf/png, S1E.pdf/png regenerated and synced to docs/assets.")
