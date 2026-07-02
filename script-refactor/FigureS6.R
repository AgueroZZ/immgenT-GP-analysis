# Figure S6. GP loadings recover protein-gated populations across GPs.
#
# Extension of Figure 6c-f to all well-aligned GPs (see
# figures/final-selected/bits/Figure S6/FigureS6_caption.md), shown as two
# gallery pages (s6-1 and s6-2). For each GP, cells are shown twice on the
# same MDE embedding: left, cells passing the GP's curated protein gate;
# right, an equally sized set of cells with the highest GP loading.
#
# Source: ported from script/gated_protein_loading_plot.R's live gallery
# section (the ~450 preceding lines of commented-out single-GP exploratory
# calls are dropped -- they never produced a saved output).

library(ggplot2)
library(dplyr)
library(patchwork)
library(Matrix)

data_path <- "data/"
figure_path <- "figure-refactor/Figure S6/"
source("code-refactor/R/gated_protein_helpers.R")
source("code-refactor/R/citeseq_shared_setup.R")

# GPs shown in main Figure 6 (panels 6c-6f) are excluded from this gallery.
GPs_fig6 <- c("GP171", "GP23", "GP12", "GP80")
other_GPs <- setdiff(well_aligned_gps, GPs_fig6)

# A few GPs get slightly larger highlighted points for visibility even at
# high cell counts.
enlarge_gps <- c("GP8", "GP30", "GP170", "GP107")

plots_all <- lapply(other_GPs, function(gp) {
  k_name <- paste0("K", sub("^GP", "", gp))
  plot_gated_gp_vs_protein(
    gp_name = k_name,
    df_markers = df_markers2,
    protein_mat = protein_mat_normalized_lognorm,
    loading_mat = L_pm_for_gating,
    mde_emb = mde_result,
    missing_threshold_action = "skip",
    threshold_df = threshold_results_subset_manual,
    exclude_cells = c(thymocyte_cells, proliferating_cells, miniverse_cells),
    selected_proteins = select_proteins,
    loading_q = NULL,
    min_pointsize = if (gp %in% enlarge_gps) 3L else 0L
  )
})

v_bar <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "grey60", color = NA))

# Sorted gallery: GPs in numerical order, one PDF per page (2 pages here),
# 2 GP-units per row x 6 rows = 12 GPs per page; panels labeled a, b, c...
sort_idx <- order(as.integer(sub("^GP", "", other_GPs)))
plots_sorted <- plots_all[sort_idx]

n_pages_sorted <- ceiling(length(plots_sorted) / 12)
for (page_idx in seq_len(n_pages_sorted)) {
  idx_start <- (page_idx - 1) * 12 + 1
  idx_end <- min(page_idx * 12, length(plots_sorted))
  page_plots <- plots_sorted[idx_start:idx_end]
  n_on_page <- length(page_plots)

  labeled_page_plots <- lapply(seq_len(n_on_page), function(i) {
    unit <- page_plots[[i]]
    p1_labeled <- unit[[1]] + labs(tag = letters[i]) + theme(plot.tag = element_text(size = 14, face = "bold"))
    p1_labeled + unit[[2]] + plot_layout(ncol = 2)
  })

  rows_list <- lapply(1:6, function(r) {
    i_left <- (r - 1) * 2 + 1
    i_right <- (r - 1) * 2 + 2
    gp_left <- if (i_left <= n_on_page) labeled_page_plots[[i_left]] else plot_spacer()
    gp_right <- if (i_right <= n_on_page) labeled_page_plots[[i_right]] else plot_spacer()
    (gp_left | v_bar | gp_right) + plot_layout(widths = c(1, 0.03, 1))
  })
  combined <- wrap_plots(rows_list, ncol = 1)

  # cairo_pdf() does not reliably truncate/overwrite an existing file of a
  # different size in place -- remove any stale output first.
  out_page <- paste0(figure_path, sprintf("s6-%d.pdf", page_idx))
  if (file.exists(out_page)) unlink(out_page)

  graphics.off()
  cairo_pdf(out_page, width = 13, height = 18)
  showtext::showtext_begin()
  print(combined)
  showtext::showtext_end()
  dev.off()
}
