# Verify the published protein-gating panels against the claims made about them.
#
#   Rscript script/verify_gating_gps.R      # exits non-zero on any mismatch
#
# Figure 6e-6j and Figure S6c-S6f are ten protein-gate vs. GP-loading panels,
# drawn by two different scripts from one curated pool (well_aligned_gps). Four
# things are asserted in prose -- in the two captions, in Extended Data Table 7's
# caption, and in code/R/citeseq_shared_setup.R's comment -- and nothing in the
# build checks any of them:
#
#   1. the ten GPs come from the curated well_aligned_gps pool,
#   2. the two figures show disjoint GPs (no GP is published twice),
#   3. every panel letter a script declares has a PDF on disk, and no orphan
#      panel PDF is left over from an earlier lettering,
#   4. every marker in those ten signatures is actually gated on -- it has a
#      column in the ADT matrix AND a manually reviewed threshold -- so no panel
#      subtitle can advertise a marker the gate silently skipped.
#
# The GP-to-letter maps are read out of the figure scripts rather than retyped,
# so an edit there cannot drift away from this check.

library(Matrix) # protein matrices are dgCMatrix; must be attached for `[` to dispatch
suppressPackageStartupMessages(library(dplyr)) # the setup uses %>%

data_path <- "data/"
source("code/R/citeseq_shared_setup.R")

failures <- character(0)
check <- function(ok, msg) if (!isTRUE(ok)) failures <<- c(failures, msg)

# Read a named GP -> letter vector out of a script by the name it is assigned to.
read_letter_map <- function(file, varname) {
  lines <- readLines(file, warn = FALSE)
  hit <- grep(paste0("^", varname, " <- "), lines, value = TRUE)
  if (length(hit) != 1) {
    stop(sprintf("%s is not assigned exactly once in %s", varname, file))
  }
  eval(parse(text = sub(paste0("^", varname, " <- "), "", hit)))
}

fig6 <- read_letter_map("script/Figure6.R", "fig6_gating")
figs6 <- read_letter_map("script/FigureS6.R", "figs6_gating")
gps <- c(names(fig6), names(figs6))

cat("=== 0. the panels each script declares ===\n")
cat("Figure 6 :", paste(sprintf("%s = %s", fig6, names(fig6)), collapse = ", "), "\n")
cat("Figure S6:", paste(sprintf("%s = %s", figs6, names(figs6)), collapse = ", "), "\n")

cat("\n=== 1. every gated GP comes from the curated well_aligned_gps pool ===\n")
outside <- setdiff(gps, well_aligned_gps)
cat(sprintf("pool: %d GPs | published: %d | outside the pool: %s\n",
            length(well_aligned_gps), length(gps),
            if (length(outside)) paste(outside, collapse = ", ") else "none"))
check(length(outside) == 0,
      sprintf("gated GP(s) not in well_aligned_gps: %s", paste(outside, collapse = ", ")))

cat("\n=== 2. Figure 6 and Figure S6 show disjoint GPs ===\n")
both <- intersect(names(fig6), names(figs6))
cat("in both figures:", if (length(both)) paste(both, collapse = ", ") else "none", "\n")
check(length(both) == 0, sprintf("GP(s) published in both figures: %s", paste(both, collapse = ", ")))
check(!anyDuplicated(gps), "a GP is lettered twice within one figure")

cat("\n=== 3. every declared letter has a panel PDF, and there are no orphans ===\n")
for (spec in list(list(map = fig6, dir = "figures/final-selected/Figure 6/", pat = "^6[e-z]\\.pdf$"),
                  list(map = figs6, dir = "figures/final-selected/Figure S6/", pat = "^s6[c-z]\\.pdf$"))) {
  want <- paste0(spec$map, ".pdf")
  missing <- want[!file.exists(file.path(spec$dir, want))]
  on_disk <- list.files(spec$dir, pattern = spec$pat)
  orphan <- setdiff(on_disk, want)
  cat(sprintf("%-34s declared %d | missing: %-12s | orphan: %s\n", spec$dir, length(want),
              if (length(missing)) paste(missing, collapse = ",") else "none",
              if (length(orphan)) paste(orphan, collapse = ",") else "none"))
  check(length(missing) == 0, sprintf("%s: no PDF for %s", spec$dir, paste(missing, collapse = ", ")))
  check(length(orphan) == 0,
        sprintf("%s: %s is not produced by any current panel -- stale lettering?",
                spec$dir, paste(orphan, collapse = ", ")))
}

cat("\n=== 4. no marker in these ten signatures is silently skipped ===\n")
# Mirrors plot_gated_gp_vs_protein()'s gate exactly: a marker is applied only if
# it survives select_proteins, has an ADT column, and has a manual threshold.
adt_cols <- colnames(protein_mat_normalized_lognorm)
have_threshold <- threshold_results_subset_manual$Protein
dropped_any <- character(0)
for (gp in gps) {
  requested <- unlist(strsplit(c(df_markers2[gp, "Positive"], df_markers2[gp, "Negative"]), ", "))
  requested <- intersect(requested[nzchar(requested)], select_proteins)
  no_column <- setdiff(requested, adt_cols)
  no_thresh <- setdiff(intersect(requested, adt_cols), have_threshold)
  dropped <- c(no_column, no_thresh)
  dropped_any <- c(dropped_any, dropped)
  cat(sprintf("%-6s %2d marker(s) requested | dropped: %s\n", gp, length(requested),
              if (length(dropped)) paste(dropped, collapse = ", ") else "none"))
  check(length(dropped) == 0,
        sprintf("%s gates on %s, which is skipped (no ADT column or no manual threshold) -- the caption must say so",
                gp, paste(dropped, collapse = ", ")))
}
cat("markers dropped across all ten panels:",
    if (length(dropped_any)) paste(unique(dropped_any), collapse = ", ") else "none", "\n")

cat("\n==================== RESULT ====================\n")
if (length(failures)) {
  cat("MISMATCH:\n", paste0("  - ", failures, collapse = "\n"), "\n", sep = "")
  quit(status = 1)
}
cat(sprintf("PASS: %d gating panels, all from well_aligned_gps, disjoint across the two\n", length(gps)))
cat("      figures, all present on disk, and gating on every marker they name.\n")
