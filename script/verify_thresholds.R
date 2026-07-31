# Verify that the PUBLISHED protein thresholds are exactly the ones the figures
# gate on.
#
#   Rscript script/verify_thresholds.R      # exits non-zero on any mismatch
#
# The manual positivity thresholds are hand-curated in
# data/Thresholds_Selected_Proteins.csv (see code/pipeline/03_protein_thresholds.R
# for why that file is an input, not an output). Three things therefore have to
# agree, and nothing in the build enforces it:
#
#   the curated CSV  ->  threshold_results_subset_manual   (code/R/citeseq_shared_setup.R)
#                    ->  the threshold_df that Figure6.R and FigureS6.R gate on
#                    ->  Extended Data Table 7, which publishes the values
#
# This script builds the middle one by sourcing the real setup -- not a copy of
# its logic -- and diffs it against the published table. It also reports the
# reverse gap: markers that are gated on but have no threshold, and so are
# silently skipped (missing_threshold_action = "skip" in the gating helpers).

library(Matrix) # protein matrices are dgCMatrix; must be attached for `[` to dispatch
suppressPackageStartupMessages(library(dplyr)) # the setup uses %>%

data_path <- "data/"
source("code/R/citeseq_shared_setup.R")

runtime <- threshold_results_subset_manual # what the gates actually receive
published <- read.csv(
  "figures/final-selected/ExtendedDataTable7_protein_thresholds.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)
names(published) <- c("Protein", "Threshold")

failures <- character(0)
check <- function(ok, msg) if (!isTRUE(ok)) failures <<- c(failures, msg)

cat("=== 1. protein set ===\n")
cat(sprintf("runtime: %d   published: %d\n", nrow(runtime), nrow(published)))
cat("in runtime but not published:", paste(setdiff(runtime$Protein, published$Protein), collapse = ", "), "\n")
cat("published but not in runtime:", paste(setdiff(published$Protein, runtime$Protein), collapse = ", "), "\n")
check(setequal(runtime$Protein, published$Protein), "protein sets differ")
check(!any(duplicated(runtime$Protein)), "duplicate protein in the runtime table")
check(!any(duplicated(published$Protein)), "duplicate protein in the published table")

cat("\n=== 2. values, protein by protein ===\n")
both <- merge(runtime, published, by = "Protein", suffixes = c("_runtime", "_published"), all = TRUE)
both$diff <- both$Threshold_runtime - both$Threshold_published
bad <- both[!is.finite(both$diff) | both$diff != 0, ]
if (nrow(bad)) print(bad, row.names = FALSE) else cat(sprintf("all %d values identical (diff == 0)\n", nrow(both)))
check(nrow(bad) == 0, "threshold values differ between runtime and published table")

cat("\n=== 3. is every published threshold reachable by the gates? ===\n")
# a marker is only gated on if it survives intersect(markers, select_proteins)
# and exists as a column of the ADT matrix
unreachable_sel <- setdiff(published$Protein, select_proteins)
unreachable_mat <- setdiff(published$Protein, colnames(protein_mat_normalized_lognorm))
cat("published but dropped by select_proteins:", paste(unreachable_sel, collapse = ", "), "\n")
cat("published but absent from the ADT matrix:", paste(unreachable_mat, collapse = ", "), "\n")
check(length(unreachable_sel) == 0, "published a threshold for a protein select_proteins drops")
check(length(unreachable_mat) == 0, "published a threshold for a protein absent from the ADT matrix")

cat("\n=== 4. reverse gap: markers gated on with no published threshold ===\n")
markers <- unique(trimws(unlist(strsplit(c(df_markers2$Positive, df_markers2$Negative), ","))))
markers <- intersect(markers[nzchar(markers)], select_proteins)
skipped <- sort(setdiff(markers, published$Protein))
cat(sprintf("markers in play: %d | silently skipped: %s\n", length(markers), paste(skipped, collapse = ", ")))
# Not a failure -- these 5 were deliberately dropped from the curated CSV. Kept
# visible because the GP marker signatures still name them (Extended Data
# Table 7's caption lists them). None of them appears in the 10 GP signatures
# that Figure 6e-6j / S6c-S6f actually gate on -- script/verify_gating_gps.R
# checks that, so the panel subtitles cannot be hiding a dropped marker.

cat("\n=== 5. Figure 6b/6c KLRG1 cutoff (Figure6.R, the 6b/6c block) ===\n")
klrg1 <- c(runtime = runtime$Threshold[runtime$Protein == "KLRG1"],
           published = published$Threshold[published$Protein == "KLRG1"])
print(klrg1)
check(identical(as.numeric(klrg1[["runtime"]]), as.numeric(klrg1[["published"]])), "KLRG1 cutoff mismatch")

cat("\n=== 6. curated CSV -> published table round trip ===\n")
curated <- read.csv(paste0(data_path, "Thresholds_Selected_Proteins.csv"), stringsAsFactors = FALSE)
n_curated <- sum(!is.na(curated$Threshold_manual))
joined <- merge(curated[, c("Protein", "Threshold_manual")], published, by = "Protein")
cat(sprintf("matched %d of %d curated rows; max abs diff %g\n",
            nrow(joined), n_curated, max(abs(joined$Threshold_manual - joined$Threshold))))
check(nrow(joined) == n_curated, "a curated row was lost between the CSV and the published table")
check(max(abs(joined$Threshold_manual - joined$Threshold)) == 0, "a curated value was altered")
# We follow Threshold_manual, not the second reviewer's Threshold_david column.
disagree <- curated$Protein[curated$Threshold_manual != curated$Threshold_david]
cat("columns disagree for:", paste(disagree, collapse = ", "), "-- published as",
    paste(published$Threshold[published$Protein %in% disagree], collapse = ", "),
    "(Threshold_manual, deliberate)\n")

cat("\n==================== RESULT ====================\n")
if (length(failures)) {
  cat("MISMATCH:\n", paste0("  - ", failures, collapse = "\n"), "\n", sep = "")
  quit(status = 1)
}
cat("PASS: the published thresholds are exactly the ones the figures gate on.\n")
