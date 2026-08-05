# Verify that the CURATED protein thresholds are exactly the ones the figures
# gate on.
#
#   Rscript script/verify_thresholds.R      # exits non-zero on any mismatch
#
# The manual positivity thresholds are hand-curated in
# data/Thresholds_Selected_Proteins.csv (see code/pipeline/03_protein_thresholds.R
# for why that file is an input, not an output). Two things have to agree, and
# nothing in the build enforces it:
#
#   the curated CSV  ->  threshold_results_subset_manual   (code/R/citeseq_shared_setup.R)
#                    ->  the threshold_df that Figure6.R and FigureS6.R gate on
#
# This script builds the second by sourcing the real setup -- not a copy of its
# logic -- and diffs it against the curated CSV read straight off disk. It also
# reports the reverse gap: markers that are gated on but have no threshold, and
# so are silently skipped (missing_threshold_action = "skip" in the gating
# helpers).
#
# The reference the runtime table is diffed against is the CSV's
# Threshold_manual column -- not Threshold, not the second reviewer's
# Threshold_david -- plus the CD62L row that citeseq_shared_setup.R appends in
# code rather than in the CSV.

library(Matrix) # protein matrices are dgCMatrix; must be attached for `[` to dispatch
suppressPackageStartupMessages(library(dplyr)) # the setup uses %>%

data_path <- "data/"
source("code/R/citeseq_shared_setup.R")

runtime <- threshold_results_subset_manual # what the gates actually receive
curated <- read.csv(paste0(data_path, "Thresholds_Selected_Proteins.csv"),
                    stringsAsFactors = FALSE)
reference <- rbind(
  data.frame(Protein = curated$Protein, Threshold = curated$Threshold_manual),
  data.frame(Protein = "CD62L", Threshold = 3)
)

failures <- character(0)
check <- function(ok, msg) if (!isTRUE(ok)) failures <<- c(failures, msg)

cat("=== 1. protein set ===\n")
cat(sprintf("runtime: %d   curated: %d\n", nrow(runtime), nrow(reference)))
cat("in runtime but not curated:", paste(setdiff(runtime$Protein, reference$Protein), collapse = ", "), "\n")
cat("curated but not in runtime:", paste(setdiff(reference$Protein, runtime$Protein), collapse = ", "), "\n")
check(setequal(runtime$Protein, reference$Protein), "protein sets differ")
check(!any(duplicated(runtime$Protein)), "duplicate protein in the runtime table")
check(!any(duplicated(reference$Protein)), "duplicate protein in the curated table")

cat("\n=== 2. values, protein by protein ===\n")
both <- merge(runtime, reference, by = "Protein", suffixes = c("_runtime", "_curated"), all = TRUE)
both$diff <- both$Threshold_runtime - both$Threshold_curated
bad <- both[!is.finite(both$diff) | both$diff != 0, ]
if (nrow(bad)) print(bad, row.names = FALSE) else cat(sprintf("all %d values identical (diff == 0)\n", nrow(both)))
check(nrow(bad) == 0, "threshold values differ between the runtime and curated tables")

cat("\n=== 3. is every curated threshold reachable by the gates? ===\n")
# a marker is only gated on if it survives intersect(markers, select_proteins)
# and exists as a column of the ADT matrix
unreachable_sel <- setdiff(reference$Protein, select_proteins)
unreachable_mat <- setdiff(reference$Protein, colnames(protein_mat_normalized_lognorm))
cat("curated but dropped by select_proteins:", paste(unreachable_sel, collapse = ", "), "\n")
cat("curated but absent from the ADT matrix:", paste(unreachable_mat, collapse = ", "), "\n")
check(length(unreachable_sel) == 0, "curated a threshold for a protein select_proteins drops")
check(length(unreachable_mat) == 0, "curated a threshold for a protein absent from the ADT matrix")

cat("\n=== 4. reverse gap: markers gated on with no curated threshold ===\n")
markers <- unique(trimws(unlist(strsplit(c(df_markers2$Positive, df_markers2$Negative), ","))))
markers <- intersect(markers[nzchar(markers)], select_proteins)
skipped <- sort(setdiff(markers, reference$Protein))
cat(sprintf("markers in play: %d | silently skipped: %s\n", length(markers), paste(skipped, collapse = ", ")))
# Not a failure -- these 5 were deliberately dropped from the curated CSV. Kept
# visible because the GP marker signatures still name them. None of them appears
# in the 10 GP signatures that Figure 6e-6j / S6c-S6f actually gate on --
# script/verify_gating_gps.R checks that, so the panel subtitles cannot be
# hiding a dropped marker.

cat("\n=== 5. Figure 6b/6c KLRG1 cutoff (Figure6.R, the 6b/6c block) ===\n")
klrg1 <- c(runtime = runtime$Threshold[runtime$Protein == "KLRG1"],
           curated = reference$Threshold[reference$Protein == "KLRG1"])
print(klrg1)
check(identical(as.numeric(klrg1[["runtime"]]), as.numeric(klrg1[["curated"]])), "KLRG1 cutoff mismatch")

cat("\n=== 6. which curated column is in force ===\n")
# We follow Threshold_manual, not the second reviewer's Threshold_david column.
disagree <- curated$Protein[curated$Threshold_manual != curated$Threshold_david]
cat("columns disagree for:", paste(disagree, collapse = ", "), "-- gated on as",
    paste(runtime$Threshold[runtime$Protein %in% disagree], collapse = ", "),
    "(Threshold_manual, deliberate)\n")
check(all(runtime$Threshold[match(disagree, runtime$Protein)] ==
            curated$Threshold_manual[match(disagree, curated$Protein)]),
      "a disputed threshold is not gated on at its Threshold_manual value")

cat("\n==================== RESULT ====================\n")
if (length(failures)) {
  cat("MISMATCH:\n", paste0("  - ", failures, collapse = "\n"), "\n", sep = "")
  quit(status = 1)
}
cat("PASS: the curated thresholds are exactly the ones the figures gate on.\n")
