# Verify what the CD69 GP subset actually is, and that the captions say so.
#
#   Rscript script/verify_cd69_gp_ranking.R      # exits non-zero on any mismatch
#
# cd69_top_gps_subset is a hand-picked list of 10 GPs, but two captions describe
# it and the comment beside it quotes correlation ranks. Nothing in the build
# recomputes either. Four things have to agree:
#
#   the curated list  ->  cd69_top_gps_subset   (code/R/citeseq_shared_setup.R)
#                     ->  the ranks quoted in the comment beside it
#                     ->  Fig. 6d's caption     (analysis/Figure6.Rmd)
#                     ->  Fig. S6a, b's caption (analysis/FigureS6.Rmd)
#
# The list lives in the shared setup rather than in a figure script because the
# panels that use it are split across two figures: Fig. 6d draws these GPs'
# genes, Fig. S6a/S6b draw the same GPs' mean activity, on the same axis order.
# (Before 2026-07-28 all three were Fig. 6i/6j/6k in one script.)
#
# This script recomputes the Spearman correlation between CD69 protein
# expression and GP loading over ALL GPs -- on the same cells and the same
# matrices as the figures, by sourcing the real setup -- and checks that the
# curated 10 sit among the most strongly correlated GPs (which is what the
# captions claim) without being the top 10 under any single ranking (which is
# what they must NOT claim). It also greps the pages for the old
# "most associated with CD69" phrasing.

library(Matrix) # protein matrices are dgCMatrix; must be attached for `[` to dispatch
suppressPackageStartupMessages(library(dplyr)) # the setup uses %>%

data_path <- "data/"
source("code/R/citeseq_shared_setup.R")

failures <- character(0)
check <- function(ok, msg) if (!isTRUE(ok)) failures <<- c(failures, msg)

# The curated list, read out of the shared setup rather than retyped, so an edit
# there cannot drift away from this check.
setup_file <- "code/R/citeseq_shared_setup.R"
setup_lines <- readLines(setup_file, warn = FALSE)
subset_line <- grep("^cd69_top_gps_subset <-", setup_lines, value = TRUE)
check(length(subset_line) == 1,
      sprintf("cd69_top_gps_subset is not assigned exactly once in %s", setup_file))
published <- eval(parse(text = sub("^cd69_top_gps_subset <- ", "", subset_line[1])))

# Single source of truth: no figure script may keep its own copy of the list.
for (f in c("script/Figure6.R", "script/FigureS6.R")) {
  copies <- grep("cd69_top_gps_subset <-", readLines(f, warn = FALSE), value = TRUE)
  check(length(copies) == 0,
        sprintf("%s re-assigns cd69_top_gps_subset -- it must come from %s only", f, setup_file))
}

# Same cells and same CD69 vector as Figure 6d / Figure S6a-b, over all GPs.
shared_cells_cd69 <- intersect(rownames(L_pm_filtered), rownames(protein_mat_normalized_lognorm))
cd69_expr_vec <- as.numeric(protein_mat_normalized_lognorm[shared_cells_cd69, "CD69"])
L_all <- as.matrix(L_pm_filtered[shared_cells_cd69, , drop = FALSE])
corr_all <- apply(L_all, 2, function(x) cor(x, cd69_expr_vec, method = "spearman"))

signed <- sort(corr_all, decreasing = TRUE)
abs_rank <- setNames(rank(-abs(corr_all), ties.method = "first"), names(corr_all))
pos_rank <- setNames(rank(-corr_all, ties.method = "first"), names(corr_all))
n_gp <- length(corr_all)

cat(sprintf("=== 0. inputs ===\ncells: %d   GPs: %d   curated subset: %d\n",
            length(shared_cells_cd69), n_gp, length(published)))
check(length(published) == 10, "the curated subset is no longer 10 GPs")
check(all(published %in% names(corr_all)), "a curated GP is not a column of L_pm_filtered")

cat("\n=== 1. the curated 10, with their ranks ===\n")
tbl <- data.frame(GP = published, rho = corr_all[published],
                  pos_rank = pos_rank[published], abs_rank = abs_rank[published])
print(tbl[order(tbl$abs_rank), ], row.names = FALSE)

cat("\n=== 2. caption claim: drawn from among the MOST strongly correlated ===\n")
# "Among the most" is only honest if every one of the 10 is near the top by
# |rho|. 20 of 200 is the bar: the top decile.
worst <- max(abs_rank[published])
cat(sprintf("worst |rho| rank among the 10: %d of %d (bar: top 20)\n", worst, n_gp))
check(worst <= 20, sprintf("a curated GP ranks %d of %d by |rho| -- too weak for \"among the most\"", worst, n_gp))

cat("\n=== 3. the caption must NOT claim these are THE ten most ===\n")
true_top10_pos <- names(signed)[1:10]
true_top10_abs <- names(sort(abs_rank))[1:10]
cat("true top 10 by rho   :", paste(true_top10_pos, collapse = ", "), "\n")
cat("true top 10 by |rho| :", paste(true_top10_abs, collapse = ", "), "\n")
cat("curated but outside top 10 by rho  :", paste(setdiff(published, true_top10_pos), collapse = ", "), "\n")
cat("curated but outside top 10 by |rho|:", paste(setdiff(published, true_top10_abs), collapse = ", "), "\n")
# Not a failure -- this IS the finding. It fails only if the subset silently
# becomes a real top-10, in which case the hedged caption is now the wrong text.
if (setequal(published, true_top10_pos) || setequal(published, true_top10_abs)) {
  failures <- c(failures, "the curated subset now IS a true top-10 -- the hedged caption wording is stale")
}

cat("\n=== 4. both signs are represented (the caption says \"positively or negatively\") ===\n")
neg <- published[corr_all[published] < 0]
cat("negatively correlated:", paste(sprintf("%s (rho %.3f, rank %d of %d)", neg, corr_all[neg],
                                            pos_rank[neg], n_gp), collapse = ", "), "\n")
check(length(neg) > 0, "no negatively correlated GP left in the subset -- drop \"or negatively\" from the caption")

cat("\n=== 5. ranks quoted in the shared setup's comment ===\n")
# The comment beside cd69_top_gps_subset quotes both rank vectors; a drifting
# correlation would leave those numbers wrong.
quoted_pos <- sort(pos_rank[published][corr_all[published] > 0])
quoted_abs <- sort(abs_rank[published])
cat("positive-rank vector:", paste(quoted_pos, collapse = ", "), "\n")
cat("|rho|-rank vector   :", paste(quoted_abs, collapse = ", "), "\n")
# Unwrap the comment: strip each line's "# " so a rank vector that wraps across
# two lines still reads as one string.
comment_block <- paste(sub("^#\\s*", "", setup_lines[grep("^# Curated list -- NOT a top-10", setup_lines) + 0:8]),
                       collapse = " ")
for (v in list(list(quoted_pos, "positive"), list(quoted_abs, "|rho|"))) {
  want <- paste(v[[1]], collapse = ", ")
  check(grepl(want, comment_block, fixed = TRUE),
        sprintf("%s's comment does not quote the %s ranks \"%s\"", setup_file, v[[2]], want))
}
cat("skipped by the curated list, top 10 by |rho|:",
    paste(setdiff(true_top10_abs, published), collapse = "/"), "\n")
check(grepl(paste0("skipping ", paste(setdiff(true_top10_abs, published), collapse = "/")),
            comment_block, fixed = TRUE),
      sprintf("%s's comment names the wrong skipped GPs", setup_file))

cat("\n=== 6. the old overclaiming wording is gone ===\n")
# Both captions describe the subset now: Fig. 6d (the gene heatmap) and
# Fig. S6a, b (the same GPs per tissue and per lineage).
for (f in c("analysis/Figure6.Rmd", "analysis/FigureS6.Rmd", setup_file,
            "script/Figure6.R", "script/FigureS6.R")) {
  hits <- grep("(ten|10) GPs most associated", readLines(f, warn = FALSE), value = TRUE)
  cat(sprintf("%-34s %s\n", f, if (length(hits)) paste("STILL PRESENT:", hits[1]) else "clean"))
  check(length(hits) == 0, sprintf("%s still claims these are the ten GPs most associated with CD69", f))
}

cat("\n=== 7. both captions still hedge, and both name the subset ===\n")
# A caption that dropped the "hand-picked" hedge would be overclaiming again in a
# way the grep above does not catch.
for (f in c("analysis/Figure6.Rmd", "analysis/FigureS6.Rmd")) {
  txt <- paste(readLines(f, warn = FALSE), collapse = " ")
  has_hedge <- grepl("hand-picked", txt, fixed = TRUE)
  has_source <- grepl("from among", txt, fixed = TRUE) ||
    grepl("most strongly correlated", txt, fixed = TRUE)
  cat(sprintf("%-24s hand-picked: %-5s  \"from among/most strongly correlated\": %s\n",
              f, has_hedge, has_source))
  check(has_hedge, sprintf("%s no longer calls the CD69 GP subset hand-picked", f))
  check(has_source, sprintf("%s no longer says where the CD69 GP subset was drawn from", f))
}

cat("\n==================== RESULT ====================\n")
if (length(failures)) {
  cat("MISMATCH:\n", paste0("  - ", failures, collapse = "\n"), "\n", sep = "")
  quit(status = 1)
}
cat("PASS: the Fig. 6d / S6a-b subset is 10 curated GPs from among the most CD69-correlated,\n")
cat("      it is not a true top-10, and both captions and the comment say so.\n")
