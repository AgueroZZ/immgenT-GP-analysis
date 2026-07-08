# Build per-GP signature-gene blocks (Gene, Direction, Score) from a gene x GP
# factor score matrix: keep genes with |score| > cutoff, split by sign, rank
# each direction by |score| descending, and cap each direction at `cap` genes.
# Blocks are returned in GP-column order, up-regulated genes listed before
# down-regulated within each block. When `annotate_truncation` is TRUE (the
# default) and a direction has more than `cap` qualifying genes, one extra row
# is appended right after its (capped) gene rows -- e.g. row 101 after 100
# up-regulated genes -- noting how many more were left out, so the truncation
# is visible directly in the data instead of a side channel (header
# text/extra column); set it FALSE to silently cap with no such row. Used by
# Extended Data Table 8's long (annotate_truncation = FALSE) and wide
# (annotate_truncation = TRUE, the default) exports.
build_gp_gene_signature_blocks <- function(score_mat, cutoff = 0.1, cap = 100, annotate_truncation = TRUE) {
  build_direction <- function(vals, all_genes, direction_label) {
    shown <- head(all_genes, cap)
    genes <- shown
    dirs <- rep(direction_label, length(shown))
    scores <- unname(vals[shown])
    if (annotate_truncation && length(all_genes) > cap) {
      n_more <- length(all_genes) - cap
      genes <- c(genes, sprintf(
        "[%d more %s-regulated gene%s not shown]",
        n_more, direction_label, if (n_more == 1) "" else "s"
      ))
      dirs <- c(dirs, direction_label)
      scores <- c(scores, NA_real_)
    }
    list(genes = genes, dirs = dirs, scores = scores)
  }

  lapply(seq_len(ncol(score_mat)), function(i) {
    vals <- score_mat[, i]
    up_all <- names(vals)[vals > cutoff]
    up_all <- up_all[order(vals[up_all], decreasing = TRUE)]
    down_all <- names(vals)[vals < -cutoff]
    down_all <- down_all[order(abs(vals[down_all]), decreasing = TRUE)]

    up <- build_direction(vals, up_all, "up")
    down <- build_direction(vals, down_all, "down")

    data.frame(
      Gene = c(up$genes, down$genes),
      Direction = c(up$dirs, down$dirs),
      Score = c(up$scores, down$scores),
      stringsAsFactors = FALSE
    )
  })
}
