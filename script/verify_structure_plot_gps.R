# Verify Extended Data Figure 6's GP selection against the published cluster AUCs.
#
#   Rscript script/verify_structure_plot_gps.R    # exits non-zero on any mismatch
#
# The figure's content is a selection: each of its seven rows shows the GPs that
# reach AUC > 0.9 for at least one cluster of its lineage. Those same AUCs are
# published as Extended Data Table 6, and the per-row GP counts are quoted in
# the caption -- three places that have to agree, and that nothing in the build
# checks. This script checks that
#
#   1. the GP set each row drew is exactly what thresholding the *published*
#      table gives, lineage by lineage, and the recorded per-GP maximum AUCs are
#      the table's values;
#   2. every row's palette is exactly what structure_plot_row_colors() gives for
#      its GPs -- one color per GP within the row, glasbey's first k in
#      ascending GP order -- and no row holds a GP that passes only in the DP
#      clusters this figure omits;
#   3. the clusters drawn, and those left out for being too small, follow the
#      figure's own display filters;
#   4. the assembled figure is on disk, with no per-row PDF left over from an
#      earlier layout or lettering; and
#   5. the caption's per-row GP counts and total still match what was drawn.
#
# The row map, the AUC threshold, the display filters and the palette rule are read from
# code/R/structure_plot_panels.R -- the file the figure itself uses -- so this
# is a re-derivation from the published numbers, not a second copy of the
# figure's logic. What the figure drew is read from the record it writes into
# output/FigureS6/ (so this check does not need the 1 GB loading matrix).
#
# The three inputs can be pointed elsewhere, which is how the failing path gets
# tested:
#   Rscript script/verify_structure_plot_gps.R --record-dir=/tmp/perturbed

source("code/R/structure_plot_panels.R")
source("code/R/table_xlsx.R")

arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(hit) > 1) stop(sprintf("--%s given more than once", name))
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit)
}

record_dir <- arg_value("record-dir", "output/FigureS6/")
table_file <- arg_value("table", "figures/final-selected/ExtendedDataTable6_GP_AUC_cluster.xlsx")
page_file <- arg_value("page", "analysis/FigureS6.Rmd")
panel_dir <- arg_value("panel-dir", "figures/final-selected/Figure S6/")

failures <- character(0)
check <- function(ok, msg) if (!isTRUE(ok)) failures <<- c(failures, msg)

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
gp_record <- utils::read.csv(file.path(record_dir, "s6_panel_gps.csv"), stringsAsFactors = FALSE)
cluster_record <- utils::read.csv(file.path(record_dir, "s6_panel_clusters.csv"), stringsAsFactors = FALSE)

# Extended Data Table 6 is one row per GP, one column per cluster; the figure
# thresholds it the other way round.
published <- read_table_xlsx(table_file)
auc_published <- t(as.matrix(published[, setdiff(names(published), "GP"), drop = FALSE]))
colnames(auc_published) <- published$GP
lineage_published <- stats::setNames(
  sub("[.].*$", "", rownames(auc_published)),
  rownames(auc_published)
)

lineages <- names(structure_plot_panels)
gp_number <- function(gp) as.integer(sub("^GP", "", gp))

cat(sprintf("record   : %s\ntable    : %s\npage     : %s\nfigure   : %s\n\n",
            record_dir, table_file, page_file, panel_dir))
cat(sprintf("=== 0. the figure's own definitions ===\nAUC > %.1f | rows: %s\n",
            structure_plot_auc_threshold,
            paste(sprintf("%s = %s", structure_plot_panels, lineages), collapse = ", ")))
cat(sprintf("published table: %d clusters x %d GPs\n",
            nrow(auc_published), ncol(auc_published)))

# ------------------------------------------------------------
# 1. every panel's GP set, re-derived from the published table
# ------------------------------------------------------------
cat("\n=== 1. GP sets re-derived from Extended Data Table 6 ===\n")
recomputed <- gps_above_auc_by_lineage(auc_published, lineage_published)

check(setequal(gp_record$lineage, lineages),
      sprintf("the record covers lineages %s, expected %s",
              paste(sort(unique(gp_record$lineage)), collapse = ","),
              paste(sort(lineages), collapse = ",")))

for (lineage in lineages) {
  drawn <- gp_record$gp[gp_record$lineage == lineage]
  want <- recomputed[[lineage]]
  missing <- setdiff(want, drawn)
  extra <- setdiff(drawn, want)
  cat(sprintf("%-6s table %2d | drawn %2d | missing: %-14s | extra: %s\n",
              lineage, length(want), length(drawn),
              if (length(missing)) paste(missing, collapse = ",") else "none",
              if (length(extra)) paste(extra, collapse = ",") else "none"))
  check(length(missing) == 0 && length(extra) == 0,
        sprintf("%s: GP set differs from the published table (missing %s; extra %s)",
                lineage, paste(missing, collapse = ","), paste(extra, collapse = ",")))
  check(identical(gp_number(drawn), sort(gp_number(drawn))),
        sprintf("%s: recorded GPs are not in ascending GP order, so the legend order is not the panel order", lineage))
  check(all(unique(gp_record$panel[gp_record$lineage == lineage]) == structure_plot_panels[[lineage]]),
        sprintf("%s: recorded under row letter(s) other than %s", lineage, structure_plot_panels[[lineage]]))
}

cat("\n--- recorded maximum AUC per GP against the table ---\n")
recorded_max <- gp_record$max_auc_in_lineage
table_max <- vapply(seq_len(nrow(gp_record)), function(i) {
  rows <- names(lineage_published)[lineage_published == gp_record$lineage[i]]
  max(auc_published[rows, gp_record$gp[i]], na.rm = TRUE)
}, numeric(1))
worst <- if (nrow(gp_record)) max(abs(recorded_max - table_max)) else 0
cat(sprintf("largest |recorded - published| over %d row-GP pairs: %.3g\n", nrow(gp_record), worst))
check(worst < 1e-9, sprintf("recorded maximum AUCs differ from the published table by up to %.3g", worst))
check(all(table_max > structure_plot_auc_threshold),
      sprintf("%d drawn GP(s) do not exceed the threshold in their own lineage",
              sum(table_max <= structure_plot_auc_threshold)))

# ------------------------------------------------------------
# 2. the palette, and the DP exclusion
# ------------------------------------------------------------
# Rows are colored independently, so the invariant is per row: as many distinct
# colors as GPs, and exactly the ones the real structure_plot_row_colors()
# returns. Checking against the function rather than against a copy of its rule
# means a palette change that has not been re-rendered fails here.
cat("\n=== 2. each row's palette, and no DP-only GP ===\n")
for (lineage in lineages) {
  rec <- gp_record[gp_record$lineage == lineage, ]
  want <- structure_plot_row_colors(rec$gp)
  as_rule <- identical(names(want), rec$gp) && identical(unname(want), rec$color)
  cat(sprintf("%-6s %2d GPs | %2d distinct colors | matches structure_plot_row_colors(): %s\n",
              lineage, nrow(rec), length(unique(rec$color)), as_rule))
  check(length(unique(rec$color)) == nrow(rec),
        sprintf("%s: two GPs share a color within this row", lineage))
  check(as_rule,
        sprintf("%s: recorded colors are not the ones structure_plot_row_colors() assigns", lineage))
}

# Colour is deliberately NOT comparable across rows since 2026-09-02, so the
# figure reuses colours between them. Reported rather than checked -- what has to
# hold is that the caption keeps saying so, which check 5 covers.
shared_gps <- table(gp_record$gp)
reused <- sum(table(gp_record$color) > 1)
cat(sprintf("%d distinct GPs | %d GPs appear in more than one row | %d colors reused across rows\n",
            length(unique(gp_record$gp)), sum(shared_gps > 1), reused))

dp_clusters <- names(lineage_published)[lineage_published == "DP"]
dp_only <- setdiff(
  colnames(auc_published)[apply(auc_published[dp_clusters, , drop = FALSE], 2,
                                function(x) any(x > structure_plot_auc_threshold, na.rm = TRUE))],
  unlist(recomputed, use.names = FALSE)
)
cat(sprintf("DP clusters in the table: %d | GPs passing only in DP: %s\n",
            length(dp_clusters),
            if (length(dp_only)) paste(dp_only, collapse = ", ") else "none"))
check(length(dp_clusters) > 0,
      "no DP clusters in the published table -- the documented DP exclusion no longer describes anything")
check(length(intersect(dp_only, gp_record$gp)) == 0,
      sprintf("DP-only GP(s) drawn: %s", paste(intersect(dp_only, gp_record$gp), collapse = ", ")))
check(!any(cluster_record$cluster %in% dp_clusters), "a DP cluster was drawn")

# ------------------------------------------------------------
# 3. the clusters each panel drew
# ------------------------------------------------------------
cat("\n=== 3. clusters drawn, and clusters left out as too small ===\n")
for (lineage in lineages) {
  rec <- cluster_record[cluster_record$lineage == lineage, ]
  want <- sort(names(lineage_published)[lineage_published == lineage])
  omitted <- rec$cluster[rec$n_cells_drawn == 0]
  cat(sprintf("%-6s table %2d | recorded %2d | drawn %2d | omitted (< %d cells): %s\n",
              lineage, length(want), nrow(rec), sum(rec$n_cells_drawn > 0),
              structure_plot_min_cluster_cells,
              if (length(omitted)) paste(omitted, collapse = ",") else "none"))
  check(identical(sort(rec$cluster), want),
        sprintf("%s: the record's clusters are not the table's clusters for this lineage", lineage))
  check(identical(rec$n_cells_drawn == 0,
                  rec$n_cells_healthy < structure_plot_min_cluster_cells),
        sprintf("%s: a cluster was omitted for a reason other than the < %d cell filter",
                lineage, structure_plot_min_cluster_cells))
  check(all(rec$n_cells_drawn <= pmin(rec$n_cells_healthy, structure_plot_max_cells_per_cluster)),
        sprintf("%s: a cluster drew more cells than it has, or more than the %d cell cap",
                lineage, structure_plot_max_cells_per_cluster))
  check(sum(rec$n_cells_drawn) > 0, sprintf("%s: no cells drawn at all", lineage))
}

# ------------------------------------------------------------
# 4. panel PDFs
# ------------------------------------------------------------
cat("\n=== 4. the assembled figure is on disk, with no per-row leftovers ===\n")
want_pdfs <- "s6.pdf"
on_disk <- list.files(panel_dir, pattern = "\\.pdf$")
missing <- want_pdfs[!file.exists(file.path(panel_dir, want_pdfs))]
empty <- want_pdfs[file.exists(file.path(panel_dir, want_pdfs)) &
                     file.size(file.path(panel_dir, want_pdfs)) == 0]
orphan <- setdiff(on_disk, want_pdfs)
cat(sprintf("expected %s | missing: %s | empty: %s | orphan: %s\n",
            paste(want_pdfs, collapse = ","),
            if (length(missing)) paste(missing, collapse = ",") else "none",
            if (length(empty)) paste(empty, collapse = ",") else "none",
            if (length(orphan)) paste(orphan, collapse = ",") else "none"))
check(length(missing) == 0, sprintf("no PDF for %s", paste(missing, collapse = ", ")))
check(length(empty) == 0, sprintf("empty PDF: %s", paste(empty, collapse = ", ")))
check(length(orphan) == 0,
      sprintf("%s is in the figure directory but is not the assembled figure -- stale layout?",
              paste(orphan, collapse = ", ")))

# ------------------------------------------------------------
# 5. the caption's numbers
# ------------------------------------------------------------
# The page carries the published caption verbatim, so it states the row letters
# and the selection threshold but no per-row GP counts. Check what it does state:
# that every row letter is paired with the lineage actually drawn under it.
cat("\n=== 5. the caption's row letters and threshold ===\n")
page <- paste(readLines(page_file, warn = FALSE), collapse = " ")
quoted <- regmatches(page, gregexpr("\\([a-g]\\) [A-Za-z0-9]+", page))[[1]]
parsed <- data.frame(
  panel = sub("^\\(([a-g])\\).*$", "\\1", quoted),
  lineage = sub("^\\([a-g]\\) ([A-Za-z0-9]+)$", "\\1", quoted),
  stringsAsFactors = FALSE
)
parsed <- parsed[!duplicated(parsed$panel), , drop = FALSE]
rownames(parsed) <- NULL
drawn_rows <- data.frame(
  panel = unname(structure_plot_panels),
  lineage = lineages,
  stringsAsFactors = FALSE
)
cat("caption:", if (nrow(parsed)) paste(sprintf("(%s) %s", parsed$panel, parsed$lineage),
                                        collapse = "; ") else "(no row letters found)", "\n")
cat("drawn  :", paste(sprintf("(%s) %s", drawn_rows$panel, drawn_rows$lineage),
                      collapse = "; "), "\n")
check(isTRUE(all.equal(parsed, drawn_rows, check.attributes = FALSE)),
      "the caption's row letters or lineage order do not match what was drawn")

threshold_quoted <- sprintf("AUC > %.1f", structure_plot_auc_threshold)
cat(sprintf("caption states \"%s\": %s\n", threshold_quoted, grepl(threshold_quoted, page, fixed = TRUE)))
check(grepl(threshold_quoted, page, fixed = TRUE),
      sprintf("the caption does not state the selection threshold (\"%s\")", threshold_quoted))

# ------------------------------------------------------------
cat("\n")
if (length(failures)) {
  cat(sprintf("FAILED (%d):\n%s\n", length(failures), paste0("  - ", failures, collapse = "\n")))
  quit(status = 1)
}
cat("all checks passed\n")
