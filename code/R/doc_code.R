# Helpers for embedding script code into the workflowr pages.
#
# The Figure*.Rmd pages show, for each panel, the exact block of script/*.R
# that produced it. Those blocks used to be selected by hard-coded line
# numbers (readLines(f)[52:100]), which silently went stale every time a
# script was edited -- a page could end up showing panel 4c's code under the
# heading for 4d. code_for() selects blocks by *name* instead, so the docs
# follow the code.
#
# A block is delimited by an anchor. Two anchor styles are recognised:
#
#   # ============================================================
#   # 4a: Standardized mean difference, activated vs resting
#   # ============================================================
#
#     the banner comments the scripts already use. Its name is the banner
#     text, so `id = "4a"` matches this one (ids are matched as a prefix of
#     the banner text).
#
#   # --- doc:3f ---
#
#     an explicit sub-anchor, for the cases where one banner section builds
#     several panels that the docs present separately. Sub-anchor lines are
#     stripped from the returned code.
#
# A block runs from its anchor to the line before the next anchor of either
# kind (or end of file). Pass several ids to concatenate their blocks:
# code_for(f, c("^", "Load data")). The id "^" is the file header -- everything
# above the first anchor.
#
# An id that matches no anchor, or more than one, is an error -- better a
# failed build than a page quietly showing the wrong code. Disambiguate by
# making the prefix longer: "4d prep" and "4d:" pick out the two blocks that
# a bare "4d" would match.
#
# Anything between
#
#   # --- internal ---
#   ...
#   # --- end internal ---
#
# is stripped from the returned code. The published pages are read by people
# outside the project, so provenance notes that only make sense internally --
# how a panel was re-lettered, which pre-refactor script it came from, what an
# earlier version of the figure showed -- live between these markers: kept in
# the source for us, kept off the page. The markers are deliberately not
# `doc:`-prefixed, so they do not register as anchors and can sit inside a
# block (including the file header) without splitting it.

# Lines like "# ==========" or "# ----------" that open/close a banner.
.doc_is_rule <- function(x) grepl("^#\\s*[=-]{5,}\\s*$", x)

# Lines like "# --- doc:3f ---".
.doc_is_subanchor <- function(x) grepl("^#\\s*-{2,}\\s*doc:\\S+\\s*-{2,}\\s*$", x)

# Lines like "# --- internal ---" / "# --- end internal ---".
.doc_is_internal_open <- function(x) grepl("^#\\s*-{2,}\\s*internal\\s*-{2,}\\s*$", x)
.doc_is_internal_close <- function(x) grepl("^#\\s*-{2,}\\s*end internal\\s*-{2,}\\s*$", x)

# Drop every "# --- internal --- ... # --- end internal ---" run, markers
# included. An unclosed marker drops the rest of the block, which is the safe
# direction to fail: too little on the page, never too much.
.doc_drop_internal <- function(x) {
  open <- .doc_is_internal_open(x)
  close <- .doc_is_internal_close(x)
  if (!any(open)) {
    if (any(close)) {
      stop("code_for(): '# --- end internal ---' with no opening marker.", call. = FALSE)
    }
    return(x)
  }
  inside <- cumsum(open) > cumsum(close)
  x[!(inside | close)]
}

.doc_subanchor_id <- function(x) sub("^#\\s*-+\\s*(doc:\\S+)\\s*-+\\s*$", "\\1", x)

# All anchors in `lines`, in file order: one row per anchor with the line it
# starts on, the line its code body starts on, and its name.
.doc_anchors <- function(lines) {
  n <- length(lines)
  rule <- .doc_is_rule(lines)
  sub <- .doc_is_subanchor(lines)

  starts <- integer(0)
  bodies <- integer(0)
  ids <- character(0)

  i <- 1L
  while (i <= n) {
    if (sub[i]) {
      starts <- c(starts, i)
      bodies <- c(bodies, i + 1L)
      ids <- c(ids, .doc_subanchor_id(lines[i]))
      i <- i + 1L
    } else if (rule[i]) {
      # Banner: opening rule, one or more comment lines, closing rule.
      j <- i + 1L
      while (j <= n && !rule[j] && grepl("^#", lines[j])) j <- j + 1L
      if (j <= n && rule[j] && j > i + 1L) {
        starts <- c(starts, i)
        bodies <- c(bodies, j + 1L)
        ids <- c(ids, sub("^#\\s*", "", lines[i + 1L]))
        i <- j + 1L
      } else {
        i <- i + 1L
      }
    } else {
      i <- i + 1L
    }
  }
  data.frame(start = starts, body = bodies, id = ids, stringsAsFactors = FALSE)
}

.doc_trim <- function(x) {
  x <- x[!.doc_is_subanchor(x)]
  x <- .doc_drop_internal(x)
  keep <- which(nzchar(trimws(x)))
  if (!length(keep)) character(0) else x[min(keep):max(keep)]
}

.doc_block <- function(lines, anchors, id) {
  if (identical(id, "^")) {
    if (!nrow(anchors)) return(.doc_trim(lines))
    return(.doc_trim(lines[seq_len(anchors$start[1L] - 1L)]))
  }
  # Sub-anchor names match exactly; banner ids match as a prefix of the
  # banner text, so "4a" finds "4a: Standardized mean difference, ...".
  hit <- if (startsWith(id, "doc:")) {
    which(anchors$id == id)
  } else {
    which(substr(anchors$id, 1L, nchar(id)) == id)
  }
  if (length(hit) == 0L) {
    stop("code_for(): no anchor matching '", id, "'. Anchors found: ",
         paste(sQuote(anchors$id), collapse = ", "), call. = FALSE)
  }
  if (length(hit) > 1L) {
    stop("code_for(): '", id, "' matches ", length(hit), " anchors (",
         paste(sQuote(anchors$id[hit]), collapse = ", "),
         "); use a longer prefix.", call. = FALSE)
  }
  from <- anchors$start[hit]
  # The block ends just before the next anchor, whichever kind.
  nxt <- anchors$start[anchors$start > anchors$body[hit] - 1L &
                         seq_len(nrow(anchors)) > hit]
  to <- if (length(nxt)) min(nxt) - 1L else length(lines)
  # Trim blank lines at either end so consecutive blocks don't gain padding.
  .doc_trim(lines[from:to])
}

#' Code block(s) for a documented panel
#'
#' @param file Path to the script, relative to the Rmd (e.g.
#'   "../script/Figure4.R").
#' @param ids One or more anchor names; blocks are returned in the order given.
#' @return A character vector of source lines, for a chunk's `code=` option.
code_for <- function(file, ids) {
  lines <- readLines(file, warn = FALSE)
  anchors <- .doc_anchors(lines)
  blocks <- lapply(ids, function(id) .doc_block(lines, anchors, id))
  unlist(Reduce(function(a, b) c(a, "", b), blocks), use.names = FALSE)
}
