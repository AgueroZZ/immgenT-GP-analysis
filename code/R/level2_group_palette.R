# The annotation_level2_group order and palette.
#
# ZemmourLib carries canonical palettes for annotation_level1 and
# annotation_level2 but none for annotation_level2_group, so the order and the
# colours are defined here. They are deliberately a different colour family
# from the level1 primaries, so a level1 bar and a level2_group bar drawn one
# above the other cannot be confused, and none is near-white (against a white
# heatmap body a pale category would read as missing data rather than as a
# level).
#
# Shared by script/Figure1.R (panel 1D, cells on columns) and
# script/FigureS2d.R (Extended Data Figure 2d, level2 clusters on columns), so
# the two heatmaps cannot drift apart on what "activated" is coloured.
#
# EXCLUDE_LEVEL2_GROUPS is the group both panels drop: "miniverse" is exactly
# the seven ".wM" level2 clusters (CD8.wM, CD4.wM, Treg.wM, gdT.wM, CD8aa.wM,
# Tz.wM, DN.wM).

LEVEL2_GROUP_ORDER <- c("resting", "activated", "proliferating", "miniverse", "other")
LEVEL2_GROUP_COLORS <- c(
  resting       = "#80cdc1",
  activated     = "#b2182b",
  proliferating = "#542788",
  miniverse     = "#8c510a",
  other         = "#bdbdbd"
)
EXCLUDE_LEVEL2_GROUPS <- c("miniverse")

# One level2_group per level2 cluster, in the order `groups` is given. Errors
# if a cluster carries more than one group label, or a label outside
# LEVEL2_GROUP_ORDER -- both panels rely on the map being a function of the
# cluster alone.
level2_to_group_map <- function(meta, groups) {
  mapping <- unique(data.frame(
    group = as.character(meta$annotation_level2),
    level2_group = as.character(meta$annotation_level2_group),
    stringsAsFactors = FALSE
  ))
  if (anyDuplicated(mapping$group)) {
    stop("Each annotation_level2 label must map to exactly one annotation_level2_group label.")
  }

  cluster_group <- mapping$level2_group[match(groups, mapping$group)]
  names(cluster_group) <- groups
  if (anyNA(cluster_group) || any(!cluster_group %in% LEVEL2_GROUP_ORDER)) {
    stop("Every displayed level2 cluster must carry a known annotation_level2_group.")
  }
  cluster_group
}

# Figure 1d's column order, applied to level2 clusters: lineage first, then
# level2_group as a contiguous block ranked by the alphabetically-lowest
# cluster it contains, then the clusters alphabetically inside each block.
#
# The point of the middle key is the ".P" (proliferating) cluster: it is the
# only one whose letter would otherwise drop it mid-alphabet and split its
# lineage's activated/other run in two. Ranking the block by its first cluster
# puts "proliferating" (whose first and only cluster is ".P") at the end of the
# lineage instead. LEVEL2_GROUP_ORDER fixes the *legend* order, not this one.
#
# Returns a permutation of seq_along(groups).
#
# --- internal ---
# script/Figure1.R applies the same rule inline, to cells rather than clusters,
# and ranks the blocks over the cell set it actually sampled rather than over a
# cluster list -- which is why it does not call this. Keep the two in step.
# --- end internal ---
level2_group_block_order <- function(groups, group_level1, group_l2group, level1_order) {
  lineage <- unname(group_level1[groups])
  l2group <- unname(group_l2group[groups])
  if (anyNA(lineage) || anyNA(l2group)) {
    stop("Every group needs a level1 and an annotation_level2_group label.")
  }
  if (!all(lineage %in% level1_order) || !all(l2group %in% LEVEL2_GROUP_ORDER)) {
    stop("Block ordering saw a level1 or annotation_level2_group label it does not know.")
  }

  lineage_rank <- match(lineage, level1_order)
  # Rank of each cluster under (lineage, alphabetical) -- the order the blocks
  # are then built from.
  cluster_rank <- order(order(lineage_rank, groups))
  block <- paste(lineage_rank, l2group, sep = "|")
  block_rank <- as.integer(tapply(cluster_rank, block, min)[block])

  ordering <- order(lineage_rank, block_rank, cluster_rank)
  if (anyDuplicated(rle(block[ordering])$values)) {
    stop("A lineage x level2_group block is not contiguous after ordering.")
  }
  ordering
}
