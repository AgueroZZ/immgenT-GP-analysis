### What cells are removed from the iterative filtering process?

library(ggplot2)
library(dplyr)
library(fastTopics)
library(qs)
library(cowplot)
library(ggrepel)
data_path <- "./data/"
code_path <- "./code/"
source(paste0(code_path, "compute_cosine_sim.R"))
source(paste0(code_path, "filtering_membership.R"))
seurat_meta <- readRDS(paste0(data_path, "seurat_meta.rds"))
flashier_snmf_summary <- readRDS(paste0(data_path, "flashier_snmf_summary.rds"))

L_pm <- flashier_snmf_summary$L_pm
d <- apply(L_pm,2,max)
L_pm <- scale_cols(L_pm,1/d)

F_pm <- flashier_snmf_summary$F_pm
F_pm <- scale_cols(F_pm,d)

# filtering membership
filtering_membership <- filter_cells_by_total_membership(L_pm, numiter = 12)
L_pm_filtered <- L_pm[filtering_membership,]
d_filtered <- apply(L_pm_filtered,2,max)
L_pm_filtered <- scale_cols(L_pm_filtered,1/d_filtered)

F_pm_filtered <- F_pm
F_pm_filtered <- scale_cols(F_pm_filtered,d_filtered)

# what cells are removed?
removed_cells <- rownames(L_pm)[-filtering_membership]
removed_cells_meta <- seurat_meta[removed_cells,]


# produce a table of annotation_level1 proportions and counts for removed cells and all cells
removed_cells_meta %>%
  group_by(annotation_level1) %>%
  summarise(count = n()) %>%
  mutate(proportion_removed = count / sum(count), count_removed = count) %>%
  arrange(desc(proportion_removed)) %>%
  left_join(seurat_meta %>%
              group_by(annotation_level1) %>%
              summarise(count = n()) %>%
              mutate(proportion_all = count / sum(count), count_all = count) %>%
              arrange(desc(proportion_all)),
            by = "annotation_level1") %>%
  select(annotation_level1, proportion_removed, proportion_all, count_removed, count_all) %>%
  arrange(desc(proportion_removed)) %>%
  print(n = Inf)


# produce a table of annotation_level2 proportions and counts for removed cells and all cells
removed_cells_meta %>%
  group_by(annotation_level2_group) %>%
  summarise(count = n()) %>%
  mutate(proportion_removed = count / sum(count), count_removed = count) %>%
  arrange(desc(proportion_removed)) %>%
  left_join(seurat_meta %>%
              group_by(annotation_level2_group) %>%
              summarise(count = n()) %>%
              mutate(proportion_all = count / sum(count), count_all = count) %>%
              arrange(desc(proportion_all)),
            by = "annotation_level2_group") %>%
  select(annotation_level2_group, proportion_removed, proportion_all, count_removed, count_all) %>%
  arrange(desc(proportion_removed)) %>%
  print(n = Inf)

removed_cells_meta %>%
  group_by(annotation_level2) %>%
  summarise(count = n()) %>%
  mutate(proportion_removed = count / sum(count), count_removed = count) %>%
  arrange(desc(proportion_removed)) %>%
  left_join(seurat_meta %>%
              group_by(annotation_level2) %>%
              summarise(count = n()) %>%
              mutate(proportion_all = count / sum(count), count_all = count) %>%
              arrange(desc(proportion_all)),
            by = "annotation_level2") %>%
  select(annotation_level2, proportion_removed, proportion_all, count_removed, count_all) %>%
  arrange(desc(proportion_removed)) %>%
  print(n = Inf)

# CD4_cl13 and gdT_cl10 are very enriched in the removed cells



# produce a table of organ proportions and counts for removed cells and all cells
removed_cells_meta %>%
  group_by(organ) %>%
  summarise(count = n()) %>%
  mutate(proportion_removed = count / sum(count), count_removed = count) %>%
  arrange(desc(proportion_removed)) %>%
  left_join(seurat_meta %>%
              group_by(organ) %>%
              summarise(count = n()) %>%
              mutate(proportion_all = count / sum(count), count_all = count) %>%
              arrange(desc(proportion_all)),
            by = "organ") %>%
  select(organ, proportion_removed, proportion_all, count_removed, count_all) %>%
  arrange(desc(proportion_removed)) %>%
  print(n = Inf)


# produce a table of condition_detailed proportions and counts for removed cells and all cells
removed_cells_meta %>%
  group_by(condition_detailed) %>%
  summarise(count = n()) %>%
  mutate(proportion_removed = count / sum(count), count_removed = count) %>%
  arrange(desc(proportion_removed)) %>%
  left_join(seurat_meta %>%
              group_by(condition_detailed) %>%
              summarise(count = n()) %>%
              mutate(proportion_all = count / sum(count), count_all = count) %>%
              arrange(desc(proportion_all)),
            by = "condition_detailed") %>%
  select(condition_detailed, proportion_removed, proportion_all, count_removed, count_all) %>%
  arrange(desc(proportion_removed)) %>%
  print(n = Inf)


# produce a table of IGTHT proportions and counts for removed cells and all cells
removed_cells_meta %>%
  group_by(IGTHT) %>%
  summarise(count = n()) %>%
  mutate(proportion_removed = count / sum(count), count_removed = count) %>%
  arrange(desc(proportion_removed)) %>%
  left_join(seurat_meta %>%
              group_by(IGTHT) %>%
              summarise(count = n()) %>%
              mutate(proportion_all = count / sum(count), count_all = count) %>%
              arrange(desc(proportion_all)),
            by = "IGTHT") %>%
  select(IGTHT, proportion_removed, proportion_all, count_removed, count_all) %>%
  arrange(desc(proportion_removed)) %>%
  print(n = Inf)


# produce a table of lab proportions and counts for removed cells and all cells
removed_cells_meta %>%
  group_by(lab) %>%
  summarise(count = n()) %>%
  mutate(proportion_removed = count / sum(count), count_removed = count) %>%
  arrange(desc(proportion_removed)) %>%
  left_join(seurat_meta %>%
              group_by(lab) %>%
              summarise(count = n()) %>%
              mutate(proportion_all = count / sum(count), count_all = count) %>%
              arrange(desc(proportion_all)),
            by = "lab") %>%
  select(lab, proportion_removed, proportion_all, count_removed, count_all) %>%
  arrange(desc(proportion_removed)) %>%
  print(n = Inf)





