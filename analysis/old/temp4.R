set.seed(1)
cond_factors <- c("k10","k58")
cond_detailed <- sample_info$condition_detailed
cells <- which(sample_info$annotation_level1 == "CD8" & 
               is.element(sample_info$condition_broad,
                          c("parasite","bacteria","virus")))
cells <- which(with(sample_info,
                    annotation_level1 == "CD8" & 
                    !(grepl("LCM",condition_detailed,fixed = TRUE) |
                      grepl("MCMV",condition_detailed,fixed = TRUE) |
                      grepl("Toxo",condition_detailed,fixed = TRUE)) &
                    is.element(condition_broad,
                               c("parasite","bacteria","virus"))))
p <- structure_plot(L[cells,cond_factors],gap = 10,n = 2000,
                    grouping=factor(sample_info$annotation_level2[cells])) +
  labs(y = "membership",color = "",fill = "") +
  theme(legend.position = "bottom",legend.direction = "vertical")

