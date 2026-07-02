library(Matrix)
library(fastTopics)
library(flashier)
library(ggplot2)
library(cowplot)
load("../data/flashier_snmf_ldf.RData")
L <- fl_snmf_ldf$L
F <- fl_snmf_ldf$F
rm(fl_snmf_ldf)
L <- apply(L,2,function (x) x/quantile(x,0.999))
set.seed(1)
# cell_type_factors <- c("k1","k3","k62")
cell_type_factors <- c("k7","k29","k58","k68")
cell_type <- sample_info$annotation_level1
cells <- which(cell_type == "CD4" | cell_type == "CD8")
cells <- sample(cells,1e5)
cells <- sort(c(cells,which(cell_type != "CD4" & cell_type != "CD8")))
p <- structure_plot(L[cells,cell_type_factors],gap = 10,n = 2000,
                    grouping = cell_type[cells]) +
  labs(y = "membership",color = "",fill = "") +
  theme(legend.position = "bottom",legend.direction = "vertical")

cell_type <- sample_info$annotation_level1
y <- rank_effects_by_group_contrasts(L,cell_type)
k_set <- order(y,decreasing = TRUE)
k <- 200
p <- qplot(1:k,y[k_set]) +
  labs(x = "rank",y = "value") +
  theme_cowplot(font_size = 12)
