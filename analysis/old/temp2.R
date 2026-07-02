load("fl_nmf_k=40.RData")
L <- fl_nmf_ldf$L
colnames(L) <- paste0("k",1:40)
B <- model.matrix(~0 + annotation_level1,sample_info)
colnames(B) <- levels(sample_info$annotation_level1)
R <- cor(L,B)
heatmap(R)

set.seed(1)
cell_type_factors <- c("k181","k10","k58")
cell_type <- sample_info$annotation_level1
cells <- which(cell_type == "CD4" | cell_type == "CD8")
cells <- sample(cells,1e5)
cells <- sort(c(cells,which(cell_type != "CD4" & cell_type != "CD8")))
p <- structure_plot(L[cells,cell_type_factors],gap = 10,n = 2000,
                    grouping = cell_type[cells]) +
  labs(y = "membership",color = "",fill = "") +
  theme(legend.position = "bottom",legend.direction = "vertical")
