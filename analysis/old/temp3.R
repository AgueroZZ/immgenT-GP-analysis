colnames(fl_nmf_ldf$L) <- paste0("k",1:20)
R <- cor(L,fl_nmf_ldf$L)
set.seed(1)
cond_factors <- c("k10","k12","k25","k68","k181")
cond_broad <- sample_info$condition_broad
cells <- which(cell_type == "CD4" | cell_type == "CD8")
cells <- sample(cells,1e5)
cells <- sort(c(cells,which(cell_type != "CD4" & cell_type != "CD8")))
p <- structure_plot(L[cells,cond_factors],gap = 10,n = 2000,
                    grouping = cond_broad[cells]) +
  labs(y = "membership",color = "",fill = "") +
  theme(legend.position = "bottom",legend.direction = "vertical")

p2 <- structure_plot(L[cells,cond_factors],gap = 10,n = 2000,
                    grouping = cell_type[cells]) +
  labs(y = "membership",color = "",fill = "") +
  theme(legend.position = "bottom",legend.direction = "vertical")
plot_grid(p,p2,nrow = 2,ncol = 1)
