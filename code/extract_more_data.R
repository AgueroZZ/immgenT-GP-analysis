# module load R/4.1.0-no-openblas
# > .libPaths()[1]
# "/home/pcarbo/R_libs_4_10_no_openblas"
library(flashier)
library(tools)
k <- 200
datadir <- "/project2/mstephens/immgent"
load("../data/gene_info.RData")

# Extract the sample data that is useful to us.
#
#   IGT = batch
#   IGTHT = sample key
#   spleen_standard = flags spike-in cells to study batch effect
#
sample_info <- readRDS(file.path(datadir,"seurat_meta.rds"))
sample_info <- sample_info[c("sample_code","organ_simplified",
                             "IGT","IGTHT","spleen_standard",
                             "annotation_level1","annotation_level2",
                             "condition_broad","condition_detailed_organ",
                             "condition_detailed",
                             "condition_detailed_simplified")]
sample_info <- transform(sample_info,
                         sample_code       = factor(sample_code),
                         IGT               = factor(IGT),
                         IGTHT             = factor(IGTHT),
                         annotation_level1 = factor(annotation_level1),
                         annotation_level2 = factor(annotation_level2))
sample_info <-
  transform(sample_info,
    condition_broad               = factor(condition_broad),
    condition_detailed            = factor(condition_detailed),
    condition_detailed_organ      = factor(condition_detailed_organ),
    condition_detailed_simplified = factor(condition_detailed_simplified))

# Load the flashier results, then align the L and F matrices to
# sample_info and gene_info.
fit <- readRDS(file.path(datadir,"flashier_snmf.rds"))
fl_snmf_ldf <- ldf(fit,type = "i")
colnames(fl_snmf_ldf$L) <- paste0("k",1:200)
colnames(fl_snmf_ldf$F) <- paste0("k",1:200)
ids <- rownames(sample_info)
fl_snmf_ldf$L <- fl_snmf_ldf$L[ids,]
save(list = c("sample_info","gene_info","fl_snmf_ldf"),
     file = "flashier_snmf_ldf.RData")
resaveRdaFiles("flashier_snmf_ldf.RData")
