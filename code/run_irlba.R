# sinteractive -p mstephens --time=48:00:00 --mem=64G -c 8
# module load R/4.2.0
# export OMP_NUM_THREADS=8
# > .libPaths()[1]
# [1] "/home/pcarbo/R_libs_4_20"
library(qs)
library(irlba)
library(flashier)
library(Matrix)
library(tools)
set.seed(1)
datadir <- "/project2/mstephens/immgent"
Y <- qread(file.path(datadir,"shifted_log_counts.qs"))
print(dim(Y))
t0 <- proc.time()
# It looks like this step takes about an hour so.
# (I'm not sure using >1 thread is helping.)
out <- irlba(Y,nv = 250,maxit = 1000,verbose = TRUE)
t1 <- proc.time()
print(t1 - t0)
save(list = "out",file = "irlba-k=250.RData")
resaveRdaFiles("irlba-k=250.RData")
