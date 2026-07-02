# Fit an NMF to the shifted log counts using flashier.
#
# sinteractive -p mstephens --time=24:00:00 --mem=32G -c 8
# module load R/4.2.0
# export OMP_NUM_THREADS=8
# > .libPaths()[1]
# [1] "/home/pcarbo/R_libs_4_20"
#
# Notes:
# - For k = 20, I requested 24 h and 32 GB memory.
# - For k = 25 and 30, I requested 36 h and 48 GB memory.
# - For k = 40, I requested 48 h and 64 GB memory.
#
library(tools)
library(qs)
library(Matrix)
library(ebnm)
library(flashier)
k <- 20
outfile <- sprintf("fl_nmf_k=%d.RData",k)
cat("k =",k,"\n")
cat("outfile =",outfile,"\n")
set.seed(1)
datadir <- "/project2/mstephens/immgent"
load(file.path(datadir,"flashier_snmf_ldf.RData"))
rm(fl_snmf_ldf)
shifted_log_counts <- qread(file.path(datadir,"shifted_log_counts.qs"))
n <- nrow(shifted_log_counts)
cat("n =",n,"\n")

# Set a lower bound on the variances.
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))

# Initialize an NMF using flashier's greedy initialization.
# With k = 20, this step took about 2 h.
# With k = 30, this step took about 4 h.
# With k = 40, this step took about 6 h.
t0 <- proc.time()
fl_nmf <- flash(shifted_log_counts,
                ebnm_fn = ebnm_point_exponential,
                greedy_Kmax = k,var_type = 2,S = s1,
                backfit = FALSE,verbose = 2)
t1 <- proc.time()
print(t1 - t0)

# Refine the NMF with several backfitting iterations.
# With k = 20, this step took about 10 h.
# With k = 30, this step took about 20 h.
# With k = 40, this step took about 24 h.
t0 <- proc.time()
fl_nmf <- flash_backfit(fl_nmf,extrapolate = FALSE,maxiter = 100,verbose = 2)
fl_nmf <- flash_backfit(fl_nmf,extrapolate = TRUE,maxiter = 100,verbose = 2)
t1 <- proc.time()
print(t1 - t0)

# Save the model fits to an .Rdata file.
session_info <- sessionInfo()
fl_nmf_ldf <- ldf(fl_nmf,type = "i")
save(list = c("sample_info","gene_info","fl_nmf_ldf","session_info"),
     file = outfile)
resaveRdaFiles(outfile)
