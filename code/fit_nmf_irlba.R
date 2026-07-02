# 
# TO DO:
# - Add a brief description of this script at the top.
# - Update the sinteractive call.
#
# sinteractive -p mstephens --time=48:00:00 --mem=64G -c 8
# module load R/4.2.0
# export OMP_NUM_THREADS=8
# > .libPaths()[1]
# [1] "/home/pcarbo/R_libs_4_20"
library(qs)
library(Matrix)
library(ebnm)
library(flashier)
k <- 20
outfile <- "fl_nmf_full_k=20.RData"
cat("k =",k,"\n")
cat("outfile =",outfile,"\n")
set.seed(1)
datadir <- "/project2/mstephens/immgent"
# load(file.path(datadir,"irlba-k=250.RData"))
# n <- nrow(out$u)
out <- qread(file.path(datadir,"shifted_log_counts.qs"))
n <- nrow(out)
# Testing:
# i <- sample(n,1e4)
# n <- 1e4
# out$u <- out$u[i,]
cat("n =",n,"\n")

# Set a lower bound on the variances.
x  <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))

# Initialize an NMF using flashier's greedy initialization.
# With k = 20, this step took about 1 h.
t0 <- proc.time()
fl_nmf <- flash(out,greedy_Kmax = k,var_type = 2,S = s1,backfit = FALSE,
                verbose = 2)
t1 <- proc.time()
print(t1 - t0)

# Refine the NMF with several backfitting iterations.
# With k = 20, this step took about 2 h.
t0 <- proc.time()
fl_nmf <- flash_backfit(fl_nmf,extrapolate = FALSE,maxiter = 100,verbose = 2)
fl_nmf <- flash_backfit(fl_nmf,extrapolate = TRUE,maxiter = 100,verbose = 2)
t1 <- proc.time()
print(t1 - t0)

# Save the model fits to an .Rdata file.
session_info <- sessionInfo()
fl_nmf_ldf <- ldf(fl_nmf,type = "i")
save(list = c("fl_nmf_ldf","session_info"),
     file = outfile)
resaveRdaFiles(outfile)
