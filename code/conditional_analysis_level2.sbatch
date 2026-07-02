#!/bin/bash

#SBATCH --account=pi-mstephens
#SBATCH --time=16:00:00
#SBATCH --ntasks=1                         # Run on a single task
#SBATCH --cpus-per-task=1                  # Number of CPU cores (threads) per task
#SBATCH --mem=50G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ;
#SBATCH --mail-user=ziangzhang@uchicago.edu   # Email to which notifications will be sent

# Get the input data file as specified by the command-line argument.
# INPUT=${1}

# For reproducibility, add version numbers to the "module load" calls.
module load R/4.4.1

# Run the R script.
# echo "Running conditional_analysis_level2.R"
# Rscript conditional_analysis_level2.R > conditional_analysis_level2.Rout
echo "Running conditional_analysis_level2.R"
Rscript conditional_analysis_level2.R > conditional_analysis_level2.Rout
