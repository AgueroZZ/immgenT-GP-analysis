#!/bin/bash

#SBATCH --account=pi-mstephens
#SBATCH --time=30:00:00
#SBATCH --ntasks=1                         # Run on a single task
#SBATCH --cpus-per-task=4                  # Number of CPU cores (threads) per task
#SBATCH --mem=180G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=ziangzhang@uchicago.edu   # Email to which notifications will be sent

# Get the input data file as specified by the command-line argument.
# INPUT=${1}

# For reproducibility, add version numbers to the "module load" calls.
module load R/4.4.1

# Run the R script.
# echo "Running flashier_longer.R"
# Rscript flashier_longer.R > flashier_longer.Rout
echo "Running flashier_longer.R"
Rscript flashier_longer.R > flashier_longer.Rout
