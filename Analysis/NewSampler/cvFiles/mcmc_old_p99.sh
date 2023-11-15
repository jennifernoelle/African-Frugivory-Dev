#!/bin/bash
#SBATCH --output=cv_old_p99.out
#SBATCH -J old_p99
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 30
#SBATCH --mem-per-cpu=5GB

module load R/3.6.0
R CMD BATCH 3a_cv_old_p99.R



