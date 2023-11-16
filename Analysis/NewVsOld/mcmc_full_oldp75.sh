#!/bin/bash
#SBATCH --output=full_old_p75.out
#SBATCH -J full_oldp75
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=5GB

module load  R/3.6.0
R CMD BATCH 2b_analysis_full_oldp75.R
