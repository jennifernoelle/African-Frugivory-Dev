#!/bin/bash
#SBATCH --output=full_new_p75.out
#SBATCH -J full_newp75
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=5GB

module load R/4.1.1-rhel8 
R CMD BATCH 2b_analysis_full_newp75.R
