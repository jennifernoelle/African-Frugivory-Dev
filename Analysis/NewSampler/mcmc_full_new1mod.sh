#!/bin/bash
#SBATCH --output=full_new_1modalt.out
#SBATCH -J full_new1mod
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=5GB

module load R/4.1.1-rhel8 
R CMD BATCH 2b_analysis_full_new1mod.R
