#!/bin/bash
#SBATCH --output=fu_new1mod_MO.out
#SBATCH -J fuN1modMO
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=5GB

module load R/4.1.1-rhel8 
R CMD BATCH 2b_analysis_full_new1mod_MO.R
