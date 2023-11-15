#!/bin/bash
#SBATCH --output=fu_new75_MO.out
#SBATCH -J fuN75MO
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=5GB

module load R/4.1.1-rhel8 
R CMD BATCH 2b_analysis_full_newp75_MO.R
