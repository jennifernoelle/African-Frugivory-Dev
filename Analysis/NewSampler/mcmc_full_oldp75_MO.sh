#!/bin/bash
#SBATCH --output=fu_old_p75_MO.out
#SBATCH -J fuO75MO
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=5GB

module load  R/3.6.0
R CMD BATCH 2b_analysis_full_oldp75_MO.R
