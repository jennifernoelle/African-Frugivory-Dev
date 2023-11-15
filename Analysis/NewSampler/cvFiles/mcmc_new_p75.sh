#!/bin/bash
#SBATCH --output=cv_new_p75.out
#SBATCH -J newp75
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 30
#SBATCH --mem-per-cpu=5GB

module load R/4.1.1-rhel8 
R CMD BATCH 3a_cv_new_p75.R

