#!/bin/bash
#SBATCH --output=cv_new_1mod.out
#SBATCH -J new1mod
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/4.1.1-rhel8 
R CMD BATCH 3a_cv_new_1mod.R

