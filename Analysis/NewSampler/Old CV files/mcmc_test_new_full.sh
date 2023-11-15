#!/bin/bash
#SBATCH --output=cv_test_full.out
#SBATCH -J new_full
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/4.1.1-rhel8 
R CMD BATCH 3a_cv_test_bg_full.R

