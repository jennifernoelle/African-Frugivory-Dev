#!/bin/bash
#SBATCH --output=cv_test_old_full.out
#SBATCH -J t_old_full
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/3.6.0
R CMD BATCH 3a_cv_old_bg_full.R



