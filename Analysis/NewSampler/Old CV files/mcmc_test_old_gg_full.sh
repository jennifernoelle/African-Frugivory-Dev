#!/bin/bash
#SBATCH --output=cv_test_old_full_gg.out
#SBATCH -J gg_old
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/3.6.0
R CMD BATCH 3a_cv_old_gg_full.R



