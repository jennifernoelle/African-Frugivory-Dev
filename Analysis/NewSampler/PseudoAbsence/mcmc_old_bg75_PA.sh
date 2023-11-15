#!/bin/bash
#SBATCH --output=cv_test_old_bg75_pa2.out
#SBATCH -J old_bg75_pa2
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/3.6.0
R CMD BATCH 3a_cv_old_bg75_PA2.R



