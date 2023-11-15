#!/bin/bash
#SBATCH --output=cv_test_new.out
#SBATCH -J test_new
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/4.1.1-rhel8 
R CMD BATCH 3a_cv_test_badguess.R

