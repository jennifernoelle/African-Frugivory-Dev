#!/bin/bash
#SBATCH --output=cv_test_1mod_PA2b.out
#SBATCH -J ns_1mod_PA2b
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=40GB

module load R/4.1.1-rhel8 
R CMD BATCH 3a_cv_test_1mod_PA2b.R

