#!/bin/bash
#SBATCH --output=cv_test_bg5_PA.out
#SBATCH -J ns_bg5_PA
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/4.1.1-rhel8 
R CMD BATCH 3a_cv_test_bg5_PA.R

