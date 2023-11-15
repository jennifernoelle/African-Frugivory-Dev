#!/bin/bash
#SBATCH --output=cv_test_pseudo_ab.out
#SBATCH -J pseudo_ab
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/4.1.1-rhel8 
R CMD BATCH 3a_cv_test_pseudoabsence.R

