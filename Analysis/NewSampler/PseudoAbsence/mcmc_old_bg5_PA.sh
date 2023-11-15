#!/bin/bash
#SBATCH --output=cv_test_old_bg5_pa.out
#SBATCH -J old_bg5_pa
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/3.6.0
R CMD BATCH 3a_cv_old_bg5_PA.R



