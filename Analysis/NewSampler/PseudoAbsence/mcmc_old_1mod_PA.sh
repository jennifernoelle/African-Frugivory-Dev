#!/bin/bash
#SBATCH --output=cv_old_1mod_pa.out
#SBATCH -J old_1mod_pa
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/3.6.0
R CMD BATCH 3a_cv_old_exp1mod_PA2.R



