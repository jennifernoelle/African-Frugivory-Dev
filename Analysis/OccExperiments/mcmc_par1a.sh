#!/bin/bash
#SBATCH --output=mcmc_p1a.out
#SBATCH -J mcmc_p1a
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/3.6.0
R CMD BATCH 3a_cv_full_par1a.R

