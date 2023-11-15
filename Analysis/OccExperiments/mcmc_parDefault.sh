#!/bin/bash
#SBATCH --output=mcmc_Default.out
#SBATCH -J mcmc_default
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/3.6.0
R CMD BATCH 3a_cv_full_parDefault.R

