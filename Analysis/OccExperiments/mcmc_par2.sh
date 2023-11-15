#!/bin/bash
#SBATCH --output=mcmc_par1.out
#SBATCH -J mcmc_p2
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/3.6.0
R CMD BATCH 3a_cv_full_par2.R

