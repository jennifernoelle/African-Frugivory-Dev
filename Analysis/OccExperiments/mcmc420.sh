#!/bin/bash
#SBATCH --output=mcmc420.out
#SBATCH --job-name=mcmc420
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=20GB

module load R/3.6.0
R CMD BATCH 2b_analysis_subset.R
R CMD BATCH 3a_cross_validation.R
R CMD BATCH 4_trait_matching.R