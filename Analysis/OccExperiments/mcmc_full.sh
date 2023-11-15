#!/bin/bash
#SBATCH --output=mcmcfull.out
#SBATCH --job-name=mcmcfull
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=50GB

module load R/3.6.0
R CMD BATCH 2b_analysis_full_trim.R
R CMD BATCH 3a_cv_full_trim.R
R CMD BATCH 4b_trait_matching_full.R