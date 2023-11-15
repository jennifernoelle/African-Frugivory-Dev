#!/bin/bash
#SBATCH --output=mcmc_3_20.out
#SBATCH --job-name=mcmc_3_20
#SBATCH --partition=dunsonlab --account=dunsonlab
#SBATCH --mem=10GB

module load R/3.6.0
R CMD BATCH 2d_analysis_cutfeedback_trimresults_parallel.R
