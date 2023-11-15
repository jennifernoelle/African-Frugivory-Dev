#!/bin/bash
#SBATCH --output=cv_new_1modalt.out
#SBATCH -J new_1modalt
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH --mem=20GB

module load R/4.1.1-rhel8 
R CMD BATCH 3a_cv_new_1mod.R

