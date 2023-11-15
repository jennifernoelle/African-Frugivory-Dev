# Performing cross-validation 


# --------- TO DO: set  your directories and name the current results files using the date---------#

## YOU REALLY HAVE TO SET THESE ##

# The directory where the analysis is performed:
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"

# Save results using convention: res_date_i.rda
date <- 'Dec14'

## THESE SHOULD BE THE SAME IF YOU CLONED THE REPO ##

# Where the processed data are saved:
data_path <- 'Data/'
# Where you want to save MCMC results:
save_path <- 'Results/'
# Where the functions are available:
source_path <- 'HelperScriptsPapadogeorgou/'


# ------ STEP 0: Some functions. --------- #

setwd(wd_path)
source(paste0(source_path, 'UpdExtraVar_function.R'))
source(paste0(source_path, 'UpdTraitCoef_function.R'))
source(paste0(source_path, 'UpdLatFac_function.R'))
source(paste0(source_path, 'UpdProbObs_function.R'))
source(paste0(source_path, 'UpdRho_function.R'))
source(paste0(source_path, 'OmegaFromV_function.R'))
source(paste0(source_path, 'useful_functions.R'))
source(paste0(source_path, 'CorrMat_function.R'))
source(paste0(source_path, 'MCMC_function.R'))
source(paste0(source_path, 'PredictInteractions_function.R'))
source(paste0(source_path, 'GetPredLatFac_function.R'))
source(paste0(source_path, 'GetPredWeights_function.R'))

library(doParallel)
library(foreach)


# --------------------------------------------------------------- #

# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv_phylo.dat'))
load(paste0(data_path, 'obs_A_mammals.dat'))
load(paste0(data_path, 'F_obs_m.dat'))
load(paste0(data_path, 'Obs_X.dat'))
load(paste0(data_path, 'Obs_W.dat'))
load(paste0(data_path, 'obs_OM.dat'))
load(paste0(data_path, 'obs_OP.dat'))

Cu <- Cu_phylo
Cv <- Cv_phylo

obs_A <- A.obs.m

# Sample sizes of the two sets of species:
nB <- nrow(Cu)
nP <- nrow(Cv)
nS <- dim(obs_A)[3]
print(c(nB, nP, nS))


# Getting the combined network for the interactions recorded in any study
# comb_A_ij = 1 if ij interact in any study, 0 otherwise
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1



# -------------- STEP 1: Specifications. ------------ #

bias_cor <- TRUE

Nsims <- 500 # 600, 500 reasonable
burn <- 2000 #20000, 2500 reasonable
thin <- 2 # 20, 2 reasonable 
use_H <- 10
theta_inf <- 0.01
mh_n_pis <- 100  # Parameter for proposal in Metropolis-Hastings for pi update.
mh_n_pjs <- 100
mh_n_rho <- 100

# Prior distributions:
stick_alpha <- 5
prior_theta <- c(1, 1)
prior_tau <- c(5, 5)
prior_rho <- c(5, 5)  # I do not update this for now.
prior_mu0 <- 0
prior_sigmasq0 <- 10
prior_sigmasq <- c(1, 1)

start_values <- NULL
sampling <- NULL


# ------------- STEP 2: Setting some interactions to out of sample ------------ #

# We highly recommend running the following code in parallel on 30 machines.

repetitions <- 10 #30

# Set up the parallel backend
cl <- parallel::makeCluster(repetitions)
doParallel::registerDoParallel(cl)

# Run your parallel loop
foreach(rr = 1:repetitions) %dopar% {
#for (rr in 1  : repetitions) {
  
  ## Create a filenames for each interation of the parallel loop
  each_filename_pred <- paste0('pred_', date, '_', as.character(rr), '.rda') 
  each_filepath_pred <- file.path(save_path, each_filename_pred)
  
  each_filename_cv <- paste0('cv_indices_', date, '_', as.character(rr), '.rda') 
  each_filepath_cv <- file.path(save_path, each_filename_cv)
  
  
  # If the file exists, skip to the next iteration
  if (file.exists(each_filepath_pred)) { # might need to check for cv too
    next
  }
  
  ## Otherwise, run your code
  set.seed(rr)
  
  # Matrix that chooses 100 recorded interactions to be held-out.
  set_out <- matrix(0, nrow = nB, ncol = nP)
  set_out[sample(which(comb_A == 1), 100)] <- 1
  
  # Zero-ing out corresponding entries in A.
  use_A <- obs_A
  use_A[which(set_out == 1)] <- 0  
  
  # Getting the indices that were zero-ed out.
  cv_indices <- matrix(NA, nrow = 100, ncol = 2)
  wh <- which(set_out == 1)
  for (ii in 1 : 100) {
    row_ii <- wh[ii] %% nB
    row_ii <- ifelse(row_ii == 0, nB, row_ii)
    col_ii <- ceiling(wh[ii] / nB)
    cv_indices[ii, ] <- c(row_ii, col_ii)
  }
  
  # Zero-ing out corresponding entries in A.
  use_A <- obs_A
  for (ii in 1 : 100) {
    use_A[cv_indices[ii, 1], cv_indices[ii, 2], ] <- 0
  }
  
  # Running the MCMC with the new recorded interaction matrix:
  set.seed(rr)
  mcmc <- MCMC(obs_A = use_A, focus = F_obs_m, occur_B = obs_OM, occur_P = obs_OP,
               obs_X = Obs_X, obs_W = Obs_W, Cu = Cu, Cv = Cv,
               Nsims = Nsims, burn = burn, thin = thin,
               use_H = use_H, bias_cor = bias_cor,
               theta_inf = theta_inf, mh_n_pis = mh_n_pis,
               mh_n_pjs = mh_n_pjs, mh_n_rho = mh_n_rho,
               stick_alpha = stick_alpha, prior_theta = prior_theta,
               prior_tau = prior_tau, prior_rho = prior_rho,
               prior_mu0 = prior_mu0, prior_sigmasq0 = prior_sigmasq0,
               prior_sigmasq = prior_sigmasq, start_values = start_values,
               sampling = sampling)
  
  # Saving the predictions:
  pred <- apply(mcmc$Ls, c(2, 3), mean)
  save(pred, file = each_filepath_pred)
  rm(pred)
  
  # Saving indices of the test data in each cv fold: 
  save(cv_indices, file = each_filepath_cv)
  rm(cv_indices)
}
  

