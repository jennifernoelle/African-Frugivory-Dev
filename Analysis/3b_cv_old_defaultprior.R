# Performing cross-validation using the proposed bird-plant interaction model.

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
result_path <- 'Results/'
# Where the functions are available:
source_path <- 'HelperScriptsPapadogeorgouWithFixes/'

#Save results using convention: res_date_i.rda
date <- 'Old_default_prior'

# ------ STEP 0: Some functions. --------- #

source(paste0(source_path, 'UpdExtraVar_function.R'))
source(paste0(source_path, 'UpdTraitCoef_function.R'))
source(paste0(source_path, 'UpdLatFac_function.R'))
source(paste0(source_path, 'UpdProbObs_function.R'))
source(paste0(source_path, 'UpdOccur_function.R'))
source(paste0(source_path, 'UpdRho_function.R'))
source(paste0(source_path, 'OmegaFromV_function.R'))
source(paste0(source_path, 'useful_functions.R'))
source(paste0(source_path, 'CorrMat_function.R'))
source(paste0(source_path, 'MCMC_function_trim.R'))
source(paste0(source_path, 'PredictInteractions_function.R'))
source(paste0(source_path, 'GetPredLatFac_function.R'))
source(paste0(source_path, 'GetPredWeights_function.R'))

library(parallel)
library(abind)
library(tidyverse)
#library(msm) # this one has log option for truncnorm

# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv_phylo.dat'))
load(paste0(data_path, 'obs_A_mammals.dat'))
load(paste0(data_path, 'F_obs_default.dat'))
load(paste0(data_path, 'Obs_X.dat')) # mammal traits
load(paste0(data_path, 'Obs_W.dat')) # plant traits
load(paste0(data_path, 'obs_OM.dat')) # mammal occurrences: many species have only one study
load(paste0(data_path, 'obs_OP.dat')) # plant occurrences: most species have many studies, but some have 0?
load(paste0(data_path, 'traits_p_709_clean.dat')) # plant occurrences: most species have many studies, but some have 0?
load(paste0(data_path, 'Subset_no_baboons.dat')) # list with appropriately subset mammals, plants, studies

## Rename for convenience
Cu <- Cu_phylo
Cv <- Cv_phylo
obs_A <- A.obs.m
obs_X <- Obs_X
obs_W <- Obs_W

# Try simplifying traits and getting rid of NAs
obs_W <- obs_W[, 1:2]

## Subset to remove baboons
wh_keep_m <- which(rownames(obs_A) %in% nobab.list[[1]])
#common.plants <- names(which(apply(obs_A, 2, sum)>10))
wh_keep_p <- which(colnames(obs_A) %in% nobab.list[[2]])
wh_keep_s <- which(unlist(dimnames(obs_A)[3]) %in% nobab.list[[3]])

obs_A <- obs_A[wh_keep_m, wh_keep_p, wh_keep_s]
obs_F <- obs_F[wh_keep_m, wh_keep_p, wh_keep_s]
obs_OP <- obs_OP[wh_keep_p, wh_keep_s]
obs_OM <- obs_OM[wh_keep_m, wh_keep_s]
obs_X <- Obs_X[wh_keep_m, ]
obs_W <- Obs_W[wh_keep_p, ]
Cu <- Cu[wh_keep_m, wh_keep_m]
Cv <- Cv[wh_keep_p, wh_keep_p]


## The following assignments were used in creating obs_OP
# Same study: 1, same site: 0.75
# Same country and habitat: 0.5, same region and habitat: 0.45, same habitat only: 0.25,
# Same country not habitat: 0.1, same region not habitat: 0.05

## Improved default guess: 0.75/1

obs_OP <- ifelse(obs_OP == 1, 1, 0.75)

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

# Useful values
nP <- ncol(obs_A)
nM <- nrow(obs_A)
nS <- dim(obs_A)[3]


# -------------- STEP 1: Specifications. ------------ #

bias_cor <- TRUE  # Performing bias correction.

Nsims <- 5 #10000 recommended
burn <-  1 #10000 recommended
thin <-  1 #10-40 recommended
use_H <- 10 
theta_inf <- 0.01
mh_n_pis <- 70  # Parameter for proposal in Metropolis-Hastings for pi update.
mh_n_pjs <- 70
mh_n_rho <- 100

# Prior distributions:
stick_alpha <- 5
prior_theta <- c(1, 1)
prior_tau <- c(5, 5)
prior_rho <- c(5, 5)  
prior_mu0 <- 0
prior_sigmasq0 <- 10
prior_sigmasq <- c(1, 1)

sampling <- NULL
start_values <- NULL


# ------------- STEP 2: Parallel setup --------------- #


# Define the cv function to run in parallel 

mcmc.cv.parallel <- function(rr, n.cv, 
                             obs_A, focus, occur_B, occur_P, obs_X, obs_W, Cu, Cv,
                             Nsims, burn, thin, use_H = 10, use_shrinkage = TRUE,
                             bias_cor = TRUE, theta_inf = 0.01,
                             mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100,
                             stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
                             prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
                             prior_sigmasq = c(1, 1), start_values = NULL,
                             sampling = NULL, cut_feed = FALSE, trim = TRUE){
  
  set.seed(rr)
  
  # Matrix that chooses n.cv recorded interactions to be held-out.
  set_out <- matrix(0, nrow = nM, ncol = nP)
  set_out[sample(which(comb_A == 1), n.cv)] <- 1
  
  
  # Getting the indices that were zero-ed out.
  cv_indices <- matrix(NA, nrow = n.cv, ncol = 2)
  wh <- which(set_out == 1)
  for (ii in 1 : n.cv) {
    row_ii <- wh[ii] %% nM
    row_ii <- ifelse(row_ii == 0, nM, row_ii)
    col_ii <- ceiling(wh[ii] / nM)
    cv_indices[ii, ] <- c(row_ii, col_ii)
  }
  
  # Zero-ing out corresponding entries in A. 
  use_A <- obs_A
  for (ii in 1 : n.cv) {
    use_A[cv_indices[ii, 1], cv_indices[ii, 2], ] <- 0
  }
  
  # Shouldn't we also set Focus to zero for these observations?
  use_F <- obs_F
  for (ii in 1 : n.cv) {
    use_F[cv_indices[ii, 1], cv_indices[ii, 2], ] <- 0
  }
  
  # Running the MCMC with the new recorded interaction matrix:
  mcmc <- MCMC.trim(obs_A = use_A, focus = use_F, occur_B = obs_OM, occur_P = obs_OP,
               obs_X = obs_X, obs_W = obs_W, Cu = Cu, Cv = Cv,
               Nsims = Nsims, burn = burn, thin = thin,
               use_H = use_H, bias_cor = bias_cor,use_shrinkage = use_shrinkage,
               theta_inf = theta_inf, mh_n_pis = mh_n_pis,
               mh_n_pjs = mh_n_pjs, mh_n_rho = mh_n_rho,
               stick_alpha = stick_alpha, prior_theta = prior_theta,
               prior_tau = prior_tau, prior_rho = prior_rho,
               prior_mu0 = prior_mu0, prior_sigmasq0 = prior_sigmasq0,
               prior_sigmasq = prior_sigmasq, start_values = start_values,
               sampling = sampling, cut_feed = cut_feed, trim = trim)
  
  # Saving the predictions:
  pred <- apply(mcmc$Ls, c(2, 3), mean) # Why not save the interaction probs instead of the predictions?
  save(pred, file = paste0(result_path, 'pred_', date, '_', rr, '.dat'))
  rm(pred)
  
  save(cv_indices, file = paste0(result_path, 'cv_indices_', date, '_', rr, '.dat'))
  rm(cv_indices)

}

#---------------------- STEP 3: RUN THE SAMPLER -------------------------------------------#

# We highly recommend running the following code in parallel on 30 machines.
repetitions <- 30
n.cv <- 100

t1 <- Sys.time()
mclapply(1:repetitions, function(i) mcmc.cv.parallel(rr=i, n.cv = n.cv, 
                                               obs_A = obs_A, focus = obs_F, occur_B = obs_OM, occur_P = obs_OP,
                                               obs_X = obs_X, obs_W = obs_W, Cu = Cu, Cv = Cv,
                                               Nsims = Nsims, burn = burn, thin = thin,
                                               use_H = use_H, bias_cor = bias_cor,use_shrinkage = TRUE,
                                               theta_inf = theta_inf, mh_n_pis = mh_n_pis,
                                               mh_n_pjs = mh_n_pjs, mh_n_rho = mh_n_rho,
                                               stick_alpha = stick_alpha, prior_theta = prior_theta,
                                               prior_tau = prior_tau, prior_rho = prior_rho,
                                               prior_mu0 = prior_mu0, prior_sigmasq0 = prior_sigmasq0,
                                               prior_sigmasq = prior_sigmasq, start_values = start_values,
                                               sampling = sampling, cut_feed = FALSE),
         mc.cores = repetitions)



Sys.time() - t1


