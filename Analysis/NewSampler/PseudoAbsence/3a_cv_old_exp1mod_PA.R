# Performing cross-validation using the proposed bird-plant interaction model.

# The directory where the analysis is performed:
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"
# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
result_path <- 'Results/NewSampler/'
# Where the functions are available:
source_path <- 'HelperScriptsNew/'

setwd(wd_path)

#Save results using convention: res_date_i.rda
date <- 'July18_old_1mod_PA'

# ------ STEP 0: Some functions. --------- #

source(paste0(source_path, 'UpdExtraVar_function.R'))
source(paste0(source_path, 'UpdTraitCoef_function.R'))
source(paste0(source_path, 'UpdLatFac_function.R'))
source(paste0(source_path, 'UpdProbObs_function.R'))
source(paste0(source_path, 'UpdOccur_function.R'))
#source(paste0(source_path, 'UpdOccurP_function.R'))
source(paste0(source_path, 'UpdRho_function.R'))
source(paste0(source_path, 'OmegaFromV_function.R'))
source(paste0(source_path, 'useful_functions.R'))
source(paste0(source_path, 'CorrMat_function.R'))
source(paste0(source_path, 'MCMC_function.R'))
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

## Improved guess: Expert 1 modified
# Same study: 1, same site: 0.75 -> 0.85
# Same country and habitat: 0.5 -> 0.65, same region and habitat: 0.45 -> 0.35, same habitat only: 0.25,
# Same country not habitat: 0.1, same region not habitat: 0.05 

obs_OP <- ifelse(obs_OP == 0.75, 0.85, obs_OP)
obs_OP <- ifelse(obs_OP == 0.5, 0.65, obs_OP)
obs_OP <- ifelse(obs_OP == 0.45, 0.35, obs_OP)
#obs_OP <- ifelse(obs_OP == 0.25, 0.2, obs_OP)
#obs_OP <- ifelse(obs_OP == 0.05, 0.01, obs_OP)

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

# Useful values
nP <- ncol(obs_A)
nM <- nrow(obs_A)
nS <- dim(obs_A)[3]

## Create psuedoabsence data
O_ijs <- array(data = NA, dim = c(nM, nP, nS))

# Use binary version of occurrence probs such that they are present if at site, absent otherwise
obs_OPb <- ifelse(obs_OP ==1, 1, 0)
obs_OMb <- ifelse(obs_OM ==1, 1, 0)

dimnames(O_ijs) <- list(rownames(obs_A), colnames(obs_A), unlist(dimnames(obs_A)[3]))
for(s in 1:nS){
  for(m in 1:nM){
    for(p in 1:nP){
      O_ijs[m,p,s] <- obs_OMb[m,s] * obs_OPb[p,s] 
    }
  }
}

# Create array of indicators for observability
OF <- O_ijs * obs_F
comb_OF <- apply(OF, c(1,2), sum)
comb_OF <- (comb_OF > 1)*1 # Species with at least 2 observable interaction opportunities

# Define psuedoabsence: pairs with at least one observable interaction but no interactions observed
comb_A_neg <- ifelse(comb_A == 0, -1, comb_A)
comb_A_noint <- comb_A_neg * comb_OF
sum(comb_A_noint == -1) # only 350, so maybe only holdout 25
#View(comb_A_noint) # Pan troglodytes x Caloncaba glauca: seeds are poisonous, makes sense


# -------------- STEP 1: Specifications. ------------ #
repetitions <- 5
n.cv1 <- 100
n.cv2 <- 25

bias_cor <- TRUE  # Performing bias correction.

Nsims <- 500 #1000 #5000 # original 10000, reasonable 2500
burn <-  2000 #5000 # 2000 # 22000 # original 40000, reasonable 2500
thin <-  10 #10 #5 # original 40Nsims <-   250 #5000 # original 10000, reasonable 2500
use_H <- 10 # original 10
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

mcmc.cv.parallel <- function(rr, n.cv1, n.cv2,
                             obs_A, focus, occur_B, occur_P, obs_X, obs_W, Cu, Cv,
                             Nsims, burn, thin, use_H = 10, use_shrinkage = TRUE,
                             bias_cor = TRUE, theta_inf = 0.01,
                             mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100,
                             stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
                             prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
                             prior_sigmasq = c(1, 1), start_values = NULL,
                             sampling = NULL, cut_feed = FALSE){
  
  set.seed(rr)
  
  # Matrix that chooses n.cv1 recorded interactions to be held-out.
  set_out1 <- matrix(0, nrow = nM, ncol = nP)
  set_out1[sample(which(comb_A == 1), n.cv1)] <- 1
  
  # Matrix that chooses n.cv.2 recorded non-interactions to be held-out.
  set_out2 <- matrix(0, nrow = nM, ncol = nP)
  set_out2[sample(which(comb_A_noint == -1), n.cv2)] <- 1
  
  
  # Getting the indices that were zero-ed out: heldout interactions
  cv_indices1 <- matrix(NA, nrow = n.cv1, ncol = 2)
  wh1 <- which(set_out1 == 1)
  for (ii in 1 : n.cv1) {
    row_ii <- wh1[ii] %% nM
    row_ii <- ifelse(row_ii == 0, nM, row_ii)
    col_ii <- ceiling(wh1[ii] / nM)
    cv_indices1[ii, ] <- c(row_ii, col_ii)
  }
  
  # Getting the indices that were zero-ed out: heldout non-interactions
  cv_indices2 <- matrix(NA, nrow = n.cv2, ncol = 2)
  wh2 <- which(set_out2 == 1)
  for (ii in 1 : n.cv2) {
    row_ii <- wh2[ii] %% nM
    row_ii <- ifelse(row_ii == 0, nM, row_ii)
    col_ii <- ceiling(wh2[ii] / nM)
    cv_indices2[ii, ] <- c(row_ii, col_ii)
  }
  
  # Zero-ing out corresponding entries in A and F for heldout interactions
  use_A <- obs_A
  use_F <- obs_F
  for (ii in 1 : n.cv1) {
    use_A[cv_indices1[ii, 1], cv_indices1[ii, 2], ] <- 0
    use_F[cv_indices1[ii, 1], cv_indices1[ii, 2], ] <- 0
  }
  
  # Zero-ing out corresponding entries in F for heldout non-interactions
  for (ii in 1 : n.cv2) {
    #use_A[cv_indices[ii, 1], cv_indices[ii, 2], ] <- 0 # A is already 0
    use_F[cv_indices2[ii, 1], cv_indices2[ii, 2], ] <- 0
  }
  
  # Running the MCMC with the new recorded interaction matrix:
  mcmc <- MCMC(obs_A = use_A, focus = use_F, occur_B = obs_OM, occur_P = obs_OP,
               obs_X = obs_X, obs_W = obs_W, Cu = Cu, Cv = Cv,
               Nsims = Nsims, burn = burn, thin = thin,
               use_H = use_H, bias_cor = bias_cor,use_shrinkage = TRUE,
               theta_inf = theta_inf, mh_n_pis = mh_n_pis,
               mh_n_pjs = mh_n_pjs, mh_n_rho = mh_n_rho,
               stick_alpha = stick_alpha, prior_theta = prior_theta,
               prior_tau = prior_tau, prior_rho = prior_rho,
               prior_mu0 = prior_mu0, prior_sigmasq0 = prior_sigmasq0,
               prior_sigmasq = prior_sigmasq, start_values = start_values,
               sampling = sampling, cut_feed = FALSE)
  
  # Saving the predictions:
  pred <- apply(mcmc$Ls, c(2, 3), mean) # Why not save the interaction probs instead of the predictions?
  save(pred, file = paste0(result_path, 'pred_', date, '_', rr, '.dat'))
  rm(pred)
  
  # Saving the cv indices
  save(cv_indices1, file = paste0(result_path, 'cv_indices1_', date, '_', rr, '.dat'))
  save(cv_indices2, file = paste0(result_path, 'cv_indices2_', date, '_', rr, '.dat'))
  rm(cv_indices1)
  rm(cv_indices2)

}

#---------------------- STEP 3: RUN THE SAMPLER -------------------------------------------#

# We highly recommend running the following code in parallel on 30 machines.
repetitions <- 5

t1 <- Sys.time()
mclapply(1:repetitions, function(i) mcmc.cv.parallel(rr=i, n.cv1 = n.cv1, n.cv2 = n.cv2,
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


