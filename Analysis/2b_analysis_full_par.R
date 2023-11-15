# This version doesn't save full probabilities but rather running mean
# Doesn't save marginal probabilities or latent factors at all

# --------- TO DO: set  your directories and name the current results files using the date---------#

## YOU REALLY HAVE TO SET THESE ##

# The directory where the analysis is performed:
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"
setwd(wd_path)

# Save results using convention: res_date_i.rda
date <- 'May10_FullOcc4_par'
date <- 'TEST_par'
batchname <- "mcmc_full"


## THESE SHOULD BE THE SAME IF YOU CLONED THE REPO ##

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
save_path <- 'Results/'
# Where the functions are available:
source_path <- 'HelperScriptsNew/'

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

#library(doParallel)
library(parallel)
library(foreach)
library(abind)
library(magrittr)


# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv_phylo.dat'))
load(paste0(data_path, 'obs_A_mammals.dat'))
load(paste0(data_path, 'F_obs_default.dat'))
load(paste0(data_path, 'Obs_X.dat')) # mammal traits
load(paste0(data_path, 'Obs_W.dat')) # plant traits
load(paste0(data_path, 'obs_OM.dat')) # site level obs_OM
load(paste0(data_path, 'obs_OP.dat')) # site level obs_OP

# Data for subsetting

## Rename for convenience
Cu <- Cu_phylo
Cv <- Cv_phylo
obs_A <- A.obs.m
obs_W <- Obs_W[, 1:2]
obs_X <- Obs_X

# Useful values
nP <- ncol(obs_A)
nM <- nrow(obs_A)
nS <- dim(obs_A)[3]

## The following assignments were used in creating obs_OP
# Same study: 1, same site: 0.75
# Same country and habitat: 0.5, same region and habitat: 0.45, same habitat only: 0.25,
# Same country not habitat: 0.1, same region not habitat: 0.05

# # Create alternate expert-defined occurrence
# obs_OP <- ifelse(obs_OP == 0.75, 0.85, obs_OP) # boost same site P(occ) = 0.85 #2
# obs_OP <- ifelse(obs_OP == 0.5, 0.65, obs_OP) # boost country-habitat P(occ) = 0.65 #2
# obs_OP <- ifelse(obs_OP == 0.25, 0.35, obs_OP) # boost same habitat P(occ) = 0.45 #3 is 0.35 for #4
# obs_OP <- ifelse(obs_OP == 0, 0.1, obs_OP) #4

# -------------- STEP 1: Specifications. ------------ #

bias_cor <- TRUE # Performing bias correction.

Nsims <- 500 #1000 #5000 # original 10000, reasonable 2500
burn <- 2000 #5000 # 2000 # 22000 # original 40000, reasonable 2500
thin <-  10 #10 #5 # original 40
use_H <- 10 # original 10
theta_inf <- 0.01
mh_n_pis <- 70  # Parameter for proposal in Metropolis-Hastings for pi update.
mh_n_pjs <- 70
mh_n_rho <- 100

# Prior distributions:
stick_alpha <- 5
prior_theta <- c(1, 1)
prior_tau <- c(5, 5)
prior_rho <- c(5, 5)  # I do not update this for now.
prior_mu0 <- 0
prior_sigmasq0 <- 10
prior_sigmasq <- c(1, 1)

sampling <- NULL
start_values <- NULL

# --------------- STEP 2: SETUP PARALLEL ----------------- #

# Set up function to execute in parallel
n.chains <- 4
mcmc.parallel <- function(cc, obs_A, focus, occur_B, occur_P, obs_X, obs_W, Cu, Cv,
                               Nsims, burn, thin, use_H = 10, use_shrinkage = TRUE,
                               bias_cor = TRUE, theta_inf = 0.01,
                               mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100,
                               stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
                               prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
                               prior_sigmasq = c(1, 1), start_values = NULL,
                               sampling = NULL, cut_feed = FALSE){
  
  set.seed(cc) 
  
  ## Create a unique filename for each interation of the parallel loop
  each_filename <- paste0('testruns/res_', date, '_', as.character(cc), '.dat') 
  each_filepath <- file.path(save_path, each_filename)
  
  mcmc <- MCMC.trim(obs_A = obs_A, focus = obs_F, occur_B = obs_OM, occur_P = obs_OP,
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
  
  
  # Binding different predictions of interest: Posterior samples of the
  # interaction indicators, the linear predictor of the interaction model,
  all_pred <- abind::abind(pred_L = mcmc$Ls, probL = mcmc$mod_pL1s, along = 4)
  
  # Phylogenetic correlation parameter for bird and plant correlation matrices.
  correlations <- cbind(U = mcmc$rU, V = mcmc$rV)
  
  # Running mean of detection probabilities
  p_detect <- list(pis = mcmc$pi_mean, pjs = mcmc$pj_mean)
  
  # Running mean of latent factors
  factors <- list(U = mcmc$U_mean, V = mcmc$V_mean)
  
  # Combining the results we are interested in to a list and saving:
  res <- list(all_pred = all_pred, correlations = correlations, p_detect = p_detect, factors = factors)
  save(res, file = each_filepath)
  
  rm(res)
}
# 
# s1 <- system.time(mapply(function(i) mcmc.parallel(cc=i,
#                                                    obs_A = obs_A, focus = obs_F, occur_B = obs_OM, occur_P = obs_OP,
#                                                    obs_X = obs_X, obs_W = obs_W, Cu = Cu, Cv = Cv,
#                                                    Nsims = Nsims, burn = burn, thin = thin,
#                                                    use_H = use_H, bias_cor = bias_cor,use_shrinkage = TRUE,
#                                                    theta_inf = theta_inf, mh_n_pis = mh_n_pis,
#                                                    mh_n_pjs = mh_n_pjs, mh_n_rho = mh_n_rho,
#                                                    stick_alpha = stick_alpha, prior_theta = prior_theta,
#                                                    prior_tau = prior_tau, prior_rho = prior_rho,
#                                                    prior_mu0 = prior_mu0, prior_sigmasq0 = prior_sigmasq0,
#                                                    prior_sigmasq = prior_sigmasq, start_values = start_values,
#                                                    sampling = sampling, cut_feed = FALSE),
#                          1:3))


#---------------------- STEP 3: RUN THE SAMPLER -------------------------------------------#
t1 <- Sys.time()
mclapply(1:n.chains, function(i) mcmc.parallel(cc=i,
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
                           mc.cores = 3)



Sys.time() - t1
detectCores()
