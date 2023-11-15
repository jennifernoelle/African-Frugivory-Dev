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
#date <- 'April19_cutFeedbackTrimmed'
date <- 'May9_OccExpert1'
batchname <- "mcmc_4_13"


## THESE SHOULD BE THE SAME IF YOU CLONED THE REPO ##

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
save_path <- 'Results/'
# Where the functions are available:
source_path <- 'HelperScriptsNew/'

#setwd(wd_path)

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
source(paste0(source_path, 'MCMC_function_trimResults.R'))
source(paste0(source_path, 'MCMC_function.R'))
source(paste0(source_path, 'PredictInteractions_function.R'))
source(paste0(source_path, 'GetPredLatFac_function.R'))
source(paste0(source_path, 'GetPredWeights_function.R'))

library(doParallel)
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
load(paste0(data_path, 'subset1_studies.dat')) 
load(paste0(data_path, 'subset1_mammals.dat')) 
load(paste0(data_path, 'subset1_plants.dat')) 

## Rename for convenience
Cu <- Cu_phylo
Cv <- Cv_phylo
obs_A <- A.obs.m

wh_keep_m <- which(rownames(obs_A) %in% good.mammals)
wh_keep_p <- which(colnames(obs_A) %in% good.plants)
wh_keep_s <- which(unlist(dimnames(obs_A)[3]) %in% good.studies)

obs_A <- obs_A[wh_keep_m, wh_keep_p, wh_keep_s]
F_obs <- obs_F[wh_keep_m, wh_keep_p, wh_keep_s]
obs_OP <- obs_OP[wh_keep_p, wh_keep_s]
obs_OM <- obs_OM[wh_keep_m, wh_keep_s]
obs_X <- Obs_X[wh_keep_m, ]
obs_W <- Obs_W[wh_keep_p, ]
Cu <- Cu[wh_keep_m, wh_keep_m]
Cv <- Cv[wh_keep_p, wh_keep_p]

# Try simplifying traits and getting rid of NAs
# obs_X <- obs_X[, 2:3] # No missing values, try more covars
obs_W <- obs_W[, 1:2] # Dropping IUCN status because high  missingness

# Useful values
nP <- ncol(obs_A)
nM <- nrow(obs_A)
nS <- length(good.studies)

# Create alternate expert-defined occurrence
# obs_OP <- ifelse(obs_OP == 0.75, 0.85, obs_OP) # boost same site P(occ) = 0.85 #2
# obs_OP <- ifelse(obs_OP == 0.5, 0.65, obs_OP) # boost country-habitat P(occ) = 0.65 #2
# obs_OP <- ifelse(obs_OP == 0.25, 0.45, obs_OP) # boost same habitat P(occ) = 0.45 #3 is 0.35 for #4
# obs_OP <- ifelse(obs_OP == 0, 0.1, obs_OP) #4

# -------------- STEP 1: Specifications. ------------ #

bias_cor <- TRUE # Performing bias correction.

Nsims <- 5000 # original 10000, reasonable 2500
burn <-  2000 # 22000 # original 40000, reasonable 2500
thin <- 5 # original 40
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

# --------------- STEP 2: MCMC. ----------------- #

# We run 4 chains. We suggest that you run the following code in parallel instead.

nchains <- 1

t1 <- Sys.time()

for (cc in 1 : nchains) {  # Chain index:
  
  set.seed(cc)
  
  
  mcmc <- MCMC(obs_A = obs_A, focus = F_obs, occur_B = obs_OM, occur_P = obs_OP,
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

  # Attaching the results:
  #attach(mcmc)
  
  # Binding different predictions of interest: Posterior samples of the
  # interaction indicators, the linear predictor of the interaction model,
  # and the probability we use when sampling the interaction indicators.
  # Studying MCMC() will clarify the three quantities.
  all_pred <- abind::abind(pred_L = mcmc$Ls, probL = mcmc$mod_pL1s, pL1s = mcmc$pL1s, along = 4)
  
  # Phylogenetic correlation parameter for bird and plant correlation matrices.
  correlations <- cbind(U = mcmc$rU, V = mcmc$rV)
  
  p_detect <- list(pis = mcmc$pis, pjs = mcmc$pjs)
  
  # Combining the results we are interested in to a list and saving:
  #res <- list(pL1s = all_pred, correlations = correlations)
  res <- list(all_pred = all_pred, correlations = correlations, p_detect = p_detect)
  save(res, file = paste0(save_path, 'res_', date, "_", cc, '.dat'))
  
  rm(res)
  #detach(mcmc)

}
Sys.time() - t1

time_elapsed <- print(Sys.time() - t1)


# For line by line debugging
# obs_A = obs_A
# focus = F_obs
# occur_B = obs_OM
# occur_P = obs_OP
# obs_X = obs_X
# obs_W = obs_W
# Cu = Cu
# Cv = Cv
# Nsims = Nsims
# burn = burn
# thin = thin
# use_H = use_H
# bias_cor = bias_cor
# use_shrinkage = TRUE
# cut_feed = FALSE

