# Performing cross-validation using the proposed bird-plant interaction model.

# The directory where the analysis is performed:
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"
# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
result_path <- 'Results/'
# Where the functions are available:
source_path <- 'HelperScriptsNew/'

#Save results using convention: res_date_i.rda
#date <- 'April19_cutFeedbackFull'
date <- 'May1_Occ1mod'

setwd(wd_path)
# ------ STEP 0: Some functions. --------- #

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
library(abind)

# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv_phylo.dat'))
load(paste0(data_path, 'obs_A_mammals.dat'))
load(paste0(data_path, 'F_obs_default.dat'))
load(paste0(data_path, 'Obs_X.dat')) # mammal traits
load(paste0(data_path, 'Obs_W.dat')) # plant traits
#load(paste0(data_path, 'obs_OP_defaultMethod.dat'))
#load(paste0(data_path, 'obs_OM_defaultMethod.dat'))
load(paste0(data_path, 'obs_OM.dat')) # mammal occurrences: many species have only one study
load(paste0(data_path, 'obs_OP.dat')) # plant occurrences: most species have many studies, but some have 0?
load(paste0(data_path, 'traits_p_709_clean.dat')) # plant occurrences: most species have many studies, but some have 0?

## Rename for convenience
Cu <- Cu_phylo
Cv <- Cv_phylo
obs_A <- A.obs.m
# obs_OM <- obs_OM2
# obs_OP <- obs_OP2

## Subset plants to simplify for exploratory analysis
apply(obs_A, 3, sum) #Study.16 aka index 14 has only one interaction, drop it and drop the monkey
good.studies <- which(apply(obs_A, 3, sum)>1)
common.plants <- names(which(apply(obs_A[,,good.studies], 2, sum)>5))
good.mammals <- names(which(apply(obs_A[,,good.studies], 1, sum)>1))

obs_A_common <- obs_A[unlist(dimnames(obs_A)[1]) %in% good.mammals, unlist(dimnames(obs_A)[2]) %in% common.plants, good.studies ]
obs_OP_common <- obs_OP[rownames(obs_OP) %in% common.plants, ]
obs_OM_common <- obs_OM[rownames(obs_OM) %in% good.mammals,]
F_obs_common <- obs_F[rownames(obs_OM) %in% good.mammals, colnames(obs_F) %in% common.plants, good.studies]
Cu_common <- Cu[rownames(Cu) %in% good.mammals, colnames(Cu) %in% good.mammals]
Cv_common <- Cv[rownames(Cv) %in% common.plants, colnames(Cv) %in% common.plants]
obs_W_common <- Obs_W[rownames(Obs_W) %in% common.plants, ]
obs_X_common <- Obs_X[rownames(Obs_X) %in% good.mammals, ]

# Try simplifying traits and getting rid of NAs
obs_X_common <- obs_X_common[, 2:3]
obs_W_common <- obs_W_common[, 1:2]

# No missing data in the trait matrices
colSums(obs_X_common)
colSums(obs_W_common)

# Sample sizes of the two sets of species:
nM <- dim(obs_A_common)[1]
nP <- dim(obs_A_common)[2]
nS <- dim(obs_A_common)[3]
print(c(nM, nP, nS))

# Try replacing occurrence data with matrix of 1s
# obs_OP_common <- matrix(1, nrow = nP, ncol = nS)
# obs_OM_common <- matrix(1, nrow = nM, ncol = nS)

# Try replacing 0s in site-level occurrence with 0.5
obs_OP_common <- ifelse(obs_OP_common ==0, 0.5, 1)
obs_OM_common <- ifelse(obs_OP_common ==0, 0.5, 1)

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A_common, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1


# -------------- STEP 1: Specifications. ------------ #

bias_cor <- TRUE  # Performing bias correction.

Nsims <- 5000 #5000 # original 10000, reasonable 2500
burn <- 2000 # 22000 # original 40000, reasonable 2500
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

## We fix
# O_B: assume mammals are all detected because studies are animal-focused
# Since our O_B, O_P are binary, they are by default not updated
# Shrinkage parameters, theta, tau = 1


## Also try: fixing more detection parameters, rU phylo covar for animals

sampling <- list(L = TRUE, lambda = TRUE, tau = FALSE, beta = TRUE,
                 gamma = TRUE, sigmasq = TRUE, sigmasq_p = TRUE,
                 delta = TRUE, zeta = TRUE, U = TRUE, V = TRUE, v = TRUE,
                 z = TRUE, theta = FALSE, pis = FALSE, pjs = TRUE, rU = TRUE, 
                 rV = TRUE, miss_X = TRUE, miss_W = TRUE, O_B = TRUE, # don't sample animal occurrences
                 O_P = TRUE)

start_values <- list(theta = 1, tau = 1)

sampling <- NULL
start_values <- NULL


# ------------- STEP 2: Setting some interactions to out of sample ------------ #

# We highly recommend running the following code in parallel on 30 machines.

repetitions <- 10
n.cv <- 35


for (rr in 1  : repetitions) {
  
  set.seed(rr)
  
  # Matrix that chooses n.cv recorded interactions to be held-out.
  set_out <- matrix(0, nrow = nM, ncol = nP)
  set_out[sample(which(comb_A == 1), n.cv)] <- 1
  
  # # Zero-ing out corresponding entries in A. Seems redundant with later code
  # use_A <- obs_A_common
  # use_A[which(set_out == 1)] <- 0  

  
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
  use_A <- obs_A_common
  for (ii in 1 : n.cv) {
    use_A[cv_indices[ii, 1], cv_indices[ii, 2], ] <- 0
  }
  
  # Shouldn't we also set Focus to zero for these observations?
  use_F <- F_obs_common
  for (ii in 1 : n.cv) {
    use_F[cv_indices[ii, 1], cv_indices[ii, 2], ] <- 0
  }
  
  # Running the MCMC with the new recorded interaction matrix:
  set.seed(rr)
  mcmc <- MCMC(obs_A = use_A, focus = use_F, occur_B = obs_OM_common, occur_P = obs_OP_common,
               obs_X = obs_X_common, obs_W = obs_W_common, Cu = Cu_common, Cv = Cv_common,
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
  
  save(cv_indices, file = paste0(result_path, 'cv_indices_', date, '_', rr, '.dat'))
  rm(cv_indices)
}
  

