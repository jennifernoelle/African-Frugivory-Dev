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
date <- 'April3_shrink_and_feedback'
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
load(paste0(data_path, 'obs_OM.dat')) # mammal occurrences: many species have only one study
load(paste0(data_path, 'obs_OP.dat')) # plant occurrences: most species have many studies, but some have 0?
load(paste0(data_path, 'traits_p_709_clean.dat')) # plant occurrences: most species have many studies, but some have 0?

## Rename for convenience
Cu <- Cu_phylo
Cv <- Cv_phylo
obs_A <- A.obs.m

## Useful quantities: sample sizes of the two sets of species
nM <- nrow(Cu)
nP <- nrow(Cv)
nStudies <- dim(obs_A)[3]

## Sanity checks

# Construct the matrix F_ijs * O_ijs
O_ijs <- array(data = NA, dim = c(nM, nP, nStudies))
dimnames(O_ijs) <- list(rownames(obs_A), colnames(obs_A), paste0("Study_", 1:nStudies))
for(s in 1:nStudies){
  for(m in 1:nM){
    for(p in 1:nP){
      O_ijs[m,p,s] <- obs_OM[m,s] * obs_OP[p,s] 
    }
  }
}

# Interaction prevalence is very low, even conditioning on observability (O*F)
obs_A_subset <- obs_A
obs_A_subset[F_obs_m ==0] <- NA
obs_A_subset[O_ijs==0] <- NA
mean(obs_A_subset, na.rm = TRUE)
mean(obs_A)

# All species occur at least once
summary(rowSums(obs_OP))
summary(rowSums(obs_OM))

# Look at plant and mammal phylogenetic correlations
summary(c(Cu))
summary(c(Cv))

## Subset plants to simplify for exploratory analysis
# apply(obs_A, 2, sum) # note sum plants have no observed interactions but they were present, this is fine
# common.plants <- names(which(apply(obs_A, 2, sum)>5))
# 
# obs_A_common <- obs_A[, unlist(dimnames(obs_A)[2]) %in% common.plants, ]
# obs_OP_common <- obs_OP[rownames(obs_OP) %in% common.plants, ]
# F_obs_common <- F_obs_m[, colnames(F_obs_m) %in% common.plants, ]
# Cv_common <- Cv[rownames(Cv) %in% common.plants, colnames(Cv) %in% common.plants]
# obs_W_common <- Obs_W[rownames(Obs_W) %in% common.plants, ]
apply(obs_A, 3, sum) #Study.16 aka index 14 has only one interaction, drop it and drop the monkey
good.studies <- which(apply(obs_A, 3, sum)>1)
common.plants <- names(which(apply(obs_A[,,good.studies], 2, sum)>5))
good.mammals <- names(which(apply(obs_A[,,good.studies], 1, sum)>1))

obs_A_common <- obs_A[unlist(dimnames(obs_A)[1]) %in% good.mammals, unlist(dimnames(obs_A)[2]) %in% common.plants, good.studies ]
obs_OP_common <- obs_OP[rownames(obs_OP) %in% common.plants, ]
obs_OM_common <- obs_OM[rownames(obs_OM) %in% good.mammals,]
F_obs_common <- F_obs_m[rownames(obs_OM) %in% good.mammals, colnames(F_obs_m) %in% common.plants, good.studies]
Cu_common <- Cu[rownames(Cu) %in% good.mammals, colnames(Cu) %in% good.mammals]
Cv_common <- Cv[rownames(Cv) %in% common.plants, colnames(Cv) %in% common.plants]
obs_W_common <- Obs_W[rownames(Obs_W) %in% common.plants, ]
obs_X_common <- Obs_X[rownames(Obs_X) %in% good.mammals, ]

# Try simplifying traits and getting rid of NAs
obs_X_common <- obs_X_common[, 2:3]
obs_W_common <- obs_W_common[, 1:2]

O_ijs_common <- O_ijs[rownames(O_ijs) %in% good.mammals, colnames(O_ijs) %in% common.plants, good.studies]


### Make sure the data looks as expected

# Check ordering of species
sum(rownames(obs_A_common) != rownames(F_obs_common))
sum(colnames(obs_A_common) != colnames(F_obs_common))
sum(rownames(obs_X_common) != rownames(obs_A_common))
sum(rownames(obs_W_common) != colnames(obs_A_common))
sum(rownames(obs_OM_common) != rownames(obs_A_common))
sum(rownames(obs_OP_common) != colnames(obs_A_common))

# Overall interaction level is low
mean(obs_A_common)
mean(obs_A_common[F_obs_common == 1])
mean(obs_A_common[O_ijs_common*F_obs_common == 1])

# No missing data in the trait matrices
colSums(obs_X_common)
colSums(obs_W_common)


# -------------- STEP 1: Specifications. ------------ #

bias_cor <- TRUE # Performing bias correction.

Nsims <- 10000 #5000 # original 10000, reasonable 2500
burn <- 5000 # 22000 # original 40000, reasonable 2500
thin <- 5 # original 40
use_H <- 2 # original 10
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
# Continuous and binary covariates.
entries_X <- apply(obs_X_common, 2, function(x) length(unique(x[!is.na(x)])))
pM <- c(sum(entries_X > 2), sum(entries_X == 2))
entries_W <- apply(obs_W_common, 2, function(x) length(unique(x[!is.na(x)])))
pP <- c(sum(entries_W > 2), sum(entries_W == 2))

sampling <- list(L = TRUE, lambda = TRUE, tau = FALSE, beta = TRUE,
                 gamma = TRUE, sigmasq = FALSE, sigmasq_p = FALSE,
                 delta = TRUE, zeta = TRUE, U = TRUE, V = TRUE, v = TRUE,
                 z = TRUE, theta = FALSE, pis = TRUE, pjs = TRUE, rU = FALSE, 
                 rV = FALSE, miss_X = TRUE, miss_W = TRUE, O_B = TRUE, 
                 O_P = TRUE)

start_values <- list(this_theta = rep(1, use_H), this_tau_zeta = 1, this_tau_beta = matrix(1, sum(pM), use_H), 
                     this_tau_gamma = matrix(1, sum(pP), use_H), 
                     this_tau_delta = 1, this_tau_lambda = 1, this_ru = 0, this_rv = 0, 
                     this_sigmasqa_m = rep(1, pM[1]), this_sigmasq_l = rep(1, pP[1]), 
                     this_sigmasq_pB = 1, this_sigmasq_pP= 1)

sampling <- NULL
start_values <- NULL

# --------------- STEP 2: MCMC. ----------------- #

# We run 4 chains. We suggest that you run the following code in parallel instead.

nchains <- 2

#registerDoParallel(4)
#foreach(cc=1:nchains) %dopar% {
t1 <- Sys.time()

for (cc in 1 : nchains) {  # Chain index:
  
  set.seed(cc)
  
  # # Running the method (trimmed version):
  # mcmc <- MCMC.trim(obs_A = obs_A_common, focus = F_obs_common, occur_B = obs_OM, occur_P = obs_OP_common,
  #              obs_X = Obs_X, obs_W = obs_W_common, Cu = Cu, Cv = Cv_common,
  #              Nsims = Nsims, burn = burn, thin = thin,
  #              use_H = use_H, bias_cor = bias_cor,use_shrinkage = TRUE,
  #              theta_inf = theta_inf, mh_n_pis = mh_n_pis,
  #              mh_n_pjs = mh_n_pjs, mh_n_rho = mh_n_rho,
  #              stick_alpha = stick_alpha, prior_theta = prior_theta,
  #              prior_tau = prior_tau, prior_rho = prior_rho,
  #              prior_mu0 = prior_mu0, prior_sigmasq0 = prior_sigmasq0,
  #              prior_sigmasq = prior_sigmasq, start_values = start_values,
  #              sampling = sampling, cut_feed = TRUE)
  # 
  
  mcmc <- MCMC(obs_A = obs_A_common, focus = F_obs_common, occur_B = obs_OM_common, occur_P = obs_OP_common,
                             obs_X = obs_X_common, obs_W = obs_W_common, Cu = Cu_common, Cv = Cv_common,
                             Nsims = Nsims, burn = burn, thin = thin,
                             use_H = use_H, bias_cor = bias_cor,use_shrinkage = TRUE,
                             theta_inf = theta_inf, mh_n_pis = mh_n_pis,
                             mh_n_pjs = mh_n_pjs, mh_n_rho = mh_n_rho,
                             stick_alpha = stick_alpha, prior_theta = prior_theta,
                             prior_tau = prior_tau, prior_rho = prior_rho,
                             prior_mu0 = prior_mu0, prior_sigmasq0 = prior_sigmasq0,
                             prior_sigmasq = prior_sigmasq, start_values = start_values,
                             sampling = sampling, cut_feed = TRUE)

  # Attaching the results:
  #attach(mcmc)
  
  # Binding different predictions of interest: Posterior samples of the
  # interaction indicators, the linear predictor of the interaction model,
  # and the probability we use when sampling the interaction indicators.
  # Studying MCMC() will clarify the three quantities.
  all_pred <- abind::abind(pred_L = mcmc$Ls, probL = mcmc$mod_pL1s, pL1s = mcmc$pL1s, along = 4)
  #all_pred <- mcmc$pL1s_mean # trimmed version
  
  # Phylogenetic correlation parameter for bird and plant correlation matrices.
  correlations <- cbind(U = mcmc$rU, V = mcmc$rV)
  
  # Combining the results we are interested in to a list and saving:
  #res <- list(pL1s = all_pred, correlations = correlations)
  res <- list(all_pred = all_pred, correlations = correlations)
  save(res, file = paste0(save_path, 'res_', date, "_", cc, '.dat'))
  
  rm(res)
  #detach(mcmc)

}
Sys.time() - t1

time_elapsed <- print(Sys.time() - t1)

save(time_elapsed, file = paste(save_path, batchname, '_time.dat'))
# 
# # For line by line debugging
# obs_A = obs_A_common
# focus = F_obs_common
# occur_B = obs_OM
# occur_P = obs_OP_common
# obs_X = obs_X_common
# obs_W = obs_W_common
# Cu = Cu_common
# Cv = Cv_common
# Nsims = Nsims
# burn = burn
# thin = thin
# use_H = use_H
# bias_cor = bias_cor
# use_shrinkage = FALSE
# cut_feed = TRUE
