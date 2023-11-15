# This file performs a grid search for reasonable occurrence probabilities
# Note: results isn't saved properly, see 5_loadres_gridsearch.R for code to create

# --------- TO DO: set  your directories and name the current results files using the date---------#

## YOU REALLY HAVE TO SET THESE ##

# The directory where the analysis is performed:
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"
setwd(wd_path)

# Save results using convention: res_date_i.rda
date <- '6_6_gridsearch_'

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
library(pheatmap)


# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv_phylo.dat'))
load(paste0(data_path, 'obs_A_mammals.dat'))
load(paste0(data_path, 'F_obs_default.dat'))
load(paste0(data_path, 'Obs_X.dat')) # mammal traits
load(paste0(data_path, 'Obs_W.dat')) # plant traits
load(paste0(data_path, 'obs_OM.dat')) # site level obs_OM
load(paste0(data_path, 'obs_OP.dat')) # site level obs_OP
load(paste0(data_path, 'Subset_no_baboons.dat')) # list with appropriately subset mammals, plants, studies

## Rename for convenience
Cu <- Cu_phylo
Cv <- Cv_phylo
obs_A <- A.obs.m
obs_W <- Obs_W[, 1:2]
obs_X <- Obs_X

## Subset out baboons and associated plants, studies
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

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

## The following assignments were used in creating obs_OP
# Same study: 1, same site: 0.75
# Same country and habitat: 0.5, same region and habitat: 0.45, same habitat only: 0.25,
# Same country not habitat: 0.1, same region not habitat: 0.05

# Create list of probabilities to iterate over
p1.site <- c(0.5, 0.7, 0.9)
p2.hab.c <- seq(0.2, 0.8, by = 0.2)
p3.hab.r <- seq(0.1, 0.4, by = 0.2)
p4.h <- seq(0.2, 0.4, by = 0.2)
p5.c <- c(0.01, 0.05, 0.1)
p6.r <- c(0.01, 0.05, 0.1)

prob.grid <- c(1, 0.75, 0.5, 0.45, 0.25, 0.1, 0.05)
# First exercise, just create nested loops that create a grid of probabilities
for(p1 in p1.site){
  for(p2 in p2.hab.c[p2.hab.c<p1]){
    for(p3 in p3.hab.r[p3.hab.r< p2]){
      for(p4 in p4.h[p4.h<p3]){
        for(p5 in p5.c[p5.c<p4]){
          for(p6 in p6.r[p6.r<p5]){
            prob.grid <- rbind(prob.grid, c(1,p1,p2,p3,p4,p5,p6))
            
          }
        }
      }
    }
    
  }
}


# -------------- STEP 1: Specifications. ------------ #

bias_cor <- TRUE # Performing bias correction.

repetitions <- 4
Nsims <- 500 #500 #1000 #5000 # original 10000, reasonable 2500
burn <- 1000 #2000 #5000 # 2000 # 22000 # original 40000, reasonable 2500
thin <- 10  #10 #5 # original 40
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



# ------------- STEP 2: Parallel setup --------------- #

# # Debugging assignments
# occur_B <- obs_OM
# occur_P <- obs_OP
# focus <- obs_F
# prob.grid <- rbind(c(1,0,0,0,0,0,0), c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))

# Define the cv function to run in parallel 

mcmc.cv.parallel <- function(x, n.cv, prob.grid, repetitions,
                             obs_A, focus, occur_B, occur_P, obs_X, obs_W, Cu, Cv,
                             Nsims, burn, thin, use_H = 10, use_shrinkage = TRUE,
                             bias_cor = TRUE, theta_inf = 0.01,
                             mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100,
                             stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
                             prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
                             prior_sigmasq = c(1, 1), start_values = NULL,
                             sampling = NULL, cut_feed = FALSE){
  
  ### Set up the plant occurrence matrix using the current probability
  # We parallelize this part because number of probs > number of repetitions
  p <- prob.grid[x,]
  
  this_obs_OP <- occur_P
  this_obs_OP <- ifelse(this_obs_OP == 0.75, p[2], this_obs_OP)
  this_obs_OP <- ifelse(this_obs_OP == 0.5, p[3], this_obs_OP)
  this_obs_OP <- ifelse(this_obs_OP == 0.45, p[4], this_obs_OP)
  this_obs_OP <- ifelse(this_obs_OP == 0.25, p[5], this_obs_OP)
  this_obs_OP <- ifelse(this_obs_OP == 0.1, p[6], this_obs_OP)
  this_obs_OP <- ifelse(this_obs_OP == 0.05, p[7], this_obs_OP)
  
  ### Set up storage for CV results
  # Create an array to store posterior means from each cv rep
  nM <- nrow(obs_A)
  nP <- ncol(obs_A)
  our_preds <- array(NA, dim = c(nM, nP, repetitions))
  
  # Create
  perf_metrics <- matrix(NA, nrow = repetitions, ncol = 5)
  
  ### Run the CV: this is not parallelized
  
  for (rr in 1  : repetitions) {
    
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
  use_F <- focus
  for (ii in 1 : n.cv) {
    use_F[cv_indices[ii, 1], cv_indices[ii, 2], ] <- 0
  }

  
  # Running the MCMC with the new recorded interaction matrix:
  mcmc <- MCMC.trim(obs_A = use_A, focus = use_F, occur_B = obs_OM, occur_P = this_obs_OP,
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
  
  
  ### Saving things and creating the results to return
  
  # Saving the predictions:
  pred <- apply(mcmc$Ls, c(2, 3), mean) # Why not save the interaction probs instead of the predictions?
  our_preds[, , rr] <- pred # Save array with all the cv reptitions for this probability set
  cv_preds <- pred[cv_indices]
  
  save(pred, file = paste0(save_path, 'GridSearch/pred_probs_', x, '_', rr, '.dat'))
  save(cv_indices, file = paste0(save_path, 'GridSearch/cv_indices_probs_', x, '_', rr, '.dat'))
  
  # Add performance metrics for this repetition
  overall_mean <- mean(pred)
  overall_median <- median(pred)
  
  # Average and median in the held out data.
  pred_mean <- mean(cv_preds)
  pred_median <- median(cv_preds)
  
  # Performance metrics
  mean_ratio <- mean(pred_mean/overall_mean)
  median_ratio <- mean(pred_median/overall_median)
  cv_mean <- mean(pred_mean)
  p50 <- sum(pred>0.5)/length(pred) # what proportion of true interactions are predicted as "likely" >0.5
  p75 <- sum(pred>0.75)/length(pred) # what proportion of true interactions are predicted as "v likely" >0.75
  
  res <- c(mean_ratio, median_ratio, cv_mean, p50, p75)
  perf_metrics[rr, ]<- res
  
  rm(pred)
  rm(cv_indices)
  }
  
  mean_pred <- apply(our_preds, c(1,2), mean)
  png(filename = paste0(save_path, "GridSearch/heatmap_p",x, ".png"))
  heatmap(mean_pred, Rowv = NA, Colv = NA, scale = 'none')
  dev.off()
  
  # Return a vector of performance metrics averaged across cv repetitions for this run
  mean_perf <- colMeans(perf_metrics)
  res <- c(x, p, mean_perf) # Include the current prob just in case they get jumbled
}

#---------------------- STEP 3: RUN THE SAMPLER -------------------------------------------#

# We highly recommend running the following code in parallel on 30 machines.
repetitions <- 4
n.cv <- 100

t1 <- Sys.time()
res <- mclapply(1:nrow(prob.grid), function(i) mcmc.cv.parallel(x=i, n.cv = n.cv, prob.grid = prob.grid,  repetitions = repetitions,
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
         mc.cores = nrow(prob.grid))

res.df <- do.call(rbind, res)
colnames(res.df) <- c("p.set", paste0("p", 1:7), "mean.ratio", "median.ratio", "cv.mean", "cvp50", "cvp75")

write.csv(res.df, paste0(save_path, "GridSearch/results_summary.csv"), row.names = FALSE)


Sys.time() - t1