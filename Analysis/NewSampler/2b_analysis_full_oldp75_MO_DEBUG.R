# LOOKING AT STRANGE RESULTS OF PLANT OCCURRENCE INDICATORS
# This version doesn't save full probabilities but rather running mean
# Doesn't save marginal probabilities or latent factors at all
# Saves a running mean of the occurrence indicators

# --------- TO DO: set  your directories and name the current results files using the date ---------#

## YOU REALLY HAVE TO SET THESE ##

# The directory where the analysis is performed:
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"
setwd(wd_path)

# Save results using convention: res_date_i.rda
date <- 'Final_old_p75mod_MO_DEBUG'

## THESE SHOULD BE THE SAME IF YOU CLONED THE REPO ##

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
save_path <- 'Results/NewSampler/'
# Where the functions are available:
source_path <- 'HelperScriptsNew/'

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
source(paste0(source_path, 'MCMC_function_trim.R'))
#source(paste0(source_path, 'MCMC_function_trimmore.R'))
source(paste0(source_path, 'PredictInteractions_function.R'))
source(paste0(source_path, 'GetPredLatFac_function.R'))
source(paste0(source_path, 'GetPredWeights_function.R'))

library(parallel)
library(foreach)
library(abind)
library(magrittr)
library(tidyverse)


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
obs_X <- Obs_X
obs_W <- Obs_W

# Try simplifying traits and getting rid of NAs
obs_W <- obs_W[, 1:2]

## Subset to remove baboons and baboon-only plants
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

bias_cor <- TRUE # Performing bias correction.

Nsims <-  100 #10000 #500 #1000 #5000 # original 10000, reasonable 2500
burn <-  100 # 40000 #2000 #5000 # 2000 # 22000 # original 40000, reasonable 2500
thin <-  1 #10 #5 # original 40
use_H <- 10 # original 10
theta_inf <- 0.01
mh_n_pis <- 70  # Parameter for proposal in Metropolis-Hastings for pi update.
mh_n_pjs <- 70
mh_n_rho <- 100
mh_occ_sd <- 0.1
mh_occ_step <- 0.25

# Prior distributions:
stick_alpha <- 5
prior_theta <- c(1, 1)
prior_tau <- c(5, 5)
prior_rho <- c(5, 5)  # I do not update this for now.
prior_mu0 <- 0
prior_sigmasq0 <- 10
prior_sigmasq <- c(1, 1)

# sampling <- NULL
sampling <- list(L = TRUE, lambda = TRUE, tau = TRUE, beta = TRUE,
                 gamma = TRUE, sigmasq = TRUE, sigmasq_p = TRUE,
                 delta = TRUE, zeta = TRUE, U = TRUE, V = TRUE, v = TRUE,
                 z = TRUE, theta = TRUE, pis = TRUE, pjs = TRUE, rU = TRUE,
                 rV = TRUE, miss_X = TRUE, miss_W = TRUE, O_B = FALSE,
                 O_P = TRUE)
start_values <- NULL


# --------------- STEP 2: SETUP PARALLEL ----------------- #

# # Set up function to execute in parallel
n.chains <- 1
mcmc.parallel <- function(cc, obs_A, focus, occur_B, occur_P, obs_X, obs_W, Cu, Cv,
                               Nsims, burn, thin, use_H = 10, use_shrinkage = TRUE,
                               bias_cor = TRUE, theta_inf = 0.01,
                               mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100, 
                               mh_occ_sd = 0.1, mh_occ_step,
                               stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
                               prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
                               prior_sigmasq = c(1, 1), start_values = NULL,
                               sampling = NULL, cut_feed = FALSE, trim = FALSE){


set.seed(cc)

  ## Create a unique filename for each interation of the parallel loop
  each_filename <- paste0('res_', date, '_', as.character(cc), '.dat')
  each_filepath <- file.path(save_path, each_filename)
  
  mcmc <- MCMC.trim(obs_A, focus, occur_B, occur_P, obs_X, obs_W, Cu, Cv,
           Nsims, burn, thin, use_H = 10, use_shrinkage = TRUE,
           bias_cor = TRUE, theta_inf = 0.01,
           mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100,
           stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
           prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
           prior_sigmasq = c(1, 1), start_values = NULL,
           sampling = sampling, cut_feed = FALSE, trim = FALSE) 
  
  
  # Binding different predictions of interest: Posterior samples of the
  # interaction indicators, the linear predictor of the interaction model,
  all_pred <- abind::abind(pred_L = mcmc$Ls, probL = mcmc$mod_pL1s, along = 4)
  
  # Phylogenetic correlation parameter for bird and plant correlation matrices.
  correlations <- cbind(U = mcmc$rU, V = mcmc$rV)
  
  # Running mean of detection probabilities
  p_detect <- list(pis = mcmc$pi_mean, pjs = mcmc$pj_mean)
  
  # Running mean of latent factors
  factors <- list(U = mcmc$U_mean, V = mcmc$V_mean)
  
  # Occurrence probabilities for plants
  occ_plants <- list(mcmc$OP_mean, mcmc$OP_samples)
  
  # Combining the results we are interested in to a list and saving:
  res <- list(all_pred = all_pred, correlations = correlations, 
              p_detect = p_detect, factors = factors, occ_plants = occ_plants)
  save(res, file = each_filepath)
  
  # 
  # plot_pred <- apply(mcmc$Ls, c(2,3), mean)
  # 
  # # The combined network:
  # comb_A <- apply(obs_A, c(1, 2), sum)
  # comb_A <- (comb_A > 0) * 1 # 1 if any interactions observed
  # plot_pred[comb_A == 1] <- NA
  # superheat(plot_pred)
  

  rm(res)
}


#---------------------- STEP 3: RUN THE SAMPLER -------------------------------------------#
detectCores()
t1 <- Sys.time()
mclapply(1:n.chains, function(i) mcmc.parallel(cc=i,
                                                          obs_A = obs_A, focus = obs_F, occur_B = obs_OM, occur_P = obs_OP,
                                                          obs_X = obs_X, obs_W = obs_W, Cu = Cu, Cv = Cv,
                                                          Nsims = Nsims, burn = burn, thin = thin,
                                                          use_H = use_H, bias_cor = bias_cor,use_shrinkage = TRUE,
                                                          theta_inf = theta_inf, mh_n_pis = mh_n_pis,
                                                          mh_n_pjs = mh_n_pjs, mh_n_rho = mh_n_rho,
                                                          mh_occ_sd = mh_occ_sd, mh_occ_step = mh_occ_step, 
                                                          stick_alpha = stick_alpha, prior_theta = prior_theta,
                                                          prior_tau = prior_tau, prior_rho = prior_rho,
                                                          prior_mu0 = prior_mu0, prior_sigmasq0 = prior_sigmasq0,
                                                          prior_sigmasq = prior_sigmasq, start_values = start_values,
                                                          sampling = sampling, cut_feed = FALSE, trim = FALSE),
                           mc.cores = n.chains)



Sys.time() - t1


#----------------------------- LOOK AT THE RESULTS ----------------------------#
# Define a useful function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
res <- loadRData(paste0(save_path, 'res_',  date, '_', 1, '.dat')) 

# Extract results
occ_plants <- res$occ_plants
mean_OP <- occ_plants[[1]] # running mean returned by the function
OP_samples <- occ_plants[[2]]
OP_post_mean <- apply(OP_samples, c(2,3), mean) # should be the same as mean above
sum(OP_post_mean != mean_OP)

sum(OP_samples[1,,] != OP_samples[2,,])

par(mfrow = c(1,2))
#png(file = paste0(results_path, "Figures/OccPlants_", date, ".png"))
hist(mean_OP, main = "Plants occurrence indictors", xlab = "Posterior Mean")
abline(v = mean(mean_OP))
hist(obs_OP, main = "Plants occurrence indictors", xlab = "Prior Mean")
abline(v = mean(obs_OP))
#dev.off()


# Change filtering depending on which version using
df2 <- data.frame(obs_op = c(ifelse(obs_OP ==1, 1, 0.75)), post_op = c(mean_OP),
                  color = factor(obs_OP, labels = c("Same region only", "Same country only",
                                                    "Same habitat only", "Same habitat and region", 
                                                    "Same habitat and country", "Same site", "Same study")), 
                  type = "Occurrence indicators") %>% 
  filter(obs_op == 0.75)


df2 %>% group_by(color) %>% 
  mutate(weight = 1/n()) %>% 
  ggplot(aes(x =post_op)) + 
  geom_histogram() +
  ylab("") + 
  xlab("Posteior mean") + 
  guides(fill = "none")+ 
  facet_wrap(~color) + 
  theme_minimal() + 
  # scale_x_continuous(limits = c(0.65,0.8)) + 
  theme(text = element_text(family = "serif", size = 20))

df2 %>% 
  group_by(color, post_op) %>%
  summarize(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% 
  ggplot(aes(x= post_op, y = freq)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~ color)

table(df2$post_op, df2$color)

# Look at values across chains for a few plants x studies
#png(paste0(results_path, date, "_mixing.png"), width = 2000, height =1000)
par(mfrow = c(3,4))
for(p in 1:12){
  plant <- sample(1:nP,1)
  study <- sample(1:nS,1)
  plot(OP_samples[,plant,study], main = paste0(plant, " in Study ", study))
}
#dev.off()
par(mfrow = c(1,1))

## Findings

# Indicators are constant, do not change at each iteration
# We were seeing 0, 1/4, 1/2, 3/4, 1 previously because taking average over 4 chains!
# Differences are due entirely to initialization

## To do: 
# X Why is the range of the posterior so much wider than with the new sampler?
# X And why is it so similar for all of the co-occurrence categories?
# X Why are the posterior means all in (0, 1/4, 1/2, 3/4, 1)?
# X Go into the function, save all samples of obs_OP
# Why does it never update from the initialized value
# Still need to figure this out, but also consider blocking indicator and prob in new sampler
# Is it possible I'm saving the wrong info?

## Go into the function and see why it isn't updating
focus = obs_F
occur_B = obs_OM
occur_P = obs_OP
use_shrinkage = TRUE
sampling = sampling
cut_feed = FALSE
trim = FALSE
Nsims <-  10 #10000 #500 #1000 #5000 # original 10000, reasonable 2500
burn <-  10 # 40000 #2000 #5000 # 2000 # 22000 # original 40000, reasonable 2500
thin <-  1 

# Findings: it updates correctly within the function and after setting assignments, thereafter
# But when I run it from a clear workspace, it doesn't work
# So I'm guessing one of my variables is not set correctly in the parallel mcmc function
# Now when I go back into the funciton it doesn't work!
# Definitely something fishy going on
