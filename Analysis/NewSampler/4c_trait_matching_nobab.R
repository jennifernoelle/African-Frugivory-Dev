# Estimating variable importance based on interaction model
# This code is slow, run overnight and consider reducing the number of permutations for exploratory analysis

# The directory where the analysis is performed:
wd_path <- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"
# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where the MCMC results are saved and the trait matching will be saved:
save_path <- 'Results/'
# Where the functions are available:
source_path <- 'HelperScriptsNew/'

date <- 'Full_nobab'

# ------ STEP 0: Some functions. --------- #

setwd(wd_path)
source(paste0(source_path, 'TraitMatching2_function.R'))
source(paste0(source_path, 'useful_functions.R'))

# --------------------------------------------------------------- #

# Loading the data:
load(paste0(data_path, 'Obs_X.dat')) # mammal traits
load(paste0(data_path, 'Obs_W.dat')) # plant traits
load(paste0(data_path, 'Subset_no_baboons.dat')) # list with appropriately subset mammals, plants, studies

## Subset out baboons and associated plants, studies
wh_keep_m <- which(rownames(obs_A) %in% nobab.list[[1]])
wh_keep_p <- which(colnames(obs_A) %in% nobab.list[[2]])

obs_X <- Obs_X[wh_keep_m, ]
obs_W <- Obs_W[wh_keep_p, ]


# Getting the sample sizes:
nM <- nrow(Obs_X)
nP <- nrow(Obs_W)

# --------------- STEP 1: Getting the results together ----------------- #

# MCMC chains saved:
nchains <- 4

# Putting together the predictions from the chains:
all_pred <- NULL
for (ii in 1 : nchains) {
  load(paste0(save_path, 'res_', date, '_', ii, '.dat'))
  all_pred[[ii]] <- res$all_pred
}

# Number of posterior samples by chain:
Nsims <- dim(all_pred[[1]])[1]

# Using the linear predictor of the interaction model (not corrected for bias due to detection, number of studies):
mod_pL1s <- array(NA, dim = c(nchains * Nsims, nM, nP)) 
for (cc in 1 : nchains) {
  wh_entries <- Nsims * (cc - 1) + 1 : Nsims
  mod_pL1s[wh_entries, , ] <- all_pred[[cc]][, , , 2]
}

# --------------- STEP 2: Performing trait matching ----------------- #

# Dealing with extreme values for which logit(x) is infinite.
mod_pL1s[mod_pL1s > 1 - 10^{-10}] <- 1 - 10^{-10}

t1 <- Sys.time()
trait_match <- TraitMatching2(B = 100, mod_pL1s = mod_pL1s,
                              Xs = NULL, Ws = NULL,  # Imputed values not used.
                              obs_X = Obs_X, obs_W = Obs_W, obs_only = TRUE)

rsq_resampling_X <- trait_match$rsq_resampling_X
rsq_resampling_W <- trait_match$rsq_resampling_W
rsq_obs_X <- trait_match$rsq_obs_X
rsq_obs_W <- trait_match$rsq_obs_W

save(rsq_obs_X, file = paste0(save_path,  'rsq_obs_X_', date, '.dat'))
save(rsq_obs_W, file = paste0(save_path, 'rsq_obs_W_', date, '.dat'))
save(rsq_resampling_X, file = paste0(save_path, 'rsq_resampling_X_', date, '.dat'))
save(rsq_resampling_W, file = paste0(save_path, 'rsq_resampling_W_', date, '.dat'))

print(Sys.time() - t1)
