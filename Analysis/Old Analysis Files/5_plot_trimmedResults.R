# This version works with the trimmed results (for memory use)

#--------------------- TO DO --------------------------#

## YOU REALLY HAVE TO SET THESE ##

# The directory where the analysis is performed:
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
wd_path <- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"


# Results were saved using convention: res_date_i.rda
date <- 'Feb23_cutFeedbackTrimmed'
batchname <- "mcmc_3_13a"


## THESE SHOULD BE THE SAME IF YOU CLONED THE REPO ##

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
results_path <- 'Results/'
# Where the functions are available:
source_path <- 'HelperScriptsPapadogeorgou/'

setwd(wd_path)


# -------------------- LOAD PACKAGES AND FILES -----------------------#
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(superheat)
library(abind)
library(gridExtra)
library(grid)
library(dplyr)
library(caret)
library(pROC)

# Loading the data:
load(paste0(data_path, 'traits_m.dat'))
load(paste0(data_path, 'traits_p_709.dat'))
load(paste0(data_path, 'obs_A_mammals.dat'))
load(paste0(data_path, 'Obs_X.dat'))
load(paste0(data_path, 'F_obs_m.dat'))
load(paste0(data_path, 'obs_OM.dat'))
load(paste0(data_path, 'obs_OP.dat'))

source(paste0(source_path, 'draw_cm.R'))

# Define a useful function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# ------------------- SETUP --------------------------#

# Rename
obs_A <- A.obs.m

# The combined network:
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1 # 1 if any interactions observed

nM <- nrow(obs_A)
nP <- ncol(obs_A)
nS <- dim(obs_A)[3]

# Currently not using trait data, fix later
#old:  Species names are consistent with order in traits data frames
plant.names <- colnames(obs_A[,,1]) #traits.p.709$Species

m.names <- rownames(obs_A[,,1])#traits.m$Species

# Number of MCMC chains for our method and for the alternative method:
nchains <- 2

# Number of cross validation repetitions:
repetitions <- 3


# --------------- STEP 1: Getting the results together ----------------- #


### Put all chains together in a list
all_res <- NULL
for (ii in 1 : nchains) {
  res <- loadRData(paste0(results_path, 'res_',  date, '_', ii, '.dat')) 
  all_res[[ii]] <- res
}


### Bind together predicted interactions

# Number of posterior samples used:
use_Nsims <- dim(all_res[[1]]$correlations)[1]

# Creating an array to bind results across chains:
mean_pred <- mean_pred_na <- Reduce("+", lapply(all_res, function(x) x[[1]]))/nchains
dimnames(mean_pred) <- dimnames(mean_pred_na) <- list(m.names, plant.names)

# Setting the recorded interactions to NA (so that they don't overpower the colors)
mean_pred_na[comb_A == 1] <- NA


#-------------------- STEP 1B: SANITY CHECK -----------------------------#

## Sanity check: is mean posterior probability similar to sample mean?
## Compare the prevalence of interactions among F*O=1 observations to predictions


# What proportion of pairs actually have observable interactions (i.e., Focus = 1 in at least one study)
# All mammals have at least one relevant mammal-focused study, so all interactions are observable
# Note also that all studies are animal-focused, except for one which only focuses on a specific pair

# Note all interactions in Focus because all studies are animal-oriented
comb_F <- apply(F_obs_m, c(1, 2), sum)
comb_F_bin <- ifelse(comb_F > 0,1,0)
mean(comb_F_bin)

# How many possible interactions for each species
F.df <- data.frame(rowSums(comb_F)) %>% 
  mutate(Species = rownames(.)) %>% 
  rename(nFocusObs = rowSums.comb_F.) 
ggplot(F.df, aes(x = Species, y = nFocusObs)) + geom_bar(stat = "identity") + 
  coord_flip() + 
  geom_hline(yintercept = nP)

## Now consider occurrence
# Construct the matrix F_ijs * O_ijs
O_ijs <- array(data = NA, dim = c(nM, nP, nS))
dimnames(O_ijs) <- list(m.names, plant.names, paste0("Study_", 1:nS))
for(s in 1:nS){
  for(m in 1:nM){
    for(p in 1:nP){
      O_ijs[m,p,s] <- obs_OM[m,s] * obs_OP[p,s] 
    }
  }
}

# Interaction prevalence where F*O = 1 <<< mean posterior prediction
obs_A_subset <- obs_A
obs_A_subset[F_obs_m ==0] <- NA
obs_A_subset[O_ijs==0] <- NA
mean(obs_A_subset, na.rm = TRUE)
mean(mean_pred)


# How many species are the focus of each study
nFocusMammals <- apply(F_obs_m, 3, function(x) sum(rowSums(x)>0))
hist(nFocusMammals, main = "Number of Focus Mammals per Study", xlab = "Number of Mammals Followed")
mean(nFocusMammals)
table(nFocusMammals)

# Accuracy check
# Other metrics: MCC, ROC curve, PR curve
# No information rate is prevalence of the largest class

cm <- confusionMatrix(factor(rbinom(length(comb_A), 1, mean_pred)), factor(c(comb_A)), 
                      positive = "1")

#png(paste0(results_path, date, "_cm.png"), width = 500, height =500)
draw_confusion_matrix(cm)
#dev.off()


# ----------- STEP 2: PLOTTING THE HEATMAP --------------------#

# Creating the clusters that will be used
# The following two lines specify that horizontal and vertical lines in our
# plot will separate species by taxonomic families:
m_group <- traits.m$Family
p_group <- traits.p.709$Family

# Calculating the size of each cluster, will be used when plotting results
# for families of certain size:
m_size_cluster <- sapply(unique(m_group), function(x) sum(m_group == x))
p_size_cluster <- sapply(unique(p_group), function(x) sum(p_group == x))

# Set plot_pred to pred_ours for results based on our method 
plot_pred <- mean_pred_na

# Set the minimum cluster size that should be plotted. For the results of the
# manuscript, we set min_bird_size to 10, and min_plant_size to 20. Setting both
# to 0 will produce the full results.
min_m_size <- 1
min_p_size <- 20

keep_m_groups <- names(which(m_size_cluster >= min_m_size))
keep_p_groups <- names(which(p_size_cluster >= min_p_size))

keep_m_index <- which(m_group %in% keep_m_groups)
keep_p_index <- which(p_group %in% keep_p_groups)

# Plotting those with minimum size as specified:
# We replaced observed interactions with black NA
#png(paste0(results_path, date, "_heatmap.png"), width = 1000, height =1000)
superheat(X = plot_pred[keep_m_index, keep_p_index],
          membership.rows = m_group[keep_m_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 3,
          bottom.label.text.size = 3,
          bottom.label.size = 0.2, left.label.size = 0.12,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          heat.col.scheme = "grey", heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05))

# Superheat fails due to all values being too close to 1
# Note: these results are strange: interactions occur to occur with much more frequency than observed
# Lots of dark grey in rows with very little black


# --------------- STEP 3: Taxonomic correlation of latent factors ----------------- #

all_cor <- abind::abind(all_res[[1]]$correlations, all_res[[2]]$correlations, along = 3)
if(nchains>2){for (cc in 3 : nchains) {
  all_cor <- abind::abind(all_cor, all_res[[cc]]$correlations, along = 3)
}
}

# Posterior means and 95% credible intervals for the rho parameters in the
# latent factors for bird and plant species:
apply(all_cor, 2, mean)
apply(all_cor, 2, quantile, probs = c(0.025, 0.975))

# Phylogenetic correlation more important for plants than animals

# --------------- STEP 4: Cross validation results ----------------- #

# Getting the results together (held out indicies and predictions)
all_indices <- array(NA, dim = c(repetitions, 100, 2))
our_preds <- array(NA, dim = c(repetitions, nM, nP))

for (rr in 1 : repetitions) {
  cv_indices <- loadRData(paste0(result_path, 'cv_indices_', date, '_', rr, '.rda'))
  pred <- loadRData(paste0(result_path, 'pred_', date, '_', rr, '.rda'))
  all_indices[rr, , ] <- cv_indices
  our_preds[rr, , ] <- pred
}

# Predictions of the held out data from model: 100 indices held out each time 
pred <- array(NA, dim = c(repetitions, 100))
for (rr in 1 : repetitions) {
  for (ii in 1 : 100) {
    pred[rr, ii] <- our_preds[rr, all_indices[rr, ii, 1], all_indices[rr, ii, 2]]
  }
}

# Average and median probability of interaction based on the overall data
overall_mean <- cbind(apply(our_preds, 1, mean))
overall_median <- cbind(apply(our_preds, 1, median))

# Average and median in the held out data.
pred_mean <- apply(pred, 1, mean)
pred_median <- apply(pred, 1, median)

# Creating the data frame we will plot:


plot_dta <- data.frame(value = rbind(pred_mean / overall_mean, pred_median / overall_median), 
                       stat = rep(c('Prediction mean / Overall mean', 'Prediction median / Overall median'), repetitions))

# Plotting cross validation results:
ggplot(data = plot_dta) +
  geom_boxplot(aes(x = stat, y = value)) +
  theme_bw() +
  ylab('') +
  xlab('') +
  ggtitle('Out of sample performance', subtitle = 'Predictions for held-out recorded interactions') +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6)

# CV results are not great: we want to see the mean for held-out recorded interactions
# is greater than the overall mean, but in fact it is lower

# One test would be to see if predictions are better if we restrict to multi-species studies


# --------------- STEP 5: Variable importance measure ----------------- #

# JK: I haven't modified this yet, might not work. 

# ------- PART A: Plotting the variables in order of importance:

load(paste0(result_path, 'rsq_obs_X.dat'))
load(paste0(result_path, 'rsq_obs_W.dat'))
load(paste0(result_path, 'rsq_resampling_X.dat'))
load(paste0(result_path, 'rsq_resampling_W.dat'))


# Calculating the number of permuted standard deviations away from the mean.

# Starting from the bird covariates:
wh_obs <- rsq_obs_X
wh_resampling <- rsq_resampling_X
sd_awayX <- rep(NA,  length(wh_obs))
names(sd_awayX) <- good_namesX
for  (cc in 1 : length(wh_obs)) {
  sd_awayX[cc] <- (wh_obs[cc] - mean(wh_resampling[, cc])) / sd(wh_resampling[, cc])
}

# And for the plant covariates:
wh_obs <- rsq_obs_W
wh_resampling <- rsq_resampling_W
sd_awayW <- rep(NA,  length(wh_obs))
names(sd_awayW) <- good_namesW
for  (cc in 1 : length(wh_obs)) {
  sd_awayW[cc] <- (wh_obs[cc] - mean(wh_resampling[, cc])) / sd(wh_resampling[, cc])
}


# Plotting the tiles of variable importance ordering the variables in
# decreasing importance:

# For the bird species:
xx <- data.frame(value = sd_awayX, covariate = names(sd_awayX), y = 1)
xx <- xx[order(- xx$value), ]
xx$covariate <- factor(xx$covariate, levels = xx$covariate)

ggplot() + geom_tile(aes(x = covariate, y = y, fill = value), color = 'white', data = xx) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'none', axis.title = element_blank()) +
  scale_fill_gradient(low = '#BFF0B6', high = '#3B6E32') +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 8))

# For the plant species:
ww <- data.frame(value = sd_awayW, covariate = names(sd_awayW), y = 1)
ww <- ww[order(- ww$value), ]
ww$covariate <- factor(ww$covariate, levels = ww$covariate)

ggplot() + geom_tile(aes(x = covariate, y = y, fill = value), color = 'white', data = ww) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'none', axis.title = element_blank()) +
  scale_fill_gradient(low = '#D4DEF7', high = '#425075') +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 8))




# ------- PART B: Posterior probabilities based on the important covariates.

# Taking the posterior probabilities of interaction across the three chains,
# and setting the recorded interactions to NA:
use_Ls <- do.call(abind, c(lapply(all_res, function(x) x$all_pred[, , , 1]), along = 1))
use_mean_Ls <- apply(use_Ls, c(2, 3), mean)
use_mean_Ls[comb_A == 1] <- NA

# Which covariate is to be plotted. The ones we want are listed first.
wh_X <- 1
wh_W <- 1

# Showing only the species that have the covariate measured.
keep_birds <- which(!is.na(obs_X[, wh_X]))
keep_plants <- which(!is.na(obs_W[, wh_W]))
use_out <- use_mean_Ls[keep_birds, keep_plants]

# Because some species have identical values for the covariate, in order for
# plot to show all of them, we need to slightly pertube their values. That way,
# the increasing or decreasing order is not altered, but there is no overlap in
# the covariate values:
bird_cov <- obs_X[keep_birds, wh_X]
bird_cov <- bird_cov + rnorm(length(keep_birds), sd = sd(bird_cov) * 0.0001)
plant_cov <- obs_W[keep_plants, wh_W]
plant_cov <- plant_cov + rnorm(length(plant_cov), sd = sd(plant_cov) * 0.0001)

# Creating a data frame in which the species are ordered by their covariate
# values. This will allow us to plot the probability of interaction across the
# covariates in an interpretable way. We also note that we need to turn the
# covariates to factors in order for them to be plotted in the correct order.
plot_dta <- data.frame(cov_bird = rep(bird_cov, length(keep_plants)),
                       cov_plant = rep(plant_cov, each = length(keep_birds)),
                       probability = as.numeric(use_out))
plot_dta$use_cov_bird <- factor(as.numeric(factor(plot_dta$cov_bird)))
plot_dta$use_cov_plant <- factor(as.numeric(factor(plot_dta$cov_plant)))

g <- ggplot() +
  geom_raster(aes(x = use_cov_plant, y = use_cov_bird, fill = probability), data = plot_dta) +
  scale_fill_gradient(low = "#F5D4C7", high = "#02A65F", na.value = '#016B3B',
                      name = 'Posterior\ninteraction\nprobability\n', limits = c(0, 1)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 30),
        axis.title.y = element_text(vjust = - 1)) +
  ylab(expression(symbol('\256'))) + xlab(expression(symbol('\256')))

gridExtra::grid.arrange(g, left = textGrob("Bird information: Increasing Body Mass", rot = 90,
                                           x = 1.3, y = 0.57, gp = gpar(fontsize = 12)),
                        bottom = textGrob("Plant information: Increasing Fruit Diameter", 
                                          x = 0.435, y = 1.3, gp = gpar(fontsize = 12)),
                        vp=viewport(width=0.5, height=0.6))

# Why are we predicting way too many interactions?
F1 <- which(F_obs_m ==1, arr.ind = TRUE)
mean(obs_A[F1])
