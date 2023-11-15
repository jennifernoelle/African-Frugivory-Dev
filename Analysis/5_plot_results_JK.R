#--------------------- TO DO --------------------------#

## Things to try
  # X Set occurrence to 1 always: good CV performance but calibration unclear
  # Set occurrence using sites but replace 0 with 0.5: with CV
  # X Review how mammal occurrences are currently computed: exactly the same as plants for site level
  # Fix rho
  # Add in more plant species and rerun
  # Subset birds and see if similar issues occur
  # Investigate pi_i
    # The fact that we get very high interaction prevalence predicted with only moderate
    # Observed interactionp prevalence suggests to me that we might be predicting very low
    # detection probabiliites - look at how these change using the different occurrence scenarios
  # Investigate latent factors
    # How many are important in the various scenarios


## YOU REALLY HAVE TO SET THESE ##

# The directory where the analysis is performed:
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"

## SET THESE CORRESPONDING TO YOUR SAMPLER RUN ##

# Save results using convention: res_date_i.rda
date <- 'May9_OccExpert1'

# Number of MCMC chains for our method and for the alternative method:
nchains <- 1

# Number of cross validation repetitions:
repetitions <- 4

# Number of cv samples per repetition
n.cv <- 35

## THESE SHOULD BE THE SAME IF YOU CLONED THE REPO ##

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
results_path <- 'Results/'
# Where the functions are available:
source_path <- 'HelperScriptsPapadogeorgou/'

# -------------------- LOAD PACKAGES AND FILES -----------------------#


# We want the github version of superheat
if(!("superheat" %in% rownames(installed.packages()))){
  install.packages("devtools")
  remotes::install_github("rlbarter/superheat")
}

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


setwd(wd_path)

source(paste0(source_path, 'draw_cm.R'))

# Define a useful function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv_phylo.dat'))
load(paste0(data_path, 'obs_A_mammals.dat'))
load(paste0(data_path, 'F_obs_default.dat'))
load(paste0(data_path, 'Obs_X.dat')) # mammal traits
load(paste0(data_path, 'Obs_W.dat')) # plant traits
load(paste0(data_path, 'obs_OM.dat')) # site level occurrence
load(paste0(data_path, 'obs_OP.dat')) # site level occurrence
load(paste0(data_path, 'traits_p_709_clean.dat')) # plant occurrences: most species have many studies, but some have 0?
load(paste0(data_path, 'traits_m.dat'))

# Data for subsetting
load(paste0(data_path, 'subset1_studies.dat')) 
load(paste0(data_path, 'subset1_mammals.dat')) 
load(paste0(data_path, 'subset1_plants.dat')) 

# ------------------- SETUP --------------------------#

## Rename for convenience
Cu <- Cu_phylo
Cv <- Cv_phylo
obs_A <- A.obs.m

wh_keep_m <- which(rownames(obs_A) %in% good.mammals)
wh_keep_p <- which(colnames(obs_A) %in% good.plants)
wh_keep_s <- which(unlist(dimnames(obs_A)[3]) %in% good.studies)

obs_A <- obs_A[wh_keep_m, wh_keep_p, wh_keep_s]
F_obs <- obs_F[wh_keep_m, wh_keep_p, wh_keep_s]
traits.p <- traits.plants[traits.plants$Species %in% good.plants, ]
traits.m <- traits.m[traits.m$Species %in% good.mammals, ]
obs_OP <- obs_OP[wh_keep_p, wh_keep_s]
obs_OM <- obs_OM[wh_keep_m, wh_keep_s]
obs_X <- Obs_X[wh_keep_m, ]
obs_W <- Obs_W[wh_keep_p, ]


# The combined network:
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1 # 1 if any interactions observed

# Useful quantities
nM <- nrow(obs_A)
nP <- ncol(obs_A)
nS <- dim(obs_A)[3]

plant.names <- colnames(obs_A[,,1]) #traits.p.709$Species
m.names <- rownames(obs_A[,,1])#traits.m$Species



# --------------- STEP 1: Getting the results together ----------------- #


### Put all chains together in a list
all_res <- NULL
for (ii in 1 : nchains) {
  cat("\n Chain = ", ii)
  res <- loadRData(paste0(results_path, 'res_',  date, '_', ii, '.dat')) 
  cat(" Dim is ", dim(res[[1]]))
  all_res[[ii]] <- res
}

### Bind together predicted interactions

# Number of posterior samples used:
use_Nsims <- dim(all_res[[1]]$all_pred)[1]

# Creating an array to bind results across chains:
pred_ours <- array(NA, dim = c(nchains * use_Nsims, nM, nP))
p_is <- array(NA, dim = c(nchains * use_Nsims, nM))
p_js <- array(NA, dim = c(nchains * use_Nsims, nP))
for (ii in 1 : nchains) {
  # Using the posterior samples of the L matrix:
  pred_ours[1 : use_Nsims + use_Nsims * (ii - 1), , ] <- all_res[[ii]]$all_pred[, , , 1] # this is the imputed L's from the model output (included bias corr)
  p_is[1 : use_Nsims + use_Nsims * (ii - 1), ] <- all_res[[ii]]$p_detect[[1]] # mammal detection probabilities
  p_js[1 : use_Nsims + use_Nsims * (ii - 1), ] <- all_res[[ii]]$p_detect[[2]] # plant detection probabilities
}
dimnames(pred_ours)[2 : 3] <- list(mammal = rownames(obs_A), plant = colnames(obs_A))
colnames(p_is) <- rownames(obs_A)
colnames(p_js) <- colnames(obs_A)

# Calculating the posterior means across mcmc iterations:
mean_pred <- mean_pred_na <- apply(pred_ours, c(2, 3), mean)
mean_pi <- colMeans(p_is)
mean_pj <- colMeans(p_js)

#png(file = paste0(results_path, "Pi_i_", date, ".png"))
hist(mean_pi, main = "Animals Detection Probabilities \n with Expert Occ 3", xlab = "Posterior Mean")
abline(v = mean(mean_pi))
#dev.off()

#png(file = paste0(results_path, "Pi_j_", date, ".png"))
hist(mean_pj, main = "Plants Detection Probabilities \n with Expert Occ 3", xlab = "Posterior Mean")
abline(v = mean(mean_pj))
#dev.off()

# Setting the recorded interactions to NA (so that they don't overpower the colors)
mean_pred_na[comb_A == 1] <- NA
# 
# ### Bind together interaction probabilities: the bias corrected pred above are preferred
# prob_ours <- array(NA, dim = c(nchains * use_Nsims, nM, nP))
# for (ii in 1 : nchains) {
#   # Using the posterior samples of the L matrix:
#   prob_ours[1 : use_Nsims + use_Nsims * (ii - 1), , ] <- all_res[[ii]]$all_pred[, , , 2] # Original 2: mod_pL1s, 3 is pL1s
# }
# dimnames(prob_ours)[2 : 3] <- list(mammal = rownames(obs_A), plant = colnames(obs_A))
# 
# # Look at mixing for a few select probabilities 
# #png(paste0(results_path, date, "_mixing.png"), width = 2000, height =1000)
# par(mfrow = c(3,4))
# for(p in 1:12){
#   mammal <- sample(1:nM, 1)
#   plant <- sample(1:nP,1)
#   plot(prob_ours[,mammal,plant], type = "l", main = paste0(rownames(obs_A)[mammal], " x ", colnames(obs_A)[plant]))
# }
# #dev.off()
# par(mfrow = c(1,1))


#-------------------- STEP 1B: SANITY CHECK -----------------------------#

mean(mean_pred)
summary(c(mean_pred))

png(file = paste0(results_path, "Histogram", date, ".png"))
hist(c(mean_pred), main = "African Frugivory fit with Biased Network \n H = 10, with shrinkage and feedback, 
     Expert-Defined Occurrence #1", 
     xlab = "Posterior Interaction Probability")
abline(v = mean(mean_pred))
dev.off()

# ----------- STEP 2: PLOTTING THE HEATMAP --------------------#

# Creating the clusters that will be used
# The following two lines specify that horizontal and vertical lines in our
# plot will separate species by taxonomic families:
m_group <- traits.m$Family
p_group <- traits.p$Family

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
min_p_size <- 1

keep_p_groups <- names(which(p_size_cluster >= min_p_size))
keep_m_groups <- names(which(m_size_cluster >= min_m_size))

keep_m_index <- which(m_group %in% keep_m_groups)
keep_p_index <- which(p_group %in% keep_p_groups)


# Plotting those with minimum size as specified:
# We replaced observed interactions with black NA
png(paste0(results_path, "heatmap_", date, ".png"), width = 1000, height =1000)
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
          heat.pal.values = seq(0, 1, by = 0.05), 
          title = "Expert-defined occurrence #1", 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 20)
dev.off()
# 
# annotation.plants <- data.frame(p_group)
# rownames(annotation.plants) <- plant.names
# annotation.mammals <- data.frame(m_group)
# rownames(annotation.mammals) <- m.names
# 
# pheatmap(mat = plot_pred[keep_m_index, keep_p_index],
#          cluster_rows = FALSE, cluster_cols = FALSE,
#          na_col = 'black', 
#          color = colorRampPalette(c("white", "black"))(100), 
#          annotation_row = annotation.mammals, 
#          annotation_col = annotation.plants)
# 


# Note: these results are strange: interactions occur to occur with much more frequency than observed
# Lots of dark grey in rows with very little black


# --------------- STEP 3: Taxonomic correlation of latent factors ----------------- #
if(nchains > 1){
  all_cor <- abind::abind(all_res[[1]]$correlations, all_res[[2]]$correlations, along = 3)
  if(nchains>2){for (cc in 3 : nchains) {
    all_cor <- abind::abind(all_cor, all_res[[cc]]$correlations, along = 3)
    }
  }
}else{
  all_cor <- all_res[[1]]$correlations
}

# Posterior means and 95% credible intervals for the rho parameters in the
# latent factors for bird and plant species:
apply(all_cor, 2, mean)
apply(all_cor, 2, quantile, probs = c(0.025, 0.975))

plot(all_cor[,1], type = "l")

# Phylogenetic correlation more important for plants than animals

# --------------- STEP 4: Cross validation results ----------------- #

# Getting the results together (held out indicies and predictions)
all_indices <- array(NA, dim = c(repetitions, n.cv, 2))
our_preds <- array(NA, dim = c(repetitions, nM, nP))

# for (rr in 1 : repetitions) {
#   cv_indices <- loadRData(paste0(results_path, date, '_cv_indices_', rr, '.rda'))
#   pred <- loadRData(paste0(results_path, date, '_pred_', rr, '.rda'))
#   all_indices[rr, , ] <- cv_indices
#   our_preds[rr, , ] <- pred
# }

for (rr in 1 : repetitions) {
  cv_indices <- loadRData(paste0(results_path, 'cv_indices_', date, '_', rr, '.dat'))
  pred <- loadRData(paste0(results_path, 'pred_', date, '_', rr, '.dat'))
  all_indices[rr, , ] <- cv_indices
  our_preds[rr, , ] <- pred
}


# Predictions of the held out data from model: 100 indices held out each time - always interactions
pred <- array(NA, dim = c(repetitions, n.cv))
for (rr in 1 : repetitions) {
  for (ii in 1 : n.cv) {
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
                       stat = rep(c('Pred:Overall (mean)', 'Pred:Overall (median)'), repetitions))

# Plotting cross validation results:
png(filename = paste0(results_path, "cv_res_", date, ".png"))
ggplot(data = plot_dta) +
  geom_boxplot(aes(x = stat, y = value)) +
  theme_bw() +
  ylab('') +
  xlab('') +
  ggtitle('Out of sample performance', subtitle = 'Expert-defined occurrence #1') +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) + 
  theme(text = element_text(size = 20))
dev.off()




# --------------- STEP 5: Variable importance measure ----------------- #


# ------- PART A: Plotting the variables in order of importance:

load(paste0(results_path, 'rsq_obs_X.dat'))
load(paste0(results_path, 'rsq_obs_W.dat'))
load(paste0(results_path, 'rsq_resampling_X.dat'))
load(paste0(results_path, 'rsq_resampling_W.dat'))


# Covariate names that are nicer for plotting:
good_namesX <- c('Generation Length', 'Log(Body Mass)', 'Log(Brain Mass/Body Mass)', 'IUCN Endangered', 'Forest Habitat') 
good_namesW <-  c('Fruit Length (mm)', 'Mean Wood Density') 

# Here are Georgia's covariates in case we want to try to acquire additional information
## #c('Body Mass', 'Gape Size', 'Large*', 'Fruit\nDependent*', 'Endangered*')
# # c('Fruit\nDiameter', 'Fruit\nLength', 'Seed\nDiameter', 'Seed\nLength', 'Native*',
#                  'Tree*', 'Black\nFruit*', 'Red\nFruit*', 'Yellow/Orange\nFruit*', 'Green\nFruit*',
#                  'Lipid*', 'Endangered*')

# Calculating the number of permuted standard deviations away from the mean.

# Starting from the mammal covariates:
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

# For the mammal species:
xx <- data.frame(value = sd_awayX, covariate = names(sd_awayX), y = 1)
xx <- xx[order(- xx$value), ]
xx$covariate <- factor(xx$covariate, levels = xx$covariate)

png(filename = paste0(results_path, "VarImp_m_", date, ".png"))
ggplot() + geom_tile(aes(x = covariate, y = y, fill = value), color = 'white', data = xx) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'none', axis.title = element_blank()) +
  scale_fill_gradient(low = '#BFF0B6', high = '#3B6E32') +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 8))
dev.off()

# For the plant species:
ww <- data.frame(value = sd_awayW, covariate = names(sd_awayW), y = 1)
ww <- ww[order(- ww$value), ]
ww$covariate <- factor(ww$covariate, levels = ww$covariate)

png(filename = paste0(results_path, "VarImp_p_", date, ".png"))
ggplot() + geom_tile(aes(x = covariate, y = y, fill = value), color = 'white', data = ww) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'none', axis.title = element_blank()) +
  scale_fill_gradient(low = '#D4DEF7', high = '#425075') +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 8))
dev.off()


# ------- PART B: Posterior probabilities based on the important covariates.

# Taking the posterior probabilities of interaction across the three chains,
# and setting the recorded interactions to NA:
use_Ls <- do.call(abind, c(lapply(all_res, function(x) x$all_pred[, , , 1]), along = 1))
use_mean_Ls <- apply(use_Ls, c(2, 3), mean)
use_mean_Ls[comb_A == 1] <- NA

# Which covariate is to be plotted. The ones we want are listed first.
wh_X <- 3 # brain:body mass
wh_W <- 1

# Showing only the species that have the covariate measured.
keep_m <- which(!is.na(obs_X[, wh_X]))
keep_p <- which(!is.na(obs_W[, wh_W]))
use_out <- use_mean_Ls[keep_m, keep_p]

# Because some species have identical values for the covariate, in order for
# plot to show all of them, we need to slightly perturb their values. That way,
# the increasing or decreasing order is not altered, but there is no overlap in
# the covariate values:
m_cov <- obs_X[keep_m, wh_X]
m_cov <- m_cov + rnorm(length(keep_m), sd = sd(m_cov) * 0.0001)
p_cov <- obs_W[keep_p, wh_W]
p_cov <- p_cov + rnorm(length(p_cov), sd = sd(m_cov) * 0.0001)

# Creating a data frame in which the species are ordered by their covariate
# values. This will allow us to plot the probability of interaction across the
# covariates in an interpretable way. We also note that we need to turn the
# covariates to factors in order for them to be plotted in the correct order.
plot_dta <- data.frame(cov_m = rep(m_cov, length(keep_p)),
                       cov_p = rep(p_cov, each = length(keep_m)),
                       probability = as.numeric(use_out))
plot_dta$use_cov_m <- factor(as.numeric(factor(plot_dta$cov_m)))
plot_dta$use_cov_p <- factor(as.numeric(factor(plot_dta$cov_p)))

g <- ggplot() +
  geom_raster(aes(x = use_cov_p, y = use_cov_m, fill = probability), data = plot_dta) +
  scale_fill_gradient(low = "#F5D4C7", high = "#02A65F", na.value = '#016B3B',
                      name = 'Posterior\ninteraction\nprobability\n', limits = c(0, 1)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 30),
        axis.title.y = element_text(vjust = - 1)) +
  ylab(expression(symbol('\256'))) + xlab(expression(symbol('\256')))

gridExtra::grid.arrange(g, left = textGrob("Mammal information: Increasing Log Body Mass", rot = 90,
                                           x = 1.3, y = 0.57, gp = gpar(fontsize = 12)),
                        bottom = textGrob("Plant information: Increasing Fruit Length", 
                                          x = 0.435, y = 1.3, gp = gpar(fontsize = 12)),
                        vp=viewport(width=0.5, height=0.6))





