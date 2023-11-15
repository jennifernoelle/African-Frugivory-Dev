#--------------------- TO DO --------------------------#


## YOU REALLY HAVE TO SET THESE ##

# The directory where the analysis is performed:
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"

## SET THESE CORRESPONDING TO YOUR SAMPLER RUN ##

# Save results using convention: res_date_i.rda
# All of the results below had failures on scheduled cores. Try re-running without scavenger?
date <- 'Final_new_p0' # scheduled cores 3, 4, 5, 6, 8, 11, 12, 20, 22, 23 did not deliver results
date <- 'Final_new_p50' #Missing:1,4,10,25,26 ; scheduled cores 1, 2, 3, 4, 10, 11, 14, 21, 24, 25, 26 did not deliver results; 1, 4, 5, 10, 12, 13, 18, 20, 23, 25, 26

date <- 'Final_new_p75' #Missing: 1,6,8,9,12 scheduled cores 1, 4, 6, 8, 9, 12, 13, 15, 17, 28 did not deliver results; 1, 2, 5, 6, 7, 8, 9, 12, 14, 20, 24

date <- 'Final_old_p0' #Missing: 2,6,16, scheduled cores 1, 2, 4, 6, 12, 16, 19, 20, 25, 26, 29 did not deliver results; ; 2, 5, 6, 7, 8, 9, 10, 11, 16, 24, 28 did not deliver results
date <- 'Final_old_p50' #Missing: 6,11,14,21, 29 scheduled cores 1, 6, 9, 10, 11, 13, 14, 15, 20, 21, 29 did not deliver results; 6, 7, 8, 11, 14, 17, 18, 21, 24, 29, 30 did not deliver results


# Number of cross validation repetitions:
repetitions <- 30

# Number of cv samples per repetition
n.cv <- 100

## THESE SHOULD BE THE SAME IF YOU CLONED THE REPO ##

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
results_path <- 'Results/NewSampler/'
# Where the functions are available:
source_path <- 'HelperScriptsJK/'

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
library(latex2exp)


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
load(paste0(data_path, 'Subset_no_baboons.dat')) # list with appropriately subset mammals, plants, studies

# ------------------- SETUP --------------------------#

## Rename for convenience
obs_A <- A.obs.m

## Subset to for common plants and to remove baboons
wh_keep_m <- which(rownames(obs_A) %in% nobab.list[[1]])
#common.plants <- names(which(apply(obs_A, 2, sum)>10))
wh_keep_p <- which(colnames(obs_A) %in% nobab.list[[2]])
wh_keep_s <- which(unlist(dimnames(obs_A)[3]) %in% nobab.list[[3]])

obs_A <- obs_A[wh_keep_m, wh_keep_p, wh_keep_s]
obs_F <- obs_F[wh_keep_m, wh_keep_p, wh_keep_s]
obs_OP <- obs_OP[wh_keep_p, wh_keep_s]
obs_OM <- obs_OM[wh_keep_m, wh_keep_s]

## The following assignments were used in creating obs_OP
# Same study: 1, same site: 0.75
# Same country and habitat: 0.5, same region and habitat: 0.45, same habitat only: 0.25,
# Same country not habitat: 0.1, same region not habitat: 0.05

## Replace with weak guess
obs_OP <- ifelse(obs_OP == 1, 1, 0.5)

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

# Useful values
nP <- ncol(obs_A)
nM <- nrow(obs_A)
nS <- dim(obs_A)[3]

m.names <- rownames(obs_A)
p.names <- colnames(obs_A)
s.names <- unlist(dimnames(obs_A)[3])


# --------------- Cross validation results ----------------- #

# Getting the results together (held out indicies and predictions)
all_indices <- array(NA, dim = c(repetitions, n.cv, 2))
our_preds <- array(NA, dim = c(repetitions, nM, nP))

# Getting the predicitons together: this works even if some threads failed
rr.good <- 1
for (rr in 1 : repetitions) {
  tryCatch({
    print(rr)
    print(rr.good)
    cv_indices <- loadRData(paste0(results_path, 'cv_indices_', date, '_', rr, '.dat'))
    pred <- loadRData(paste0(results_path, 'pred_', date, '_', rr, '.dat'))
    all_indices[rr.good, , ] <- cv_indices
    our_preds[rr.good, , ] <- pred # all predictions
    rr.good <- rr.good + 1
  }, error = function(e){cat("Error : ", conditionMessage(e), "\n")})
}

# Remove errors
n.fail <- rr.good - 1 #sum(is.na(all_indices[,1,1]))
all_indices <- all_indices[1:(repetitions-n.fail),,]
our_preds <- our_preds[1:(repetitions-n.fail),,]

# Predictions of the held out data from model: 100 indices held out each time - always interactions
pred <- array(NA, dim = c(repetitions-n.fail, n.cv))
for (rr in 1 : (repetitions-n.fail)) {
  for (ii in 1 : n.cv) {
    pred[rr, ii] <- our_preds[rr, all_indices[rr, ii, 1], all_indices[rr, ii, 2]] # just the cv preds
  }
}



# Check: model should be outputting 1 whenever there was a non-heldout observed interaction
sum(comb_A)
apply(our_preds, 1, function(x) sum(x==1))


# Average and median probability of interaction based on the overall data
overall_mean <- cbind(apply(our_preds, 1, mean))
overall_median <- cbind(apply(our_preds, 1, median))

# Average and median in the held out data.
pred_mean <- apply(pred, 1, mean)
pred_median <- apply(pred, 1, median)

# Pseudo precision: 
mean(pred_mean/overall_mean)
mean(pred_median/overall_median)

# Pseudo accuracy: 
#mean(pred_mean)
sum(pred>0.5)/length(pred) # what proportion of true interactions are predicted as "likely" >0.5; 0.88 for new sampler - WHY IS THIS ALWAYS 0.888???
sum(pred>0.75)/length(pred) # what proportion of true interactions are predicted as "v likely" >0.75; 0.64 for new sampler


# Investigate: why are we always getting the same value of for pseudo accuracy 1:
# They aren't the same interaction pairs
# This is weird
# pred.0.new <- which(pred.new < 0.5)
# pred.0.old <- which(pred.old < 0.5)
# 
# pred.new[pred.0.new]
# pred.old[pred.0.old]
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
  ggtitle('Out of sample performance', subtitle = 'New Sampler') +
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) + 
  scale_y_continuous(limits = function(x) c(1.1, 1.6), n.breaks = 5) + 
  theme(text = element_text(size = 20))
dev.off()

# Slides plot
plot_dta <- data.frame(value = rbind(pred_mean / overall_mean, pred_median / overall_median), 
                       stat = rep(c('mean_ratio', 'median_ratio'), repetitions))

png(filename = paste0(results_path, "Figures/slides_cv_Occ0cnobab", date, ".png"), width = 750, height = 750)
ggplot(data = plot_dta) +
  geom_boxplot(aes(x = stat, y = value)) +
  theme_bw() +
  ylab('') +
  xlab('') +
  ggtitle('Model fit with Expert Occurrence') +
  theme(legend.position = 'none') +
  scale_x_discrete(labels = c('mean_ratio' = parse(text = TeX('$ \\frac{Mean(\\hat{p}_{holdout})}{Mean(\\hat{p}_{overall})}$')), 
                              'median_ratio' = parse(text = TeX('$ \\frac{Median(\\hat{p}_{holdout})}{Median(\\hat{p}_{overall}})$')))) +
  scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) + 
  theme(text = element_text(size = 35)) 
dev.off()
  

# Histogram of post probs from CV (only use this when full results are not run)
png(file = paste0(results_path, "Figures/Histogram", date, ".png"))
hist(c(our_preds), main = "CV Runs: African Frugivory fit with Biased Network, 
     Expert Occurrence 0b NB Holdouts", 
     xlab = "Posterior Interaction Probability")
dev.off()


# ----------- LOOK AT SAMPLED OCCURRENCE PROBS ----------------# 

### Put all chains together in a list
all_p <- NULL
for (rr in 1 : repetitions) {
  cat("\n Chain = ", rr)
  res <- loadRData(paste0(results_path, 'p_occ_',  date, '_', rr, '.dat')) 
  cat(" Dim is ", dim(res[[1]]))
  all_p[[rr]] <- res
}

# Number of posterior samples used:
use_Nsims <- dim(all_p[[1]][[1]])[1]

# Creating an array to bind results across chains:
p_occ <- array(NA, dim = c(repetitions*use_Nsims, nP, nS))
p_accept <- array(NA, dim = c(repetitions, nP, nS))

for (rr in 1 : repetitions) {
  # Using the posterior samples of the L matrix:
  p_occ[1 : use_Nsims + use_Nsims * (rr - 1), , ] <- all_p[[rr]][[1]]
  p_accept[rr, , ] <- all_p[[rr]][[2]]
}

# Get posterior mean occurrence probabilities
mean_p_occ <- apply(p_occ, c(2,3), mean)
mean_accept <- apply(p_accept, c(2,3), mean)

# Look at mixing for a few select probabilities
png(paste0(results_path, date, "_mixing.png"), width = 2000, height =1000)
par(mfrow = c(3,4))
plot.no <- 1
while(plot.no < 13){
  plant <- sample(1:nP,1)
  study <- sample(1:nS, 1)
  rep <- sample(1:repetitions, 1)
  this.p <-(all_p[[rep]])[[1]][, plant, study]
  
  if (mean(this.p)==1){next} # Skip occurrences probs which aren't sampled
  
  plot(this.p, type = "l", main = paste0(colnames(obs_A)[plant], " x ", s.names[study]))
  abline(h = obs_OP[plant, study], col = "blue", lwd = 3, lty = 2)
  abline(h = mean(p_occ[,plant, study]), col = "red", lwd = 3, lty = 3)
  plot.no <- plot.no + 1
}
dev.off()
par(mfrow = c(1,1))

plot(c(mean_p_occ), c(obs_OP))
data.frame(prior.mean = as.factor(c(obs_OP)), post.mean = c(mean_p_occ)) %>% 
ggplot(aes(x = prior.mean, y = post.mean)) +
   geom_boxplot()

summary(c(mean_p_occ))
hist(c(mean_p_occ))

png(paste0(results_path, date, "_pocc_hist.png"), width = 2000, height =1000)
data.frame(p_occ = c(mean_p_occ)) %>% 
  filter(p_occ != 1) %>% 
  ggplot() + 
  geom_histogram(aes(x = p_occ), color = 'black', fill = 'darkgrey') + 
  geom_vline(aes(xintercept = 0.75)) + 
  geom_text(aes(x = 0.75, label = "\nPrior Mean", y = 90, angle = 90)) + 
  xlab("Posterior Mean Occurrence Probability") + 
  ylab("") +
  ggtitle("Plant Occurrence")
  theme_minimal()
dev.off()




# ----------- STEP 2: PLOTTING THE HEATMAP --------------------#

# Create plotting data
plot_pred <- apply(our_preds, c(2,3), mean)
write.csv(plot_pred, paste0(results_path, "post_network_", date, ".csv"))
plot_pred[comb_A == 1] <- NA

# Creating the clusters that will be used
# The following two lines specify that horizontal and vertical lines in our
# plot will separate species by taxonomic families:
m_group <- traits.m$Family
p_group <- traits.p$Family

# Calculating the size of each cluster, will be used when plotting results
# for families of certain size:
m_size_cluster <- sapply(unique(m_group), function(x) sum(m_group == x))
p_size_cluster <- sapply(unique(p_group), function(x) sum(p_group == x))

# Set the minimum cluster size that should be plotted. For the results of the
# manuscript, we set min_bird_size to 10, and min_plant_size to 20. Setting both
# to 0 will produce the full results.
min_m_size <- 1
min_p_size <- 5

keep_p_groups <- names(which(p_size_cluster >= min_p_size))
keep_m_groups <- names(which(m_size_cluster >= min_m_size))

keep_m_index <- which(m_group %in% keep_m_groups)
keep_p_index <- which(p_group %in% keep_p_groups)


# Plotting those with minimum size as specified:
# We replaced observed interactions with black NA
png(paste0(results_path, "Figures/heatmap_cv_", date, ".png"), width = 1000, height =1000)
png(paste0(results_path, "slides_heatmap_cv_exp0cnb2.png"), width = 750, height =750)
superheat(X = plot_pred[keep_m_index, keep_p_index],
          membership.rows = m_group[keep_m_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 5,
          bottom.label.text.size = 5,
          bottom.label.size = 0.4, left.label.size = 0.2,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          heat.col.scheme = "grey", heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05), 
          title = "Model fit with expert occurrence", 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 13)
dev.off()


## Results summary: plot_pred has NA for observed interaction
sum(comb_A)
sum(comb_A)/length(comb_A)
(sum(plot_pred > 0.75, na.rm = TRUE) + sum(comb_A))/length(plot_pred)
sum(plot_pred > 0.75, na.rm = TRUE)
(sum(plot_pred > 0.5, na.rm = TRUE) + sum(comb_A))/length(plot_pred)
sum(plot_pred > 0.5, na.rm = TRUE)


