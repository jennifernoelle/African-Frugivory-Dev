#--------------------- TO DO --------------------------#

## YOU REALLY HAVE TO SET THESE ##

# The directory where the analysis is performed:
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"

## SET THESE CORRESPONDING TO YOUR SAMPLER RUN ##

# Save results using convention: res_date_i.rda
date <- 'Final_old_1mod_MO'

date <- 'Final_old_p75mod_MO' # Saved with bad name oh well
dir <- 'Oldp75'

# Number of MCMC chains for our method and for the alternative method:
nchains <- 4

# Number of cross validation repetitions:
repetitions <- 30

# Number of cv samples per repetition
n.cv <- 100

## THESE SHOULD BE THE SAME IF YOU CLONED THE REPO ##

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
results_path <- paste0('Results/NewSampler/Final/', dir, "/")

save_files <- TRUE


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
library(ggridges)


setwd(wd_path)

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
load(paste0(data_path, 'traits_p_709_clean.dat')) # plant occurrences: most species have many studies, but some have 0? traits.plants
load(paste0(data_path, 'traits_m.dat'))
load(paste0(data_path, 'Subset_no_baboons.dat')) # list with appropriately subset mammals, plants, studies

# ------------------- SETUP --------------------------#

## Rename for convenience
Cu <- Cu_phylo
Cv <- Cv_phylo
obs_A <- A.obs.m
obs_W <- Obs_W[, 1:2]
obs_X <- Obs_X
traits.p <- traits.plants

## The following assignments were used in creating obs_OP
# Same study: 1, same site: 0.75
# Same country and habitat: 0.5, same region and habitat: 0.45, same habitat only: 0.25,
# Same country not habitat: 0.1, same region not habitat: 0.05


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
traits.p <- traits.p[wh_keep_p, ]
traits.m <- traits.m[wh_keep_m, ]

# The combined network:
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1 # 1 if any interactions observed

# Useful quantities
nM <- nrow(obs_A)
nP <- ncol(obs_A)
nS <- dim(obs_A)[3]

p.names <- colnames(obs_A[,,1]) 
m.names <- rownames(obs_A[,,1])
s.names <- dimnames(obs_A)[[3]]


# --------------- STEP 1: Getting the results together ----------------- #


### Put all chains together in a list
all_res <- NULL
for (ii in 1 : nchains) {
  cat("\n Chain = ", ii)
  res <- loadRData(paste0(results_path, 'res_',  date, '_', ii, '.dat')) 
  cat(" Dim is ", dim(res[[1]]))
  all_res[[ii]] <- res
  rm(res)
}

### Bind together predicted interactions, occurrence probabilities, occurrence indicators

# Number of posterior samples used:
use_Nsims <- dim(all_res[[1]]$all_pred)[1]

# Creating an array to bind results across chains:
pred_ours <- array(NA, dim = c(nchains * use_Nsims, nM, nP))
pred_OP <- array(NA, dim = c(nchains, nP, nS))

# Make sure this fills in dimensions correctly!
for (ii in 1 : nchains) {
  # Using the posterior samples of the L matrix:
  pred_ours[1 : use_Nsims + use_Nsims * (ii - 1), , ] <- all_res[[ii]]$all_pred[, , , 1] # this is the imputed L's from the model output (included bias corr)
  pred_OP[ii,, ] <- all_res[[ii]]$occ_plants # plant occurrence indicators at all locations
}
dimnames(pred_ours)[2 : 3] <- list(mammal = rownames(obs_A), plant = colnames(obs_A))
dimnames(pred_OP)[[2]] <- p.names

# Calculating the posterior means across mcmc iterations:
mean_pred <- mean_pred_na <- apply(pred_ours, c(2, 3), mean)
mean_pred_na[comb_A == 1] <- NA # Setting the recorded interactions to NA (so that they don't overpower the colors)
rownames(mean_pred) <- rownames(mean_pred_na)<- m.names
colnames(mean_pred) <- colnames(mean_pred_na) <- p.names

if(save_files){write.csv(mean_pred, paste0(results_path, "PosteriorNetworks/net_", date, ".csv"))}

# Calculating the posterior mean occurrence and occurrence probability across mcmc iterations:
mean_OP <- apply(pred_OP, c(2, 3), mean)



#-------------------- STEP 1B: SANITY CHECK -----------------------------#

mean(mean_pred)
summary(c(mean_pred))

# Check: model should be outputting 1 whenever there is a known interaction 
sum(comb_A)
sum(mean_pred == 1)
mean(mean_pred[comb_A == 1])

#write.csv(mean_pred, paste0(results_path, "post_network.csv"))
# Check: model should be outputting 1 whenever there is a known interaction 
sum(comb_A)
sum(mean_pred == 1)
mean(mean_pred[comb_A == 1])

# How many new interactions are predicted
sum(mean_pred ==1) # observed interactions
sum(mean_pred == 1) /length(mean_pred) # observed prevalence

sum(mean_pred > 0.5) # observed of likely interactions
sum(mean_pred > 0.5)/length(mean_pred) # posterior prevalence
sum(mean_pred_na > 0.5, na.rm = TRUE) # likely non-observed interactions
sum(mean_pred_na > 0.5, na.rm = TRUE)/sum(!(is.na(mean_pred_na)))

sum(mean_pred > 0.75) # observed of v likely interactions
sum(mean_pred > 0.75)/length(mean_pred) # posterior prevalence
sum(mean_pred_na > 0.75, na.rm = TRUE) # likely non-observed interactions
sum(mean_pred_na > 0.75, na.rm = TRUE)/sum(!(is.na(mean_pred_na)))

# ------------- 2. EXPLORING OCCURRENCE PROBABILITIES ------------#

# Goals: record mean of occurrence indicator for each category of obs op

mean(mean_OP) # overall mean of occurrence indicator for plants
mean(mean_OP[obs_OP!=1])
mean(mean_OP[obs_OP ==1]) # this should be 1

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

# Ridges plot
png(filename = paste0(results_path, "Figures/Occpostprior_ridges_", date, ".png"), width = 1000, height = 500)
ggplot(df2, aes(x =post_op, y= color, fill = color)) + 
  geom_density_ridges(alpha = 0.75) +
  ylab("") + 
  xlab("Posteior mean") + 
  guides(fill = "none")+ 
  facet_wrap(~type) + 
  theme_minimal() + 
  # scale_x_continuous(limits = c(0.65,0.8)) + 
  theme(text = element_text(family = "serif", size = 20))
dev.off()

# Debugging: this looks wrong
  # Why is the range of the posterior so much wider than with the new sampler?
  # And why is it so similar for all of the co-occurrence categories?
  # Why are the posterior means all in (0, 1/4, 1/2, 3/4, 1)
    # Go into the function, save all samples of obs_OP
    # Still need to figure this out, but also consider blocking indicator and prob in new sampler
  # Is it possible I'm saving the wrong info?

df2 %>% group_by(color) %>% 
  mutate(weight = 1/n()) %>% 
ggplot(aes(x =post_op)) + 
  geom_histogram(aes(weight = weight), stat = 'count') +
  ylab("") + 
  xlab("Posteior mean") + 
  guides(fill = "none")+ 
  facet_wrap(~color) + 
  theme_minimal() + 
  # scale_x_continuous(limits = c(0.65,0.8)) + 
  theme(text = element_text(family = "serif", size = 20))

table(df2$post_op, df2$color)

df2 %>% 
  group_by(color, post_op) %>%
  summarize(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% 
  ggplot(aes(x= post_op, y = freq)) + 
    geom_bar(stat = "identity") + 
    facet_wrap(~ color)

# Look at values across chains for a few plants x studies
#png(paste0(results_path, date, "_mixing.png"), width = 2000, height =1000)
par(mfrow = c(3,4))
for(p in 1:12){
  plant <- sample(1:nP,1)
  study <- sample(1:nS,1)
  plot(pred_OP[,plant,study], main = paste0(p.names[plant], " in Study ", s.names[study]))
}
#dev.off()
par(mfrow = c(1,1))

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
min_p_size <- 10

keep_p_groups <- names(which(p_size_cluster >= min_p_size))
keep_m_groups <- names(which(m_size_cluster >= min_m_size))

keep_m_index <- which(m_group %in% keep_m_groups)
keep_p_index <- which(p_group %in% keep_p_groups)

# Plotting those with minimum size as specified:
# We replaced observed interactions with black NA
#png(paste0(results_path, "heatmap_", date, ".png"), width = 1000, height =1000)
png(paste0(results_path, "paper1_heatmap_imputed", date, ".png"), width = 1000, height =1000)
superheat(X = plot_pred[keep_m_index, keep_p_index],
          membership.rows = m_group[keep_m_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 8,
          bottom.label.text.size = 8,
          bottom.label.size = 0.35, left.label.size = 0.22,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          heat.col.scheme = "grey", heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05), 
          title = "(b) Imputed frugivory network", 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 15, 
          title.alignment = "center")
dev.off()

png(paste0(results_path, "paper1_heatmap_observed", date, ".png"), width = 1000, height =1000)
superheat(X = comb_A[keep_m_index, keep_p_index],
          membership.rows = m_group[keep_m_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 8,
          bottom.label.text.size = 8,
          bottom.label.size = 0.35, left.label.size = 0.22,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          heat.col.scheme = "grey", heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05), 
          title = "(a) Observed frugivory network", 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 15, 
          title.alignment = "center")
dev.off()

# 
# 
# annotation.plants <- data.frame(p_group)
# rownames(annotation.plants) <- p.names
# annotation.mammals <- data.frame(m_group)
# rownames(annotation.mammals) <- m.names
# 
# pheatmap(mat = plot_pred[keep_m_index, keep_p_index],
#          cluster_rows = FALSE, cluster_cols = FALSE,
#          na_col = 'black',
#          color = colorRampPalette(c("white", "black"))(100),
#          annotation_row = annotation.mammals,
#          annotation_col = annotation.plants)



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

#plot(all_cor[,1], type = "l")

# Phylogenetic correlation more important for plants than animals


## Overview of latent factors for animals
if(nchains > 1){
  all_u <- abind::abind(all_res[[1]]$factors[[1]], all_res[[2]]$factors[[1]], along = 3)
  all_v <- abind::abind(all_res[[1]]$factors[[2]], all_res[[2]]$factors[[2]], along = 3)  
  if(nchains>2){for (cc in 3 : nchains) {
    all_u <- abind::abind(all_u, all_res[[cc]]$factors[[1]], along = 3)
    all_v <- abind::abind(all_v, all_res[[cc]]$factors[[2]], along = 3)
  }
  }
}else{
  all_u <- all_res[[1]]$factors[[1]]
  all_v <- all_res[[1]]$factors[[2]]
}

# Look at animal latent factors
mean_u <- apply(all_u, c(1,2), mean)
u_df <- data.frame(factor = paste0("Factor ", 1:10), MeanAbsValue = colMeans(abs(mean_u)), MeanValue = colMeans(mean_u))
u_df$factor <- factor(u_df$factor, levels = u_df$factor)
ggplot(u_df, aes(x = factor, y= MeanAbsValue)) + 
  geom_bar(stat = "identity") +
  theme_minimal() + 
  ggtitle("Animal Latent Factors") + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90), text = element_text(family = "serif", size = 20)) 


u_df_animals <- data.frame(mean_u)
colnames(u_df_animals) <- paste0("Factor_", 1:10)
u_df_animals$Mammal <- m.names
ggplot(u_df_animals, aes(x = Factor_1, y = Factor_2, label = Mammal)) + 
  geom_point() + 
  geom_label(data = subset(u_df_animals, Mammal %in% c("Papio_anubis", "Papio_cynocephalus"))) + 
  theme_minimal() + 
  ggtitle("Animal Latent Factors") + 
  theme(text = element_text(family = "serif", size = 20)) 


# Look at plant latent factors
mean_v <- apply(all_v, c(1,2), mean)
v_df <- data.frame(factor = paste0("Factor ", 1:10), MeanAbsValue = colMeans(abs(mean_v)), MeanValue = colMeans(mean_v))
v_df$factor <- factor(v_df$factor, levels = v_df$factor)
ggplot(v_df, aes(x = factor, y= MeanAbsValue)) + 
  geom_bar(stat = "identity") +
  theme_minimal() + 
  ggtitle("Plant Latent Factors") + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90), text = element_text(family = "serif", size = 20)) 


v_df_plants <- data.frame(mean_v)
colnames(v_df_plants) <- paste0("Factor_", 1:10)
v_df_plants$Species <- p.names 
v_df_plants <- left_join(v_df_plants, traits.plants)
ggplot(v_df_plants, aes(x = Factor_1, y = Factor_2, color = Order)) + 
  geom_point() + 
  #geom_label(data = subset(u_df_animals, Mammal %in% c("Papio_anubis", "Papio_cynocephalus"))) + 
  theme_minimal() + 
  ggtitle("Plant Latent Factors") + 
  theme(text = element_text(family = "serif", size = 20)) 



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

# Histogram of post probs from CV (only use this when full results are not run)

png(file = paste0(results_path, "Histogram", date, ".png"))
hist(c(our_preds), main = "CV Runs: African Frugivory fit with Biased Network, 
     Expert-Defined Occurrence #4", 
     xlab = "Posterior Interaction Probability")
dev.off()

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





