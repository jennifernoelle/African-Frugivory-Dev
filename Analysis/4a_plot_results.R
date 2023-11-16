#--------------------- TO DO --------------------------#

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
results_path <- 'Results/'

## SET THESE CORRESPONDING TO YOUR SAMPLER RUN ##

# Save results using convention: res_date_i.rda
date <- 'Expert_prior'

# Number of MCMC chains for our method and for the alternative method:
nchains <- 4

# Number of cross validation repetitions:
repetitions <- 30

# Number of cv samples per repetition
n.cv <- 100

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
#library(ridgeline)
library(ggridges)


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
pred_ours <- array(NA, dim = c(nchains * use_Nsims, nM, nP)) # Saved all samples of this
pred_OP <- array(NA, dim = c(nchains, nP, nS)) # Just posterior means for each chain
pred_pOP <- array(NA, dim = c(nchains, nP, nS)) # Just posterior means for each chain

# Make sure this fills in dimensions correctly!
for (ii in 1 : nchains) {
  # Using the posterior samples of the L matrix:
  pred_ours[1 : use_Nsims + use_Nsims * (ii - 1), , ] <- all_res[[ii]]$all_pred[, , , 1] # this is the imputed L's from the model output (included bias corr)
  pred_OP[ii,, ] <- all_res[[ii]]$occ_plants # plant occurrence indicators at all locations
  pred_pOP[ii,, ] <- all_res[[ii]]$p_occ_plants$p_OPs # plant occurrence probabilities at all locations
}
dimnames(pred_ours)[2 : 3] <- list(mammal = rownames(obs_A), plant = colnames(obs_A))
dimnames(pred_OP)[[2]] <- dimnames(pred_pOP)[[2]]<-  p.names


# Calculating the posterior mean edge probability across mcmc iterations:
mean_pred <- mean_pred_na <- apply(pred_ours, c(2, 3), mean)
mean_pred_na[comb_A == 1] <- NA # Setting the recorded interactions to NA (so that they don't overpower the colors)

# Calculating the posterior mean occurrence and occurrence probability across mcmc iterations:
mean_OP <- apply(pred_OP, c(2, 3), mean)
mean_pOP <- apply(pred_pOP, c(2, 3), mean)

rownames(mean_pred) <- m.names
colnames(mean_pred) <- p.names

if(save_files){write.csv(mean_pred, paste0(results_path, "PosteriorNetwork_", date, ".csv"))}



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
mean(mean_pOP[obs_OP ==1]) # this should be 1

par(mfrow = c(1,2))
#png(file = paste0(results_path, "Figures/OccPlants_", date, ".png"))
hist(mean_OP, main = "Plants occurrence indictors", xlab = "Posterior Mean")
abline(v = mean(mean_OP))
hist(obs_OP, main = "Plants occurrence indictors", xlab = "Prior Mean")
abline(v = mean(obs_OP))
#dev.off()

par(mfrow = c(1,2))
#png(file = paste0(results_path, "Figures/pOccPlants_", date, ".png"))
hist(obs_OP, main = "Plants occurrence probabilities", xlab = "Prior Mean")
abline(v = mean(obs_OP))
hist(mean_pOP, main = "Plants occurrence probabilities", xlab = "Posterior Mean")
abline(v = mean(mean_pOP))
#dev.off()

# Code needs to be changed depending on occurrence scenario
data.frame(obs_op = c(ifelse(obs_OP ==1, 1, 0.75)), post_op = c(mean_OP)) %>%
  ggplot(aes(x = as.factor(obs_op), y = post_op)) + 
  geom_boxplot()

df1 <- data.frame(obs_op = c(ifelse(obs_OP ==1, 1, 0.75)), post_op = c(mean_pOP), 
                  color = factor(obs_OP, labels = c("Same region only", "Same country only",
                                                    "Same habitat only", "Same habitat and region", 
                                                    "Same habitat and country", "Same site", "Same study")), 
                  type = "Occurrence baseline probabilities") %>% 
  filter(obs_op == 0.75)
df2 <- data.frame(obs_op = c(ifelse(obs_OP ==1, 1, 0.75)), post_op = c(mean_OP),
                  color = factor(obs_OP, labels = c("Same region only", "Same country only",
                                                    "Same habitat only", "Same habitat and region", 
                                                    "Same habitat and country", "Same site", "Same study")), 
                  type = "Occurrence indicators") %>% 
  filter(obs_op == 0.75)
# Ridges plot
png(filename = paste0(results_path, "Figures/Occpostprior_ridges_", date, ".png"), width = 1000, height = 500)
rbind(df1, df2)%>%
  transform(type = factor(type, levels = c("Occurrence baseline probabilities", "Occurrence indicators"))) %>% 
  ggplot(aes(x =post_op, y= color, fill = color)) + 
  geom_density_ridges(alpha = 0.75) +
  ylab("") + 
  xlab("Posteior mean") + 
  guides(fill = "none")+ 
  facet_wrap(~type) + 
  theme_minimal() + 
  scale_x_continuous(limits = c(0.65,0.8)) + 
  theme(text = element_text(family = "serif", size = 20))
dev.off()

# Scatterplot
png(filename = paste0(results_path, "Figures/Occpostprior_", date, ".png"), width = 1000, height = 500)
rbind(df1, df2)%>%
  ggplot(aes(x =obs_op, y = post_op, color = color)) + 
  geom_point(alpha = 0.25) + 
  geom_abline(slope = 1, intercept = 0) + 
  xlab("Prior mean") + 
  ylab("Posterior mean") + 
  #ggtitle("Occurence probabilities")  +
  guides(color = guide_legend(title = "Plant previously observed in", 
                              override.aes = list(alpha = 1)))+
  facet_wrap(~type) + 
  theme_minimal()+ 
  theme(text = element_text(family = "serif", size = 20))
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
min_p_size <- 10

keep_p_groups <- names(which(p_size_cluster >= min_p_size))
keep_m_groups <- names(which(m_size_cluster >= min_m_size))

keep_m_index <- which(m_group %in% keep_m_groups)
keep_p_index <- which(p_group %in% keep_p_groups)

# Plotting those with minimum size as specified:
# We replaced observed interactions with black NA
#png(paste0(results_path, "Figures/heatmap_imputed", date, ".png"), width = 1000, height = 800)
superheat(X = plot_pred[keep_m_index, keep_p_index],
          membership.rows = m_group[keep_m_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 9,
          bottom.label.text.size = 10,
          bottom.label.size = 0.45, left.label.size = 0.25,
          #bottom.label = "none",
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          #heat.col.scheme = "grey", heat.na.col = 'black',
          heat.pal = c("white", "black"), heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05), 
          #title = "(b) Imputed frugivory network", 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 15
          #,title.alignment = "center"
          ,legend.text.size = 25
          #legend = FALSE
)
#dev.off()


# Plot imputed network with no legend
#png(paste0(results_path, "Figures/heatmap_imputed_noleg_", date,  ".png"), width = 1000, height = 510)
superheat(X = plot_pred[keep_m_index, keep_p_index],
          membership.rows = m_group[keep_m_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 9,
          bottom.label = "none",
          bottom.label.text.size = 10,
          bottom.label.size = 0.35, left.label.size = 0.25,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          #heat.col.scheme = "grey", 
          #heat.na.col = 'black',
          heat.pal = c("white", "black"), , heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05), 
          #title = "(a) Observed frugivory network", 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 15, 
       #   title.alignment = "center",
          legend = FALSE
)
#dev.off()


# Plot observed network
comb_A_na <- ifelse(comb_A ==1, NA, 0)
#png(paste0(results_path, "Figures/heatmap_observed.png"), width = 1000, height = 510)
superheat(X = comb_A[keep_m_index, keep_p_index],
          membership.rows = m_group[keep_m_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 9,
          bottom.label = "none",
          bottom.label.text.size = 10,
          bottom.label.size = 0.35, left.label.size = 0.25,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          #heat.col.scheme = "grey", 
          #heat.na.col = 'black',
          heat.pal = c("white", "black"),
          heat.pal.values = seq(0, 1, by = 0.05), 
          #title = "(a) Observed frugivory network", 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 15, 
          #title.alignment = "center",
          legend = FALSE
)
#dev.off()

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

# --------------- STEP 4: Cross validation results ----------------- #

# Getting the results together (held out indicies and predictions)
all_indices <- array(NA, dim = c(repetitions, n.cv, 2))
our_preds <- array(NA, dim = c(repetitions, nM, nP))

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

# Pseudo precision: 
mean(pred_mean/overall_mean)
mean(pred_median/overall_median)

# Pseudo accuracy: 
#mean(pred_mean)
sum(pred>0.5)/length(pred) # what proportion of true interactions are predicted as "likely" >0.5
sum(pred>0.75)/length(pred) # what proportion of true interactions are predicted as "v likely" >0.75


# Average and median probability of interaction based on the overall data
overall_mean <- cbind(apply(our_preds, 1, mean))
overall_median <- cbind(apply(our_preds, 1, median))

# Average and median in the held out data.
pred_mean <- apply(pred, 1, mean)
pred_median <- apply(pred, 1, median)


# Plot
plot_dta <- data.frame(value = rbind(pred_mean / overall_mean, pred_median / overall_median), 
                       stat = rep(c('mean_ratio', 'median_ratio'), repetitions))

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


