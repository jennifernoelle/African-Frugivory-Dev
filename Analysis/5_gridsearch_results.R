# This views the results of a  grid search for reasonable occurrence probabilities
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
library(parallel)
library(foreach)
library(abind)
library(magrittr)
library(pheatmap)


# Loading the data:
load(paste0(data_path, 'obs_A_mammals.dat'))
load(paste0(data_path, 'Subset_no_baboons.dat')) # list with appropriately subset mammals, plants, studies

# Define a useful function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


## Rename for convenience
obs_A <- A.obs.m

## Subset out baboons and associated plants, studies
wh_keep_m <- which(rownames(obs_A) %in% nobab.list[[1]])
wh_keep_p <- which(colnames(obs_A) %in% nobab.list[[2]])
wh_keep_s <- which(unlist(dimnames(obs_A)[3]) %in% nobab.list[[3]])

obs_A <- obs_A[wh_keep_m, wh_keep_p, wh_keep_s]


# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

nM <- nrow(obs_A)
nP <- ncol(obs_A)

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
n.cv <- 100
repetitions <- 4

# ---------------- STEP 3: Load results -------------- #
res.df <- data.frame(matrix(NA, ncol = 1 + ncol(prob.grid) + 5, nrow = nrow(prob.grid)))
colnames(res.df) <- c("p.set", paste0("p", 1:7), "mean.ratio", "median.ratio", "cv.mean", "cvp50", "cvp75")

for(x in 1:nrow(prob.grid)){
  cat("\n Probability scenario: ", x)
  p <- prob.grid[x,]
  all_indices <- array(NA, dim = c(repetitions, n.cv, 2))
  our_preds <- array(NA, dim = c(repetitions, nM, nP)) # all preds
  pred <- array(NA, dim = c(repetitions, n.cv)) # cv preds

  for (rr in 1 : repetitions) {
    cat("\n Chain = ", rr)
    cv_indices <- loadRData(paste0(save_path, 'GridSearch/cv_indices_probs_', x, '_', rr, '.dat'))
    all_preds <- loadRData(paste0(save_path, 'GridSearch/pred_probs_', x, '_', rr, '.dat'))
    all_indices[rr, , ] <- cv_indices
    our_preds[rr, , ] <- all_preds
    }
    
    # Predictions of the held out data from model: 100 indices held out each time - always interactions
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
  
  # Performance metrics
  mean_ratio <- mean(pred_mean/overall_mean)
  median_ratio <- mean(pred_median/overall_median)
  cv_mean <- mean(pred_mean)
  p50 <- sum(pred>0.5)/length(pred) # what proportion of true interactions are predicted as "likely" >0.5
  p75 <- sum(pred>0.75)/length(pred) # what proportion of true interactions are predicted as "v likely" >0.75
  
  mean_perf <- c(mean_ratio, median_ratio, cv_mean, p50, p75) # mean perf across cv repetitions
  res.df[x,] <- c(x, p, mean_perf) # Include the current prob just in case they get jumbled
}


write.csv(res.df, paste0(save_path, "GridSearch/results_summary.csv"), row.names = FALSE)
