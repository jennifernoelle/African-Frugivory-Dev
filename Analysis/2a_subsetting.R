# Perform additional subsetting prior to running the mcmc

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
result_path <- 'Results/'
# Where the functions are available:
source_path <- 'HelperScriptsNew/'


library(tidyverse)

# Loading the data:
load(paste0(data_path, 'obs_A_mammals.dat'))
load(paste0(data_path, 'F_obs_default.dat'))
load(paste0(data_path, "study_info1.dat")) # useful for subsetting

## Rename for convenience
obs_A <- A.obs.m

m.names <- rownames(obs_A)
p.names <- colnames(obs_A)
s.names <- unlist(dimnames(obs_A)[3])

## Subset to remove baboons
# Select baboons from species list
baboons <- m.names[which(grepl("Papio", m.names))]
good.mammals <- m.names[-which(grepl("Papio", m.names))]
wh_keep_m <- which(rownames(obs_A) %in% good.mammals)

# See if baboon studies cover other species: they don't
baboon.studies <- study_info1 %>% filter(Species %in% baboons) %>% pull(Study_ID)
study_info1 %>% filter(Study_ID %in% baboon.studies)
good.studies <- s.names[!(s.names %in% paste0("Study.", baboon.studies))]
wh_keep_s <- which(unlist(dimnames(obs_A)[3]) %in% good.studies)

# Remove plants which only occur in the baboon studies
baboon.plants <- colnames(obs_A)[which(apply(obs_A[,, wh_keep_s], 2, sum) ==0)]
wh_keep_p <- which(!(p.names %in% baboon.plants))
good.plants <- p.names[wh_keep_p]

nobab.list <- list(good.mammals, good.plants, good.studies )

save(nobab.list, file = paste0(data_path, "Subset_no_baboons.dat"))


