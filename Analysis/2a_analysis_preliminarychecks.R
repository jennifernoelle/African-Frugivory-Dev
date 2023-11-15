# This file just looks at all of the processed data and makes sure they look reasonable


# --------- TO DO: set  your directories and name the current results files using the date---------#
## YOU REALLY HAVE TO SET THESE ##

# The directory where the analysis is performed:
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"
setwd(wd_path)

## THESE SHOULD BE THE SAME IF YOU CLONED THE REPO ##

# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
save_path <- 'Results/'
# Where the functions are available:
source_path <- 'HelperScriptsNew/'

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
#load(paste0(data_path, 'obs_OP_defaultMethod.dat')) # study level obs_OP2
#load(paste0(data_path, 'obs_OM_defaultMethod.dat')) # study level obs_OM2
load(paste0(data_path, 'obs_OM.dat')) # site level obs_OM
load(paste0(data_path, 'obs_OP.dat')) # site level obs_OP
load(paste0(data_path, 'traits_p_709_clean.dat')) # plant occurrences: most species have many studies, but some have 0?


## Rename for convenience
Cu <- Cu_phylo
Cv <- Cv_phylo
obs_A <- A.obs.m

## Useful quantities: sample sizes of the two sets of species
nM <- nrow(Cu)
nP <- nrow(Cv)
nStudies <- dim(obs_A)[3]

m.names <- rownames(obs_A)
p.names <- colnames(obs_A)
s.names <- unlist(dimnames(obs_A)[3])


# The combined network:
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1 # 1 if any interactions observed

## Sanity checks

# Construct the matrix F_ijs * O_ijs using various constructions of occurrence
O_ijs <- array(data = NA, dim = c(nM, nP, nStudies))

# Use binary version of occurrence probs such that they are present if at site, absent otherwise
obs_OPb <- ifelse(obs_OP > 0.5, 1, 0)
obs_OMb <- ifelse(obs_OM > 0.5, 1, 0)

dimnames(O_ijs) <- list(rownames(obs_A), colnames(obs_A), paste0("Study_", 1:nStudies))
for(s in 1:nStudies){
  for(m in 1:nM){
    for(p in 1:nP){
      O_ijs[m,p,s] <- obs_OMb[m,s] * obs_OPb[p,s] 
    }
  }
}


# Interactions don't occur with non-observable species
mean(obs_A[O_ijs==0])
#mean(obs_A[O_ijs2 == 0])
mean(obs_A[obs_F==0])

# Interaction prevalence is very low, but higher when conditioning on observability (O*F)
# Using default observabiilty leads to very very high prevalence, 
# But this makes sense because most studies are single-animal focused, 
# So if a plant is recorded in that study it's usually because the focus animal ate it
mean(obs_A[O_ijs==1])
mean(obs_A[obs_F==1]) 

# Most interactions aren't observable,evenly due to occurrence and focus
mean(O_ijs)
mean(obs_F)

# All species occur at least once
summary(rowSums(obs_OPb))
summary(rowSums(obs_OMb))

# Look at mammal occurrence when Focus is 1
F_m <- apply(obs_F, 3, rowSums)
F_m <- ifelse(F_m>0, 1, 0)
mean(obs_OM[F_m ==1]) # animals are always present at the study level when they are in focus

# Look at distribution of plant occurrences
O_P_f <- factor(O_P, labels = c("Recorded in region",  "Recorded in country", 
                                "Recorded in habitat but not country", 
                                "Recorded in region and habitat", "Recorded in country and habitat",
                                "Recorded at site", "Recorded in study"))
par(las = 2)
par(mar=c(5,10,4,2)) # increase y-axis margin.
barplot(prop.table(table(O_P_f)),
        xlab = "Proportion",
        horiz = TRUE, 
        main = "Plant Occurrences")

# How many studies for each mammal species
apply(obs_A, 1, function(x) sum(colSums(x)>0))

# How many species per study
apply(obs_A, 3, function(x) sum(rowSums(x)>0))
which(apply(obs_A, 3, function(x) sum(rowSums(x)>0))==0)

# Look at plant and mammal phylogenetic correlations
summary(c(Cu)) # mammals are more closely related
summary(c(Cv)) # we have more plants diversity and lower mean plant phylo correlation

# Look at traits
apply(Obs_X, 2, function(x) sum(is.na(x)))
apply(Obs_W, 2, function(x) sum(is.na(x)))

# Make sure names are all consistent
sum(row.names(Obs_X) != m.names)
sum(row.names(Obs_W) != p.names)
sum(row.names(obs_OM) != m.names)
sum(row.names(obs_OP) != p.names)
sum(row.names(Cu) != m.names)
sum(row.names(Cv) != p.names)
sum(colnames(obs_OP) != s.names)
sum(colnames(obs_OM) != s.names)

# Which species interact the most
mean.int.0 <- round(rowSums(comb_A)/nP, 4)
mean.int <- rep(NA, nM)
for(m in 1:nM){
  this.a <- obs_A[m,, ]
  this.of <- (O_ijs * obs_F)[m,, ]
  mean.int[m] <- mean(this.a[this.of ==1], na.rm = TRUE)
}

prev.0 <- data.frame(Mammal = m.names , Prevalence.0 = mean.int.0) 
prev <- data.frame(Mammal = m.names , Prevalence = mean.int) %>% 
        arrange(-Prevalence) %>% 
        left_join(., prev.0) %>% 
        pivot_longer(!Mammal, names_to = "Type", values_to = "Prevalence") %>% 
        mutate(row = row_number())  %>% # create row to preserve ordering
        mutate(Mammal = gsub("_", " ", Mammal))

ggplot(data = prev, aes(x = reorder(Mammal, row), y = Prevalence, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  labs(title = "Interaction prevalence in the data") + 
       xlab("") + 
  scale_fill_discrete(labels = c("Conditional", "Raw"), name = NULL) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90), text = element_text(family = "serif", size = 20)) 

## Investigate isolated sites

# Idea: for each study, look at the proportion of plants that are found only at that site
# And how many studies do we have for that same site?
# I would have to go back to the raw data to get that, good to do, save as a separate file

study1 <- obs_OP[,1]
sum(study1==1) # Observed in study in interactions
sum(study1==0.75) # Observed at same site in interactions
sum(study1==0.5) # Observed in same habitat and country

obs_to_site = data.frame(study = s.names, 
                    obs_to_site = apply(obs_OP, 2, function(x) sum(x == 1))/apply(obs_OP, 2, function(x) sum(x > 0.5)))
obs_to_hab = data.frame(study = s.names, 
                        obs_to_hab = apply(obs_OP, 2, function(x) sum(x == 1))/apply(obs_OP, 2, function(x) sum(x > 0)))

study.plants <- data.frame(study = s.names, 
        obs_to_habcountry = apply(obs_OP, 2, function(x) sum(x == 1))/apply(obs_OP, 2, function(x) sum(x >0.25)))%>% 
        arrange(-obs_to_habcountry) %>% 
        left_join(., obs_to_site) %>% 
        left_join(., obs_to_hab) %>% 
        pivot_longer(!study, names_to = "Type", values_to = "Ratio") %>%
        mutate(row = row_number())  

ggplot(study.plants, aes(x = reorder(study, row), y = Ratio, fill = Type)) +
      geom_bar(stat = "identity", position = position_dodge()) + 
      labs(title = "Proportion unique plants") + 
      xlab("") + 
      scale_fill_discrete(labels = c("Habitat level", "Habitat-country level", "Site level"), name = NULL) +
      theme_minimal() + 
      theme(text = element_text(family = "serif", size = 20), axis.text.x = element_text(angle = 90, size = 15)) 

# What is the focus of the dubious studies
# Should still add to this how many studies at the iffy sites
isol.studies <- study.plants %>% 
                filter(Type == "obs_to_hab", Ratio ==1) %>% 
                pull(study)
                
F.isol <- obs_F[,,which(unlist(dimnames(obs_F)[3]) %in% isol.studies)]
m.names[apply(F.isol, 3, function(x) which(rowSums(x)> 0))]


apply(F.isol,3, function(x) rowSums(x))



## Problematic entries
# Erythrocebus patas only has one interaction observed and one study with limited plant overlap
# Look at most interactive species
hungry.m <- order(mean.int, decreasing = TRUE)[1:5]
m.names[hungry.m]

baboons <- hungry.m[1:2]
baboon.studies <- unique(which(obs_F[baboons,,]==1, arr.ind = TRUE)[,3])

# For species with a very high conditional interaction prevalence
# Investigate distribution of plants: plants with no overlap in other sites cannot
# be reliablly imputed
# No "In habitat and country" plants
table(obs_OP[, baboon.studies])
O_P_baboon <- factor(obs_OP[, baboon.studies], labels = c("Not in habitat", "In habitat, not country", 
                                      "At site", "In study"))
barplot(table(O_P_baboon), horiz = TRUE, main = "Plants at baboon sites")
baboon.plants <- unique(which(obs_OPb))

## Investigating isolated plants
# Idea: compute the # of sites, etc a plant is present at

barplot(table(O_P_f), horiz = TRUE, main = "Level of occurrence info for plants in the corpus")

# How many plants have site-level occurrence information?
ifelse(obs_OP== 0.75, 1, 0) %>% 
  rowSums(.) %>% 
  data.frame(times = .) %>% 
  ggplot(.) + geom_histogram(aes(x = times))

  hist(., xlab = "Number of times plant has site-level occurrence", main = "Site level occurrence info", 
       breaks = 21)

#--------------------- Repeat sanity checks for a subset of species ------------------#

## Subsetting - Issue: one Erythrocebus patas only eats Celtis toka and nothing else eats Celtis_toka
# Remove that study too
good.studies <- which(apply(obs_A, 3, sum)>1)
common.plants <- names(which(apply(obs_A[,,good.studies], 2, sum)>5))
good.mammals <- names(which(apply(obs_A[,,good.studies], 1, sum)>1))

obs_A_common <- obs_A[unlist(dimnames(obs_A)[1]) %in% good.mammals, unlist(dimnames(obs_A)[2]) %in% common.plants, good.studies ]
F_obs_common <- obs_F[rownames(obs_OM) %in% good.mammals, colnames(obs_F) %in% common.plants, good.studies]
obs_OP_common <- obs_OP[rownames(obs_OP) %in% common.plants, good.studies ]
obs_OM_common <- obs_OM[rownames(obs_OM) %in% good.mammals, good.studies]

# The combined network:
comb_A_common <- apply(obs_A_common, c(1, 2), sum)
comb_A_common <- (comb_A_common > 0) * 1 # 1 if any interactions observed

# Useful quantities
nM.common <- nrow(obs_A_common)
nP.common <- ncol(obs_A_common)
nS.common <- dim(obs_A_common)[3]

## Sanity checks for this subset of common plants

# Construct the matrix F_ijs * O_ijs using various constructions of occurrence
O_ijs_common <- array(data = NA, dim = c(nM.common, nP.common, nS.common))
dimnames(O_ijs_common) <- list(rownames(obs_A_common), colnames(obs_A_common), paste0("Study_", 1:nS.common))
for(s in 1:nS.common){
  for(m in 1:nM.common){
    for(p in 1:nP.common){
      O_ijs_common[m,p,s] <- obs_OM_common[m,s] * obs_OP_common[p,s] 
    }
  }
}

# Interaction prevalence is very low, but higher when conditioning on observability (O*F)
# Using default observabiilty leads to very very high prevalence, 
# But this makes sense because most studies are single-animal focused, 
# So if a plant is recorded in that study it's usually because the focus animal ate it
mean(obs_A_common[O_ijs_common==1])
mean(obs_A_common[F_obs_common==1]) # this is the effective prevalence when we set occurrence = 1 always
OF_common <- O_ijs_common * F_obs_common
mean(obs_A_common[OF_common==1])

# Interpretation: our interaction prevalence is around 0.18 among observable interactions, not too high
# However, the model is also learning detection probabilities, so my guess is that when we run with these levels
# the model is learning relatively low detection probabilities and hence a high number of unobserved interactions

# Most interactions aren't observable,evenly due to occurrence and focus
mean(O_ijs_common)
mean(F_obs_common)

good.plants <- common.plants
good.studies <- names(good.studies)

# Save the species and studies used in this subset
save(file = paste0(data_path, "subset1_studies.dat"), good.studies)
save(file = paste0(data_path, "subset1_mammals.dat"), good.mammals)
save(file = paste0(data_path, "subset1_plants.dat"), good.plants)


