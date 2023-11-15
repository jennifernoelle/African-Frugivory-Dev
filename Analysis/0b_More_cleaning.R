# This file cleans the processed data from 0a and the trait data
# And creates the following inputs for the function
# Adj_matrix_new_clean.dat
# traits_verts_clean.dat
# traits_plants_clean.dat

# -------- TO DO --------- #

# Set the directories below to correspond to paths on your machine:

# The directory where the analysis is performed:
#wd_path <- "/Users/camilledesisto/Documents/GitHub/African-Frugivory"
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"

# The directory within the working directory where the data are:
data_path <- 'RawData/'
save_path <- 'ProcessedData/'

library(foreach)
library(doParallel)
library(data.table)
library(ggmap)
library(ggplot2)
library(tidyverse)
library(forcats)
library(abind)

#setwd(wd_path)
save_files <- TRUE

# ---------- PART A: Loading in the data -------------- #

# Loading in the data sets

dta <- readRDS(paste0(save_path, 'Adj_matrix_new.rds')) # created in 0a
dta.summary <- read.csv(paste0(data_path, 'frugivore_mutualism_interactions.csv'))
dist.plants.sp <- read.csv(paste0(data_path, 'distance.matrix_species_level.csv'))
dist.plants.g <- read.csv(paste0(data_path, 'distance.matrix_genus_level.csv'))
traits.plants <- read.csv(paste0(data_path, 'plant_traits_mut.csv'))[, -1]
traits.verts <- read.csv(paste0(data_path, 'animal_traits_mut.csv'))[, -1]
# load(paste0(save_path, 'OA_full.dat'))
# load(paste0(save_path, 'OP_full.dat'))
load(paste0(save_path, 'OA_full.dat'))
load(paste0(save_path, 'OP_full.dat'))

# Clean the data
row.names(dta.summary) <- dta.summary$X
dta.summary <- dta.summary[, -1]

# --------- PART B: Cleaning the names --------- # 

### Check matrix dims: several are row/matrices, one is a scalar (missing species names) Study 19
sapply(dta, function(x) dim(x))
dta[[16]]
rownames(dta[[16]])<-"Erythrocebus_patas"
colnames(dta[[16]]) <- "Celtis_toka"

### Make sure species names are consistent across adjacency matrices, traits

# Create list of vertebrate and plant names
vert.names <- sort(unique(unlist(lapply(dta, rownames ))))
plant.names.temp <- unique(unlist(lapply(dta, colnames)))

# Plant names need to be cleaned: remove hyphens and double subscripts
plant.names.clean <- gsub(" ", "", plant.names.temp)
plant.names.clean <- gsub("-", "", plant.names.clean)
plant.names.clean <- ifelse(plant.names.clean == "Vachellia xanthophloea", "Vachellia_xanthophloea", plant.names.clean)

fix.names <- function(y){
  wh.p <- sapply(dimnames(y)[[2]], function(x) which(x == plant.names.temp))
  dimnames(y)[[2]] <- plant.names.clean[wh.p]
  return(y)
}
dta.clean <- sapply(dta, fix.names)

#F.clean <- lapply(F_ijs, fix.names) # matrices same dim so sapply will vectorize
plant.names <- sort(unique(plant.names.clean))
rm(plant.names.temp)

# Clean names in plant occurrences
rownames(O_P)[!(rownames(O_P) %in% plant.names.clean)]
plant.names.clean[!(plant.names.clean) %in% rownames(O_P)]
rownames(O_P) <- gsub(" ", "", rownames(O_P))
rownames(O_P) <- gsub("-", "", rownames(O_P))
rownames(O_P) <- ifelse(rownames(O_P) == "Vachellia xanthophloea", "Vachellia_xanthophloea", rownames(O_P))

# Check for consistency with naming in traits data
sum(!(traits.verts$Species %in% vert.names))
sum(!(traits.plants$Species %in% plant.names))
traits.plants$Species[which(!(traits.plants$Species %in% plant.names))]
traits.plants$Species <- gsub(" ", "", traits.plants$Species)
traits.plants$Species <- gsub("-", "", traits.plants$Species)
sum(!(traits.plants$Species %in% plant.names))

# Check for duplicate plants in traits file: keep entry with more complete info
which(table(traits.plants$Species)>1)
traits.plants[traits.plants$Species == "Diospyros_crassiflora", ]
traits.plants <- traits.plants[-396, ]

# Look for unit errors in body mass (grams) and generation length (days)
traits.verts <- traits.verts %>% mutate(Body_Mass = ifelse(Body_Mass < 10, Body_Mass * 1000, Body_Mass)) %>% 
                 mutate(Generation_Length = ifelse(Generation_Length < 100, Generation_Length * 365, Generation_Length)) %>% 
                 mutate(Body_Mass = ifelse(Species == "Heliosciurus_rufobrachium", 355, Body_Mass))
                

# Misc cleaning
traits.verts <- traits.verts %>% 
                mutate(logBodyMass = log(Body_Mass), logBrainToBodyMass = log(Brain_Mass/Body_Mass))


# Save new and cleaned files
if (save_files) {
  save(dta.clean, file = paste0(save_path, 'Adj_matrix_new_clean.dat'))
  save(traits.verts, file = paste0(save_path, 'traits_verts_clean.dat'))
  save(traits.plants, file = paste0(save_path, 'traits_plants_clean.dat'))
  save(O_P, file = paste0(save_path, 'OP_fullclean.dat'))
  #save(F_obs, file = paste0(save_path, 'F_obs.dat'))  
}

