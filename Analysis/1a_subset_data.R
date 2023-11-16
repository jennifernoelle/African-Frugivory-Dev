# -------- TO DO --------- #

# Set the directories below to correspond to paths on your machine:

# The directory where the original data are:
data_path <- 'RawData/'
# The directory where the processed data should be saved:
save_path <- 'ProcessedData/'
# Whether the processed data should be saved or not:
save_files <- TRUE

# --------- BEGINNING -------- #


# Setting the working directory.
#setwd(wd_path)

# Loading libraries.
library(data.table)
library(gplots)
library(dplyr)
library(abind)

# ---------- PART A: Loading in the data -------------- #

# Loading in the data 
load(paste0(save_path, 'Adj_matrix_new_clean.dat'))
load(paste0(save_path, 'traits_verts_clean.dat'))
load(paste0(save_path, 'traits_plants_clean.dat'))
load(paste0(save_path, 'OA_full.dat'))
load(paste0(save_path, 'OP_fullclean.dat'))

fruit.lengths <- read.csv(paste0(data_path, "fruit_length_compilation3.csv"))


# ------ PART B: Creating Full Adjacency matrices ---------#
# The provided data only gives a subsets of the full adjacency matrix for each study
# Fill in non-edges with zeros

# Sort vertebrates and plants by family for ease of plotting later on
traits.verts <- traits.verts[order(traits.verts$Taxa, traits.verts$Family, traits.verts$Genus), ]
traits.plants <- traits.plants[order(traits.plants$Family, traits.plants$Genus), ]

vert.names <- traits.verts$Species
plant.names <- traits.plants$Species
study.names <- names(dta.clean)

# Get useful species numbers
nP <- length(plant.names)
nV <- length(vert.names)
nS <- length(dta.clean)

# Now create full adjacency matrix for each study using proper names
b <- length(dta.clean)
acomb <- function(...) abind(..., along=3)
registerDoParallel(10)

A.obs <- foreach(i=1:b, .combine = 'acomb') %dopar% {
  A <- dta.clean[[i]]
  wh.v <- sapply(dimnames(A)[[1]], function(x) which(x==vert.names))
  wh.p <- sapply(dimnames(A)[[2]], function(x) which(x==plant.names))
  
  A.full <- matrix(0, nrow = nV, ncol = nP)
  A.full[wh.v, wh.p] <- A
  A.full
}

dimnames(A.obs) <- list(vert.names, plant.names, study.names)


# ------ PART B: Subsetting the data to only mammal-plant interactions -------- #

# Get indices of mammals: we used the species order from the traits data
which.mammals <- which(traits.verts$Taxa=="Mammal")
traits.m <- traits.verts[which.mammals, ]

# Remove bird rows and all-bird studies
A.obs.m <- A.obs[which.mammals, , ]
which.bad.studies <- which(apply(A.obs.m, 3, sum)==0) 
A.obs.m <- A.obs.m[,,-which.bad.studies]

# Remove plants which are not consumed by any mammals: 711 of original 726 remain
# Also remove two conifers
which.bad.plants <- c(which(rowSums(apply(A.obs.m, 3, colSums))==0), which(colnames(A.obs.m)=="Juniperus_procera"), 
                      which(colnames(A.obs.m)=="Podocarpus_milanjianus"))
A.obs.m <- A.obs.m[, -which.bad.plants,]
traits.p.709 <- traits.plants[-which.bad.plants, ]

p.names <- traits.p.709$Species
m.names <- traits.m$Species
s.names <- unlist(dimnames(A.obs.m)[3])
s.ids <- as.numeric(gsub("Study.", "", s.names))

nM <- length(m.names)
nP <- length(p.names)
nS <- dim(A.obs.m)[3]

# Make sure names are in the right order
sum(p.names != colnames(A.obs.m))
sum(m.names != rownames(A.obs.m))

# ------ PART C: Preparing covariates for this subset -------- #

# Remove taxonomic info from traits data
# Put continuous covariates first
# Only binary and continuous covariates allowed
Obs_X <- traits.m %>%
         `rownames<-` (.[,1]) %>% # add species as rownames
         mutate(IUCN_Endangered = ifelse(IUCN_Status == "CR" | IUCN_Status=="EN", 1, 0), 
                Habitat_Forest = ifelse(Habitat == "Forest", 1, 0)) %>%
         select(c(Generation_Length, logBodyMass, logBrainToBodyMass,
                  IUCN_Endangered, Habitat_Forest))
Obs_W <- traits.p.709 %>%
          `rownames<-` (.[,1]) %>% 
          mutate(IUCN = ifelse(IUCN == "DD", NA, IUCN), 
                 IUCN_Endangered = ifelse(IUCN == "CR" | IUCN=="EN", 1, 0)) %>% 
          select(meanWD, IUCN_Endangered)
        # mutate(IUCN = as.factor(IUCN))

# Impute missing fruit lengths by family
fruit.lengths.df <- fruit.lengths %>% 
                        rename(Species=Plant_Species_Corrected, FruitLength= mFruitLength_sp) %>%
                        group_by(genus) %>% # First get genus mean fruit widths
                        mutate(GenusMean.FL = mean(FruitLength, na.rm = TRUE)) %>%
                        ungroup() %>%
                        group_by(family) %>% # Next get family mean fruit widths
                        mutate(FamilyMean.FL = mean(FruitLength, na.rm = TRUE)) %>%
                        mutate(FruitLength = ifelse(is.na(FruitLength), ifelse(is.na(GenusMean.FL), FamilyMean.FL, GenusMean.FL), FruitLength)) %>%
                        select(FruitLength, Species)
                        
fruit.lengths.ordered <- left_join(traits.p.709, fruit.lengths.df, by="Species")
Obs_W$FruitLength <- as.numeric(fruit.lengths.ordered$FruitLength)


if (save_files) {
  save(A.obs, file = paste0(save_path, 'obs_A_full.dat'))
  save(A.obs.m, file = paste0(save_path, 'obs_A_mammals.dat'))
  save(traits.p.709, file = paste0(save_path, 'traits_p_709.dat'))
  save(Obs_X, file = paste0(save_path, 'Obs_X.dat'))
  save(Obs_W, file = paste0(save_path, 'Obs_W_partiallycleaned.dat'))
}

# ------ PART D: Subsetting occurrence -------- #
# Camille has prepared external occurrence data which assumes that if a plant species is present 
# country A, in one study, it is always present in country A
# We are also going to prepare co-occurrence using Georgia's definition to double check

# Using site-level occurrences: just need to reorder and subset the data created in 0a
obs_OM <- O_A[rownames(O_A) %in% m.names, colnames(O_A) %in% s.names]
use_order_r <- sapply(m.names, function(r) which(rownames(obs_OM) == r))
use_order_c <- sapply(s.names, function(r) which(colnames(obs_OM) == r))
obs_OM <- obs_OM[use_order_r, use_order_c]

# Using site-level occurrences: just need to reorder and subset the data created in 0a
obs_OP <- O_P[rownames(O_P) %in% p.names, colnames(O_P) %in% s.names]
use_order_r <- sapply(p.names, function(r) which(rownames(obs_OP) == r))
use_order_c <- sapply(s.names, function(r) which(colnames(obs_OP) == r))
obs_OP <- obs_OP[use_order_r, use_order_c]

sum(rownames(obs_OM) != m.names)
sum(colnames(obs_OM) != s.names)
sum(rownames(obs_OP) != p.names)
sum(colnames(obs_OP) != s.names)



# ------ PART E: Preparing and subsetting focus -------- #

### USING DEFAULT FOCUS

# Creating two arrays that correspond to whether the studies would have
# recorded an observed interaction or not, and whether species co-occur in the
# study area. We will assume that these arrays are 0/1. For the focus array,
# this makes sense as we have access to whether the study is animal or plant
# focused or a network study. For the species occurrence, we know for a fact
# that a species with a recorded interaction occurs in the area. The only
# simplifying assumption is therefore when we set probability of occurrence to
# zero for species without any recorded interaction.


# We have to categorize studies first: all animal-focused except for maybe one (Study 13)
obs_F <- array(0, dim = c(nM, nP, nS))
dimnames(obs_F) <- list(m.names, p.names, s.names)
methods <- data.frame(study = s.ids, method = 'Animal-oriented') 

for (ss in 1 : nS) {
  this_study <- s.ids[ss]
  this_A <- A.obs.m[,,ss]
  these_mammals <- which(rowSums(this_A)>0)
  obs_F[these_mammals, , ss] <- 1 # because animal focused, all interactions would be recorded
}


if (save_files) {
  save(obs_OM, file = paste0(save_path, 'obs_OM.dat'))
  save(obs_OP, file = paste0(save_path, 'obs_OP.dat')) 
  save(obs_F, file = paste0(save_path, 'F_obs_default.dat'))
  save(traits.m, file = paste0(save_path, 'traits_m.dat'))
}
