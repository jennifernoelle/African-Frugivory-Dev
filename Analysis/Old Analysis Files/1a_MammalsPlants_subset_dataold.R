# -------- TO DO --------- #

# Set the directories below to correspond to paths on your machine:

# The directory where the analysis is performed:
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"

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
load(paste0(save_path, 'F_obs.dat'))
load(paste0(save_path, 'traits_verts_clean.dat'))
load(paste0(save_path, 'traits_plants_clean.dat'))

fruit.lengths <- read.csv(paste0(data_path, "fruit_length_compilation2.csv"))
obs_OP_df <- read.csv(paste0(data_path, "occurence3.csv")) # should add code for generating this eventually 


# ------ PART B: Creating Full Adjacency matrices ---------#
# The provided data only gives a subsets of the full adjacency matrix for each study
# Fill in non-edges with zeros

# Sort vertebrates and plants by family for ease of plotting later on
traits.verts <- traits.verts[order(traits.verts$Taxa, traits.verts$Family, traits.verts$Genus), ]
traits.plants <- traits.plants[order(traits.plants$Family, traits.plants$Genus), ]

vert.names <- traits.verts$Species
plant.names <- traits.plants$Species

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

dimnames(A.obs) <- list(vert.names, plant.names, paste0("Study.", 1:nS))


# ------ PART B: Subsetting the data to only mammal-plant interactions -------- #

# Get indices of mammals: we used the species order from the traits data
which.mammals <- which(traits.verts$Taxa=="Mammal")
traits.m <- traits.verts[which.mammals, ]
m.names <- traits.m$Species



# Remove bird rows and all-bird studies
A.obs.m <- A.obs[which.mammals, , ]
which.bad.studies <- which(apply(A.obs.m, 3, sum)==0) 
A.obs.m <- A.obs.m[,,-which.bad.studies]

# Remove plants which are not consumed by any mammals: 711 of original 726 remain
# Also remove two conifers
which.bad.plants <- which(rowSums(apply(A.obs.m, 3, colSums))==0)
which.bad.plants2 <- which(colnames(apply(A.obs.m, 3,colSums ))=="Podocarpus_milanjianus")
A.obs.m <- A.obs.m[,-c(which(colnames(A.obs.m)=="Juniperus_procera"),which(colnames(A.obs.m)=="Podocarpus_milanjianus")),]
A.obs.m <- A.obs.m[, -which.bad.plants, ]
traits.p.709 <- traits.plants[-which.bad.plants, ]
traits.p.709  <- traits.p.709[-c(which(traits.p.709$Species=="Juniperus_procera"),which(traits.p.709$Species=="Podocarpus_milanjianus")),]
p.names <- traits.p.709$Species

# Make sure names are in the right order
sum(p.names != colnames(A.obs.m))
sum(m.names != rownames(A.obs.m))

nM <- length(m.names)
nP <- length(p.names)
nS <- dim(A.obs.m)[3]


# ------ PART C: Preparing covariates for this subset -------- #

# Remove taxonomic info from traits data
# Put continuous covariates first
# Only binary and continuous covariates allowed
Obs_X <- traits.m %>%
         `rownames<-` (.[,1]) %>% # add species as rownames
         mutate(IUCN_Endangered = ifelse(IUCN_Status == "CR" | IUCN_Status=="EN", 1, 0), 
                Habitat_Forest = ifelse(Habitat == "Forest", 1, 0)) %>%
         select(c(Generation_Length, logBodyMass, logBrainMass, 
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
                        group_by(family) %>% # First get family mean fruit widths
                        mutate(FamilyMean.FL = mean(FruitLength, na.rm = TRUE)) %>%
                        mutate(FruitLength = ifelse(is.na(FruitLength), FamilyMean.FL, FruitLength)) %>%
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

# ------ PART D: Preparing and subsetting occurrence -------- #

## For occurrence of network in study
obs_OM <- array(0, dim = c(nM, nS))
obs_OP <- array(0, dim = c(nP, nS))
dimnames(obs_OM) <- list(m.names, paste0("study", 1:nS))
dimnames(obs_OP) <- list(p.names, paste0("study.", 1:nS))

## Check that species name order matches Adjacency matrix order
sum(row.names(A.obs.m) != m.names)
sum(colnames(A.obs.m) != p.names)

# Define occurrence based on whether or not the species was observed in the study
for (ss in 1 : nS) {
  these_mammals <- which(rowSums(A.obs.m[,,ss])!=0)
  these_plants <- which(colSums(A.obs.m[,,ss])!=0)
  obs_OM[these_mammals, ss] <- 1
  obs_OP[these_plants, ss] <- 1
}

## Combine the occurrence matrix above with one using expert knowledge in loaded data (obs_OP_df)
# Make the occurrence matrix to be based on herbivory interactions/ country-level data 
rownames(obs_OP_df) <-obs_OP_df$X
obs_OP_df <- obs_OP_df[,-1]
obs_OP_df1 <- obs_OP_df
obs_OP2 <- t(data.frame(matrix(unlist(obs_OP_df), ncol = max(lengths(obs_OP_df)), byrow = TRUE)))
rownames(obs_OP2) <- rownames(obs_OP_df1)
obs_OP2 <- obs_OP2[-c(which(rownames(obs_OP2)=="Juniperus_procera"),which(rownames(obs_OP2)=="Podocarpus_milanjianus")),]

# Check for inconsistencies: if a plant was included in a study, it must be present in the area
which(obs_OP2 == 0 & obs_OP>0, arr.ind = TRUE)
obs_OP2 <- ifelse(obs_OP2 == 0 & obs_OP > 0, obs_OP, obs_OP2)
obs_OP <- obs_OP2




# ------ PART E: Preparing and subsetting focus -------- #

# Subset
F_obs_m <- F_obs[which(rownames(F_obs) %in% m.names), which(colnames(F_obs) %in% p.names), 
                 -which.bad.studies]
# Reorder focus
use_order_r <- sapply(m.names, function(r) which(rownames(F_obs_m) == r))
use_order_c <- sapply(p.names, function(r) which(colnames(F_obs_m) == r))

F_obs_m <- F_obs_m[use_order_r, use_order_c,]



if (save_files) {
  save(obs_OM, file = paste0(save_path, 'obs_OM.dat'))
  save(obs_OP, file = paste0(save_path, 'obs_OP.dat'))
  save(F_obs_m, file = paste0(save_path, 'F_obs_m.dat'))
  save(traits.m, file = paste0(save_path, 'traits_m.dat'))
}
