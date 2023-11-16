# This file processes the raw data to create the following key data sources: 
# Adj_matrix_new.rds, O_A, O_P
# Saves useful info about studies

# -------- TO DO --------- #

# Set the directories below to correspond to paths on your machine:

# The directory within the working directory where the data are:
data_path <- 'RawData/'
save_path <- 'ProcessedData/' # Create this directory

library(foreach)
library(doParallel)
library(data.table)
library(ggmap)
library(ggplot2)
library(tidyverse)
library(forcats)
library(abind)
library(readxl)

save_files <- TRUE

# ---------- PART A: Loading in the data -------------- #

frugivore_mut  <- read.csv(paste0(data_path, "Frugivory.csv"))
frug_fol_mut <- read.csv(paste0(data_path, "Frugivory_folivory.csv"))[, -c(1:2)] # frugivory and folivory data combined: use this for occurrences
study_dta <- read.csv(paste0(data_path, "Site_metadata.csv")) %>% 
                      mutate(Country = ifelse(Country == "Republic of Congo", "Congo", Country))
study_locs <- read_excel(paste0(data_path, "Study_locations.xlsx"))

# ------------- PART B: CLEAN STUDY META DATA ---------------------#

## EDA: check for multi-site studies in the combined data
frugivore_mut %>% group_by(Study_ID) %>% 
  summarize(n_Sites = length(unique(Study_Site))) %>% 
  filter(n_Sites > 1)

## Clean site names in frugivory and folivory data and make sure they are consistent with study metadata

# Frugivory + folivory
frug_fol_mut <- frug_fol_mut %>% 
  mutate(Study_Site = ifelse(Study_Site == "Ivindo", "Makokou", Study_Site)) %>%
  mutate(Study_Site = ifelse(Study_Site == "Lopé Reserve" | Study_Site == "Lope Reserve" | Study_Site == "Lope reserve", 
         "Lope Reserve", Study_Site)) %>%
  mutate(Study_Site = ifelse(Study_ID ==80, "Nyungwe National Park", Study_Site)) %>%
  mutate(Study_ID = ifelse(Study_ID == 78 & Study_Site == "Makokou", 78.1, Study_ID)) %>% 
  mutate(Study_ID = ifelse(Study_ID ==78 & Study_Site == "Wonga Wongue", 78.2, Study_ID)) %>% 
  mutate(Study_Site = ifelse(Study_Site == "Tiwai Island", "Tiwai", Study_Site)) %>% 
  mutate(Study_Site = ifelse(Study_Site == "Mingingi" | Study_Site == "Mingingi ", "Mingingi-Moke", Study_Site)) %>%
  mutate(Study_Site = ifelse(Study_Site == "Dja" | Study_Site == "Dja Biosphere Reserve" | 
                               Study_Site == "Bouamir Research Station, Dja Reserve", "Dja Reserve", Study_Site)) %>%
  mutate(Study_Site = ifelse(Study_Site == "Kibale Forest Reserve", "Kibale National Park", Study_Site)) %>%
  mutate(Study_Site = ifelse(Study_Site == "Budongo Forest Reserve", "Budongo Forest", Study_Site)) %>%
  mutate(Study_Site = ifelse(Study_Site == "Comoé National Park", "Comoe National Park", Study_Site)) %>%
  mutate(Study_Site = ifelse(Study_Site == "Djaloumbe ", "Djaloumbe", Study_Site)) %>%
  mutate(Study_Site = ifelse(Study_Site == "Issa Valley ", "Issa Valley", Study_Site)) %>%
  mutate(Study_Site = ifelse(Study_Site == "Djangui" | Study_Site == "Ndangaye" | Study_Site == "Petite Savane" |
                               Study_Site == "Djaloumbe", "Lobeke National Park", Study_Site)) %>% 
  mutate(Study_Site = ifelse(Study_Site == "Mabale" | Study_Site == "Guga" | Study_Site == "Bakussa" |
                               Study_Site == "Bonye 2" | Study_Site == "Mingingi-Moke" |
                               Study_Site =="Mokele" | Study_Site == "Nouabalé-Ndoki National Park",  
                              "Nouabale-Ndoki National Park", Study_Site)) %>% 
  mutate(Study_Site = ifelse(Study_Site == "Mondika Research Center", "Nouabale-Ndoki National Park", Study_Site)) %>%
  mutate(Plant_Species_Corrected = ifelse(Plant_Species_Corrected == "Diospyros_crassiflora ", "Diospyros_crassiflora", Plant_Species_Corrected))


# Clean names in frugivory data and create unique identifer for one multisite study
# Note there are other multi-site studies in the combined data, but we will not be matching on study id between these two data sets
# So we leave as is above
frugivore_mut <- frugivore_mut %>% mutate(Study_Site = ifelse(Study_Site == "Ivindo", "Makokou", Study_Site)) %>% 
  mutate(Study_Site = ifelse(Study_Site == "Bouamir Research Station, Dja Reserve" | 
                               Study_Site == "Dja Biosphere Reserve", "Dja Reserve", Study_Site)) %>% 
  mutate(Study_Site = ifelse(Study_Site == "Lopé Reserve" | Study_Site == "Lope reserve", "Lope Reserve", Study_Site)) %>% 
  mutate(Study_Site = ifelse(Study_Site == "Budongo Forest Reserve", "Budongo Forest", Study_Site)) %>% 
  mutate(Study_Site = ifelse(Study_Site == "Tiwai Island", "Tiwai", Study_Site)) %>% 
  mutate(Study_Site = ifelse(Study_Site == "Kibale Forest Reserve", "Kibale National Park", Study_Site)) %>%
  mutate(Study_ID = ifelse(Study_ID == 78 & Study_Site == "Makokou", 78.1, Study_ID)) %>% 
  mutate(Study_ID = ifelse(Study_ID ==78 & Study_Site == "Wonga Wongue", 78.2, Study_ID)) %>% 
  mutate(Study_Site = ifelse(Study_Site == "Djangui" | Study_Site == "Ndangaye" | Study_Site == "Petite Savane" |
                               Study_Site == "Djaloumbe", "Lobeke National Park", Study_Site)) %>% 
  mutate(Study_Site = ifelse(Study_Site == "Mabale" | Study_Site == "Guga" | Study_Site == "Bakussa" |
                               Study_Site == "Bonye 2" | Study_Site == "Mingingi-Moke" |
                               Study_Site =="Mokele" | Study_Site == "Nouabalé-Ndoki National Park", 
                             "Nouabale-Ndoki National Park", Study_Site)) %>% 
  mutate(Study_Site = ifelse(Study_Site == "Mondika Research Center", "Nouabale-Ndoki National Park", Study_Site)) %>%
  mutate(Plant_Species_Corrected = ifelse(Plant_Species_Corrected == "Diospyros_crassiflora ", "Diospyros_crassiflora", Plant_Species_Corrected))

# Make sure site names are consistent
unique(frugivore_mut$Study_Site)[!(unique(frugivore_mut$Study_Site) %in% unique(frug_fol_mut$Study_Site))]
unique(frugivore_mut$Study_Site)[!(unique(frugivore_mut$Study_Site) %in% unique(study_dta$Study_Site))]

# Add country and habitat information to frugivory-folivory data
frug_fol_mut <- left_join(frug_fol_mut, study_dta, by = c('Study_Site', 'Study_ID')) %>%
                mutate(Region = ifelse(Country == "Cameroon" | Country == "Central African Republic" | Country == "Congo" | 
                            Country == "DRC" | Country == "Gabon", "Central Africa", 
                           ifelse(Country == "Ethiopia" | Country == "Tanzania"  | Country == "Kenya" | Country =="Uganda" | 
                           Country == "Rwanda", "East Africa", "West Africa")))

# See if any sites are missing habitat info
unique(frug_fol_mut[is.na(frug_fol_mut$Habitat), c('Study_Site', 'Study_ID')])

# Extract plant and animal list for frugivory only data
animal_list <- frugivore_mut %>%
  select(Animal_Species_Corrected)%>%
  unique()%>%
  filter(Animal_Species_Corrected!="") %>%
  pull(Animal_Species_Corrected)

plant_list <- frugivore_mut %>%
  select(Plant_Species_Corrected)%>%
  unique() %>%
  pull(Plant_Species_Corrected)

study_list <- frugivore_mut %>%
  select(Study_ID)%>%
  unique() %>%
  pull(Study_ID)

study_names <- paste0("Study.", study_list)

nA <- length(animal_list)
nP <- length(plant_list)
nS <- length(study_list)

# How many times is each study site repeated? 
unique_locs <- frug_fol_mut %>% select(Study_ID, Study_Site, Region, Country) %>% 
                 unique() %>% 
                 left_join(., study_locs, by = "Study_ID")  %>% 
                 group_by(Study_Site, Region, Country) %>% 
                 summarize(n_studies = n(), Lat = mean(Lat, na.rm = TRUE), Long = mean(Long, na.rm = TRUE)) %>% 
                 filter(Study_Site %in% unique(frugivore_mut$Study_Site)) %>% 
                 as.data.frame(.) %>% 
                 mutate(Lat = ifelse(Study_Site == "Wonga Wongue", 0.5186, Lat), 
                        Long = ifelse(Study_Site == "Wonga Wongue", 9.4562, Long), 
                        Lat = ifelse(Study_Site == "Kakamega Forest", 0.2913, Lat), 
                        Long = ifelse(Study_Site == "Kakamega Forest", 34.8565, Long)) 

## Create data frame summarizing study information
# How many times is each study site repeated? 
animal_families <- frugivore_mut %>% 
                   select(Animal_Species_Corrected, Animal_Family) %>% 
                   unique() %>% 
                   mutate(Animal_Family = ifelse(Animal_Family=="Cercopithecidae ", "Cercopithecidae", Animal_Family))
study_focus1 <- data.frame(Study_ID = double(), nfocus = integer(), Species = character())
study_focus2 <- data.frame(Study_ID = double(), nfocus = integer(), Family1 = character(), Family2 = character(),
                           Family3 = character(), 
                           Species1 = character(), 
                           Species2 = character(), Species3 = character(), Species4 = character(), 
                           Species5 = character(), Species6 = character(), Species7 = character())

counter1 <- 1
counter2 <- 1
for(s in 1:nS){
  this_study <- study_list[s]
  these_obs <- frugivore_mut[frugivore_mut$Study_ID == this_study, ]
  these_animals <- unique(these_obs$Animal_Species_Corrected)
  these_families <- unique(these_obs$Animal_Family)
  for(a in 1:length(these_animals)){
      study_focus1[counter1,] <- c(this_study, length(these_animals), these_animals[a])
      counter1 <- counter1 + 1
    }
  study_focus2[counter2, ] <- c(this_study, length(these_animals), these_families, rep(NA, 3 - length(these_families)), 
                                these_animals, rep(NA, 7-length(these_animals)))
  counter2 <- counter2 + 1
}
study_focus1$Study_ID <- as.numeric(study_focus1$Study_ID)
study_focus2$Study_ID <- as.numeric(study_focus2$Study_ID)
study_focus1 <- left_join(study_focus1, animal_families, join_by("Species" == "Animal_Species_Corrected"))

unique_studies <- frug_fol_mut %>% 
  select(Study_ID, Study_Site, Plant_Species_Corrected, Habitat) %>% 
  left_join(., unique_locs, by = "Study_Site")  %>% # getting lat and long
  filter(Study_ID %in% study_list) %>% 
  group_by(Study_ID, Study_Site, Habitat, Country, Region) %>% 
  summarize(n_obs = n(), Lat = mean(Lat, na.rm = TRUE), Long = mean(Long, na.rm = TRUE)) %>% 
  as.data.frame(.)

# Save two versions: one with one row per study, two with one per per study-animal
study_info2 <- left_join(unique_studies, study_focus2, by = "Study_ID")
study_info1 <- left_join(study_focus1, unique_studies, by = "Study_ID")


# ------------- PART C: USE COMBINED FRUGIVORY AND FOLIVORY TO MAKE OCCURRENCE DATA ---------------------#

## Create occurence matrices using data FROM FOLIVORY AND FRUGIVORY
O_P <- matrix(0, nrow = nP, ncol = nS) 
rownames(O_P) <- plant_list
O_A <- matrix(0, nrow = nA, ncol = nS) 
rownames(O_A) <- animal_list
colnames(O_P) <- colnames(O_A) <- study_names

## Extract occurrence information using a variety of methods, determining occurrence probability when
# Species occurs in the same study
# Species occurs at the same site, but a different study
# Species occurs in the same country, but a different site

for(ss in 1:nS){
  # Subset all studies occuring at the same site
  study_id <- study_list[ss]
  study_site <- unique(frug_fol_mut[frug_fol_mut$Study_ID == study_id, 12]) 
  study_region <- unique(frug_fol_mut[frug_fol_mut$Study_ID == study_id, 23])  
  study_country <- unique(frug_fol_mut[frug_fol_mut$Study_ID == study_id, 21 ])
  study_habitat <- unique(frug_fol_mut[frug_fol_mut$Study_ID == study_id, 22 ])
  study_dta <- frug_fol_mut[frug_fol_mut$Study_ID == study_id,] # same study
  site_dta <- frug_fol_mut[frug_fol_mut$Study_Site == study_site,] # same site, different study
  reg_hab_dta <- frug_fol_mut[frug_fol_mut$Region == study_region & frug_fol_mut$Habitat == study_habitat, ] # same region and habitat
  country_hab_dta <- frug_fol_mut[frug_fol_mut$Country == study_country & frug_fol_mut$Habitat == study_habitat, ] # same country and habitat
  hab_dta <- frug_fol_mut[frug_fol_mut$Habitat == study_habitat, ] # same habitat
  country_dta <- frug_fol_mut[frug_fol_mut$Country == study_country, ] # same country
  reg_dta <- frug_fol_mut[frug_fol_mut$Region == study_region, ] # same region
  for(p in 1:nP){
    this.plant <- plant_list[p]
    if(length(study_dta$Plant_Species_Corrected[study_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 1
    } else if(length(site_dta$Plant_Species_Corrected[site_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 0.75
    } else if(length(country_hab_dta$Plant_Species_Corrected[country_hab_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 0.5
    } else if(length(reg_hab_dta$Plant_Species_Corrected[reg_hab_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 0.45
    } else if(length(hab_dta$Plant_Species_Corrected[hab_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 0.25
    }  else if(length(country_dta$Plant_Species_Corrected[country_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 0.1
    } else if(length(reg_dta$Plant_Species_Corrected[reg_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 0.05
    } 
  }  
  for(a in 1:nA){
    this.animal <- animal_list[a]
    if(length(study_dta$Animal_Species_Corrected[study_dta$Animal_Species_Corrected==this.animal])>0){
      O_A[a,ss] <- 1
    } else if(length(site_dta$Animal_Species_Corrected[site_dta$Animal_Species_Corrected==this.animal])>0){
      O_A[a,ss] <- 0.75
    } else if(length(country_hab_dta$Animal_Species_Corrected[country_hab_dta$Animal_Species_Corrected==this.animal])>0){
      O_A[a,ss] <- 0.5
    } else if(length(hab_dta$Animal_Species_Corrected[hab_dta$Animal_Species_Corrected==this.animal])>0){
      O_A[a,ss] <- 0.25
    }
  }
}
  
# Isolated studies: most plants do not have overlap except at the country
isolated_studies <- O_P[, colnames(O_P) %in% c("Study.44", "Study.19", "Study.8", "Study.52", "Study.33")]
table(c(isolated_studies))

### EDA with occurrences 

## Distribution of occurrences for plants and animals
O_P_f <- factor(O_P, labels = c("Recorded in region",  "Recorded in country", 
                                "Recorded in habitat but not country", 
                                "Recorded in region and habitat", "Recorded in country and habitat",
                                "Recorded at site", "Recorded in study"))
O_A_f <- factor(O_A, labels = c("Not recorded in habitat",  "Recorded in habitat but not country", "Recorded in country and habitat", "Recorded at site", "Recorded in study"))

barplot(prop.table(table(O_P_f)),
        ylab = "Proportion",
        xlab = "Occurrence Status")

barplot(prop.table(table(O_A_f)),
        ylab = "Proportion",
        xlab = "Occurrence Status")




#---------------------- PART C: processing the adjacency matrices --------------------------#
# Convert raw data into a list of adjacency matrices
#Adj_matrix <- array( 0, dim= c(nrow(animal_list), nrow(plant_list), nrow(study_list)))
#interactions <- matrix(0, nrow= nrow(animal_list), ncol= nrow(plant_list))
#Adj_matrix_new <- lapply(seq_len(length(studies)), function(X) interactions)

# this should be redundant
study_list <- frugivore_mut %>%
  select(Study_ID)%>%
  unique()

frugivore_mut_new <- frugivore_mut %>% select(Animal_Species_Corrected, Plant_Species_Corrected, Study_ID)
frugivore_mut_new$Interaction <- 1

Adj_matrix_new <- list()
for(study in 1:nS){
    Adj_df <- frugivore_mut_new %>% filter(Study_ID==study_list[study,]) %>% 
      select(Animal_Species_Corrected, Plant_Species_Corrected, Interaction) %>% 
      distinct() %>%pivot_wider(names_from = Plant_Species_Corrected, values_from=Interaction)
    Adj_df  <- as.data.frame(Adj_df )
    rownames(Adj_df) <- Adj_df$Animal_Species_Corrected
    Adj_df <- Adj_df[,-1]
    Adj_df[is.na(Adj_df)] <- 0
    Adj_matrix <- as.matrix(Adj_df)
    #Adj_matrix_new[[study]]  <- Adj_matrix
    s.id <- paste0("Study.",study_list[study,])
    Adj_matrix_new[[s.id]] <- Adj_matrix
}


if(save_files){
  saveRDS(Adj_matrix_new, file = paste0(save_path, "Adj_matrix_new.rds"))
  save(O_P, file = paste0(save_path, 'OP_full.dat'))
  save(O_A, file = paste0(save_path, 'OA_full.dat'))
  save(study_info1, file = paste0(save_path, 'study_info1.dat'))
  save(study_info2, file = paste0(save_path, 'study_info2.dat'))
  save(unique_locs, file = paste0(save_path, 'loc_info.dat'))
}



