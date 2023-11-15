# EDA

# The directory where the analysis is performed:
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V2"
# Where the processed data are saved:
data_path <- 'ProcessedData/'
# Where you want to save MCMC results:
result_path <- 'Results/'
save_path <- 'EDA/'
# Where the functions are available:
source_path <- 'HelperScriptsNew/'

setwd(wd_path)
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

load(paste0(data_path, "study_info1.dat"))
load(paste0(data_path, "study_info2.dat"))
load(paste0(data_path, "loc_info.dat")) #unique_locs
load(paste0(data_path, 'obs_A_mammals.dat'))
load(paste0(data_path, 'obs_OP.dat')) # site level obs_OP

# ---------------------- SUBSET study info to exclude bird studies ----------------#

study_list <- unlist(dimnames(A.obs.m)[3])
study_ids <- as.numeric(gsub("Study.", "", study_list))
mammals <- rownames(A.obs.m)

study_info1 <- study_info1 %>% 
               filter(Study_ID %in% study_ids, Species %in% mammals) %>% 
               mutate(Species = gsub("_", " ", Species))
study_info_genus <- study_info1 %>%
               mutate(Genus = gsub(" .*", "", Species)) %>% 
               unique()


# ------------------------ MAPS ----------------------------------#
# Getting the African map:
map_af <- get_stamenmap(bbox = c(left = min(unique_locs$Long) - 2,
                                 bottom = min(unique_locs$Lat) - 2,
                                 right = max(unique_locs$Long) + 2,
                                 top = max(unique_locs$Lat) + 2),
                        zoom = 6,
                        maptype = 'terrain',
                        color = "color", force = TRUE)


# Plotting the locations of the bird and plant interactions to
# visualize geographic bias:
ggmap(map_af, extent = "device", legend = "right") +
  geom_point(data = unique_locs, aes(x = Long, y = Lat, size = n_studies), pch = 21, color = "black",
             fill = '#E86861', alpha = 0.75) + 
  scale_size_continuous(range = c(5,10), limits = c(1, 11), breaks = c(1,5,11)) + 
  labs(size = "Number of studies", title = "Frugivory study sites") + 
  theme(legend.position = "bottom")

# slides version
png(filename = paste0(save_path, "slides_sitesmap.png"), width = 1000, height = 500)
ggmap(map_af, extent = "device", legend = "right") +
  geom_point(data = unique_locs, aes(x = Long, y = Lat, size = n_studies), pch = 21, color = "black",
             fill = '#E86861', alpha = 0.75) + 
  scale_size_continuous(range = c(5,10), limits = c(1, 11), breaks = c(1,5,11)) + 
  labs(size = "Number of studies") + 
  theme(legend.position = "bottom", text = element_text(size = 25))
dev.off()

# Plotting the locations and animal familes
ggmap(map_af, extent = "device", legend = "right") +
  geom_jitter(data = study_info1, aes(x = Long, y = Lat, fill = Animal_Family), size = 5, pch = 21,
             color = "black",alpha = 0.5, width = 1.5, height = 1.5) + 
  #scale_size_continuous(range = c(5,10), limits = c(1, 11), breaks = c(1,5,11)) + 
  labs(fill = "Focus animal family", title = "Frugivory study sites") + 
  theme(legend.position = "bottom", text = element_text(family = "serif", size = 20))

# Plotting the locations and animal familes: slides version + highlight baboons
png(filename = paste0(save_path, "slides_baboonmap.png"), width = 1000, height = 500)
ggmap(map_af, extent = "device", legend = "right") +
  geom_jitter(data = study_info1, aes(x = Long, y = Lat, fill = Animal_Family), size = 5, pch = 21,
              color = "black",alpha = 0.5, width = 1.5, height = 1.5) + 
  #scale_size_continuous(range = c(5,10), limits = c(1, 11), breaks = c(1,5,11)) + 
  labs(fill = "Focus animal family", title = "") + 
  theme(legend.position = "bottom", text = element_text(size = 25)) + 
  geom_label(data = subset(study_info1, Species %in% c("Papio anubis", "Papio cynocephalus")), 
             aes(x = Long, y = Lat, label = Species)) 
dev.off()



# Plotting the locations and animal species
ggmap(map_af, extent = "device", legend = "right") +
  geom_jitter(data = study_info_genus, aes(x = Long, y = Lat, fill = Genus), size = 5, pch = 21,
             color = "black",alpha = 0.75, width = 1, height = 1) + 
  #scale_size_continuous(range = c(5,10), limits = c(1, 11), breaks = c(1,5,11)) + 
  labs(fill = "Genus of Focus Animal", title = "Frugivory study sites") + 
  theme(legend.position = "bottom")

# ----------------------- ISOLATED/UNIQUE SITES ------------------------#

## Investigate isolated sites

# Idea: for each study, look at the proportion of the plants in the study are found only at that site
# proportion of study plants
# And how many studies do we have for that same site?
# I would have to go back to the raw data to get that, good to do, save as a separate file

## Compute # plants observed in the study/ # plants observed at that site
# E.g. if ratio = 0.9, this means that 90% of plants that were observed at that site were observed being eaten in this study
# A very high ratio suggests poor study coverage of the relevant ecosystem and possible selection bias
obs_to_site = data.frame(study = study_list, 
                         obs_to_site = apply(obs_OP, 2, function(x) sum(x == 1))/apply(obs_OP, 2, function(x) sum(x > 0.5)))

## Compute # plants observed in the study/ # of plants observed in the same habitat
# E.g. if ratio = 0.9 this means that 90% of plants that were observed in that habitat were observed being eaten in this study
obs_to_hab = data.frame(study = study_list, 
                        obs_to_hab = apply(obs_OP, 2, function(x) sum(x == 1))/apply(obs_OP, 2, function(x) sum(x > 0.1)))

study.plants <- data.frame(study = study_list, 
  obs_to_habcountry = apply(obs_OP, 2, function(x) sum(x == 1))/apply(obs_OP, 2, function(x) sum(x >0.45)))%>% 
  arrange(-obs_to_habcountry) %>% 
  left_join(., obs_to_site) %>% 
  left_join(., obs_to_hab) %>% 
  pivot_longer(!study, names_to = "Type", values_to = "Ratio") %>%
  mutate(row = row_number())  

## Another way of looking at it
# Number of studies at the same site (already done with unique locs)
# Number of studies in the same habitat and region
hab_country_overlap <- study_info2 %>% select(Study_ID, Study_Site, Habitat, Country, Region) %>% 
                unique() %>% 
                group_by(Habitat, Country) %>% 
                summarize(n_studies_in_hab_country = n()) %>% 
                left_join(study_info2, .) 


## Another way of looking at it
# Number of studies at the same site (already done with unique locs)
# Number of studies in the same habitat and region
hab_region_overlap <- study_info2 %>% select(Study_ID, Study_Site, Habitat, Country, Region) %>% 
  unique() %>% 
  group_by(Habitat, Region) %>% 
  summarize(n_studies_in_hab_region = n()) %>% 
  left_join(study_info2, .)

# Overlap in habitat and region
ggmap(map_af, extent = "device", legend = "right") +
  geom_jitter(data = hab_region_overlap, aes(x = Long, y = Lat, size = n_studies_in_hab_region), pch = 21,
              fill = '#E86861', color = "black",alpha = 0.75, width = 1, height = 1) + 
  scale_size_continuous(range = c(0.5,10), limits = c(1, 30), breaks = c(1,10,20, 30)) + 
  labs(size = "Number of studies in same habitat and region", title = "Frugivory study sites") + 
  theme(legend.position = "bottom")


# Overlap in habitat and country
ggmap(map_af, extent = "device", legend = "right") +
  geom_jitter(data = hab_country_overlap, aes(x = Long, y = Lat, size = n_studies_in_hab_country), pch = 21,
              fill = '#E86861', color = "black",alpha = 0.75, width = 1, height = 1) + 
  scale_size_continuous(range = c(1,6), limits = c(1, 16), breaks = c(1,5,10,15)) + 
  labs(size = "Number of studies in same habitat and country", title = "Frugivory study sites") + 
  theme(legend.position = "bottom")



# Just visualizing the regions
ggmap(map_af, extent = "device", legend = "right") +
  geom_jitter(data = hab_country_overlap, aes(x = Long, y = Lat, color = Country, shape = Region),
              size = 4, alpha = 0.75, width = 1, height = 1) + 
  #scale_size_continuous(range = c(1,6), limits = c(1, 16), breaks = c(1,5,10,15)) + 
  labs( title = "Frugivory study sites by country and region") + 
  theme(legend.position = "bottom")


ggplot(study.plants, aes(x = reorder(study, row), y = Ratio, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  labs(title = "Proportion unique plants") + 
  xlab("") + 
  scale_fill_discrete(labels = c("Habitat level", "Habitat-country level", "Site level"), name = NULL) +
  theme_minimal() + 
  theme(text = element_text(family = "serif", size = 20), axis.text.x = element_text(angle = 90, size = 15)) 


# Plotting the locations and animal species
ggmap(map_af, extent = "device", legend = "right") +
  geom_jitter(data = filter(study.plants, Type == "obs_to_hab"), aes(x = Long, y = Lat, fill = Genus), size = 5, pch = 21,
              color = "black",alpha = 0.75, width = 1, height = 1) + 
  scale_size_continuous(range = c(5,10), limits = c(1, 11), breaks = c(1,5,11)) + 
  labs(fill = "Genus of Focus Animal", title = "Frugivory study sites") + 
  theme(legend.position = "bottom")


# Other things to do
#Plot plant uniqueness on map
# Create a table by species of number of studies, number of sites (already partly did i think?)

ggplot(data = unique_locs, mapping = aes(x = Study_Site, y = n_studies)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90), text = element_text(family = "serif", size = 20)) 

# Number of studies by family
study_info1 %>% group_by(Animal_Family) %>% 
  summarize(Effort = n()) %>%
  ggplot() + geom_col(aes(x = fct_reorder(Animal_Family, -Effort), y = Effort)) + 
  xlab("Family") + 
  theme_minimal() +
  theme(text = element_text(family = "serif", size = 20), axis.text.x = element_text(angle = 90, size = 15)) 

# Number of studies by family
study_info1 %>% group_by(Species) %>% 
  summarize(Effort = n()) %>%
  ggplot() + geom_col(aes(x = fct_reorder(Species, -Effort), y = Effort)) + 
  xlab("Species") + 
  theme_minimal() +
  theme(text = element_text(family = "serif", size = 20), axis.text.x = element_text(angle = 90, size = 15)) 

  
  traits.verts %>% group_by(Family) %>%
  summarize(Effort = sum(effort)) %>%
  ggplot() + geom_col(aes(x = fct_reorder(Family,-Effort), y = Effort)) + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90)) + 
  xlab("Family") + 
  theme(text = element_text(family = "serif", size = 20), axis.text.x = element_text(angle = 90, size = 15)) 

  ## How many studies per site? This informs how much undersampling we think there might be
# REDO this using saved data unique_locs

ggplot(unique_locs) + geom_bar(aes(x = n_studies)) +
    scale_x_continuous(breaks=0:12) + 
    theme_minimal() + 
    xlab("Number of Studies") + 
    ylab("Number of Sites") + 
    ggtitle("How many times is each site studied?")
  
  # Assessment: the vast majority of sites only have a single study associated
  # so it is very likely that there are plant species present which were not recorded 
  # at either the study level or the site-level. 
  # this motivates increasing the probability for P(present | same habitat and country)
  
  
