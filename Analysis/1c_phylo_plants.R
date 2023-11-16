# Builds plant phylogenies and also adds additional taxonomical info to traits files
# Uses taxonomical info to impute missing fruit lengths

# -------- TO DO --------- #

# The directory where the original data are:
data_path <- 'RawData/'
# The directory where the processed data should be saved:
save_path <- 'ProcessedData/'
# Whether the processed data should be saved or not:
save_files <- TRUE


# --------- BEGINNING -------- #

# Setting the working directory.
#setwd(wd_path)
if(!("V.PhyloMaker" %in% rownames(installed.packages()))){
  install.packages("remotes")
  remotes::install_github("jinyizju/V.PhyloMaker")
}

# Loading libraries.
library(ape)
library(V.PhyloMaker)
library(gplots)
library(magrittr)
library(dplyr)

# Load plant traits to get plant names
load(paste0(save_path, "traits_p_709.dat"))
load(paste0(save_path, "Obs_W_partiallycleaned.dat"))

uni_plants <- cbind(species = traits.p.709$Species, # this has the same order as the adjacency matrices
                    genus = traits.p.709$Genus, 
                    family = traits.p.709$Family)

uni_plants <- as.data.frame(unique(uni_plants))
nP <- nrow(uni_plants)

#------------------- GET PHYLOGENY ----------------------#

y <- phylo.maker(uni_plants, scenarios = "S3") # default scenario

#create a random tree to test later
y_rand <- rtree(709)

# Using the phylogenetic tree in the ape R package to get a phylogenetic
# correlation matrix.
Cv_phylo <- round(ape::vcv(y[[1]], corr = TRUE), 10)
Cv_phylo1 <- round(ape::vcv(y_rand, corr = TRUE), 10)

# Clean taxonomy 
traits.plants <- traits.p.709
traits.plants$Family[traits.plants$Family=="Maranthaceae" ] <- "Marantaceae"
traits.plants$Family[traits.plants$Family=="Thromandersiaceae" ] <- "Thomandersiaceae"
traits.plants$Family[traits.plants$Family=="Ebanaceae" ] <- "Ebenaceae"
traits.plants$Order <- traits.plants$Family
traits.plants$Order[traits.plants$Order=="Acanthaceae"|traits.plants$Order=="Bignoniaceae"| traits.plants$Order=="Lamiaceae"|traits.plants$Order=="Oleaceae"|traits.plants$Order=="Thomandersiaceae"|traits.plants$Order=="Verbenaceae"] <- "Lamiales"
traits.plants$Order[traits.plants$Order=="Achariaceae"|traits.plants$Order== "Calophyllaceae" |traits.plants$Order=="Chrysobalanaceae"|traits.plants$Order=="Clusiaceae"| traits.plants$Order=="Dichapetalaceae"|traits.plants$Order=="Erythroxylaceae"|traits.plants$Order== "Euphorbiaceae"|traits.plants$Order=="Humiriaceae"|traits.plants$Order=="Hypericaceae"|traits.plants$Order=="Irvingiaceae"|traits.plants$Order=="Linaceae"|traits.plants$Order=="Pandaceae"|traits.plants$Order=="Passifloraceae"|traits.plants$Order=="Phyllanthaceae"|traits.plants$Order=="Putranjivaceae"|traits.plants$Order=="Rhizophoraceae"|traits.plants$Order=="Salicaceae"|traits.plants$Order=="Turneraceae"] <- "Malpighiales"
traits.plants$Order[traits.plants$Order=="Anacardiaceae"|traits.plants$Order=="Burseraceae"|traits.plants$Order=="Meliaceae"|traits.plants$Order=="Rutaceae"|traits.plants$Order=="Sapindaceae"|traits.plants$Order=="Simaroubaceae"] <- "Sapindales"
traits.plants$Order[traits.plants$Order=="Anisophylleaceae"| traits.plants$Order== "Cucurbitaceae"] <- "Cucurbitales"
traits.plants$Order[traits.plants$Order=="Annonaceae"|traits.plants$Order=="Canellaceae"|traits.plants$Order=="Myristicaceae"] <- "Magnoliales"
traits.plants$Order[traits.plants$Order=="Apocynaceae" | traits.plants$Order== "Gentianaceae"|traits.plants$Order=="Loganiaceae"|traits.plants$Order=="Rubiaceae"] <- "Gentianales"
traits.plants$Order[traits.plants$Order=="Aquifoliaceae"] <- "Aquifoliales"
traits.plants$Order[traits.plants$Order=="Araliaceae"] <- "Apiales"
traits.plants$Order[traits.plants$Order=="Arecaceae"] <- "Arecales"
traits.plants$Order[traits.plants$Order=="Asparagaceae"] <- "Asparagales"
traits.plants$Order[traits.plants$Order=="Asteraceae"| traits.plants$Order== "Compositae"] <- "Asterales"
traits.plants$Order[traits.plants$Order=="Balsaminaceae"| traits.plants$Order=="Ebenaceae"|traits.plants$Order== "Lecythidaceae"|traits.plants$Order== "Primulaceae"|traits.plants$Order=="Sapotaceae"|traits.plants$Order=="Theaceae"] <- "Ericales"
traits.plants$Order[traits.plants$Order=="Boraginaceae"] <- "Boraginales"
traits.plants$Order[traits.plants$Order=="Buxaceae"] <- "Buxales"
traits.plants$Order[traits.plants$Order=="Cannabaceae"|traits.plants$Order=="Montiniaceae"|traits.plants$Order=="Moraceae"|traits.plants$Order=="Rhamnaceae"|traits.plants$Order=="Rosaceae"|traits.plants$Order=="Ulmaceae"|traits.plants$Order=="Urticaceae"] <- "Rosales"
traits.plants$Order[traits.plants$Order=="Capparaceae" |traits.plants$Order== "Caricaceae"] <- "Brassicales"
traits.plants$Order[traits.plants$Order=="Casuarinaceae" ] <- "Fagales"
traits.plants$Order[traits.plants$Order=="Celastraceae"| traits.plants$Order== "Hippocrateaceae"] <- "Celastrales"
traits.plants$Order[traits.plants$Order=="Combretaceae"|traits.plants$Order=="Melastomataceae"|traits.plants$Order=="Myrtaceae"|traits.plants$Order=="Penaeaceae"] <- "Myrtales"
traits.plants$Order[traits.plants$Order=="Commelinaceae" ] <- "Commelinales"
traits.plants$Order[traits.plants$Order=="Connaraceae" ] <- "Oxalidales"
traits.plants$Order[traits.plants$Order=="Convolvulaceae"|traits.plants$Order=="Solanaceae"] <- "Solanales"
traits.plants$Order[traits.plants$Order=="Costaceae"|traits.plants$Order=="Marantaceae"|traits.plants$Order=="Musaceae"|traits.plants$Order=="Zingiberaceae"] <- "Zingiberales"
traits.plants$Order[traits.plants$Order=="Cupressaceae"|traits.plants$Order=="Podocarpaceae"] <- "Pinales"
traits.plants$Order[traits.plants$Order=="Cyperaceae"|traits.plants$Order=="Poaceae"] <- "Poales"
traits.plants$Order[traits.plants$Order=="Dioscoreaceae" ] <- "Dioscoreales"
traits.plants$Order[traits.plants$Order=="Fabaceae"|traits.plants$Order=="Polygalaceae"] <- "Fabales"
traits.plants$Order[traits.plants$Order=="Flacourtiaceae" ] <- "Bixales"
traits.plants$Order[traits.plants$Order=="Hernandiaceae"| traits.plants$Order=="Lauraceae"|traits.plants$Order=="Monimiaceae"] <- "Laurales"
traits.plants$Order[traits.plants$Order=="Icacinaceae" ] <- "Icacinales"
traits.plants$Order[traits.plants$Order=="Loranthaceae"|traits.plants$Order=="Olacaceae"|traits.plants$Order=="Opiliaceae"|traits.plants$Order=="Santalaceae"] <- "Santalales"
traits.plants$Order[traits.plants$Order=="Malvaceae" ] <- "Malvales"
traits.plants$Order[traits.plants$Order=="Menispermaceae" ] <- "Ranunculales"
traits.plants$Order[traits.plants$Order=="Phytolaccaceae" ] <- "Caryophyllales"
traits.plants$Order[traits.plants$Order=="Piperaceae" ] <- "Piperales"
traits.plants$Order[traits.plants$Order=="Pittosporaceae" ] <- "Apiales"
traits.plants$Order[traits.plants$Order=="Smilacaceae" ] <- "Liliales"
traits.plants$Order[traits.plants$Order=="Vitaceae" ] <- "Vitales"
traits.plants$Order[traits.plants$Order=="Zygophyllaceae" ] <- "Zygophyllales"

sort(unique(traits.plants$Order))

# If we want to use taxonomic correlation instead (not currently saved)
Cv_tax <- diag(nP)
p.names <- uni_plants$species
dimnames(Cv_tax) <- list(p.names, p.names)

for (i1 in 1 : (nP - 1)) {
  for (i2 in (i1 + 1) : nP) {
    
    if (traits.plants$Genus[i1] == traits.plants$Genus[i2]) {
      Cv_tax[i1, i2] <- Cv_tax[i2, i1] <- 0.75
    } else if (traits.plants$Family[i1] == traits.plants$Family[i2]) {
      Cv_tax[i1, i2] <- Cv_tax[i2, i1] <- 0.5
    } else if (traits.plants$Order[i1] == traits.plants$Order[i2]) {
      Cv_tax[i1, i2] <- Cv_tax[i2, i1] <- 0.25
    }
    
  }
}

# Re-ordering the correlation matrix in the way that is in our data.
use_order <- sapply(as.character(uni_plants$species),
                    function(r) which(rownames(Cv_phylo) == r))
Cv_phylo <- Cv_phylo[use_order, use_order]


# Most plants are not very closely related
gplots::heatmap.2(Cv_phylo, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE,
                  labRow = FALSE, labCol = FALSE)
hist(c(Cv_phylo))


# Use new data on Order to impute remaining missing values for fruit length
Obs_W_cleaned <- Obs_W %>%
                 mutate(Species = row.names(.)) %>%
                 left_join(., traits.plants[, c(1,2,3,7)], by = "Species") %>%
                 group_by(Order) %>%
                 mutate(FruitLength.order = mean(FruitLength, na.rm = TRUE)) %>%
                 mutate(FruitLength = ifelse(is.na(FruitLength), FruitLength.order, FruitLength)) %>%
                 ungroup()
Obs_W_cleaned[is.na(Obs_W_cleaned$FruitLength), "Species"]
Obs_W_cleaned <- Obs_W_cleaned %>% # Filling in the rest by Google search and eyeballing
                  mutate(FruitLength = ifelse(Species == "Buxus_acutata", 10, FruitLength)) %>%
                  mutate(FruitLength = ifelse(Genus == "Palisota", 5, FruitLength)) %>%
                  mutate(FruitLength = ifelse(Species == "Rourea_orientalis", 11 , FruitLength)) %>% # Camille paper
                  mutate(FruitLength = ifelse(Species == "Dioscorea_sansibarensis", 15 , FruitLength)) %>%
                  mutate(FruitLength = ifelse(Species == "Aneilema_nyasense", 12 , FruitLength)) %>% #https://plants.jstor.org/compilation/aneilema.nyasense
                  as.data.frame(.) %>%
                  set_rownames(.$Species) %>%
                  select(FruitLength, meanWD, IUCN_Endangered) 

Obs_W <- Obs_W_cleaned                  


if (save_files) {
  save(Cv_phylo, file = paste0(save_path, 'Cv_phylo.dat'))
  save(Obs_W, file = paste0(save_path, 'Obs_W.dat'))
  save(traits.plants, file = paste0(save_path, 'traits_p_709_clean.dat'))
}



