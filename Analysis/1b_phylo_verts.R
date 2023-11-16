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

# Loading libraries.
library(ape)
library(data.table)
library(superheat)

load(paste0(save_path, "Obs_X.dat"))
#load(paste0(save_path, "Obs_W.dat"))
m.names <- rownames(Obs_X)
#p.names <- rownames(Obs_W)

# --------------------
# PART 1: Loading the (single) vertebrate phylogenetic trees from the donwloaded file.
# --------------------

x <- read.nexus(paste0(data_path, 'mammals2.nex'))

consensus <- consensus(x)

plot(consensus)

# There are 100 phylogenetic trees. Getting the name of the species included
# in those trees:
x_names <- sapply(x, function(y) y$tip.label) # each column is a tree, vert names in rows
table(apply(x_names, 1, function(y) length(unique(y)))) # how many different verts names at each leaf (across trees)
length(unique(as.character(x_names)))
# From these results we see that the 100 phylogenetic trees do not have the
# species in the same order, even though they include the same species overall.


# --------------------
# PART 2: Checking the overlap of species in our data and in the phylogenies.
# --------------------

# Checking how many ave species from our interaction data exist in the
# phylogenetic trees.
sort(unique(as.character(x_names)))
sort(m.names[which(!(m.names %in% x_names))])

# Which ave species are missing:
miss_m <- m.names[which(!(m.names %in% x_names[, 1]))]

# Some of these species are also known under alternative names in ecology.
# Specify these alternative names:

alt_names <- rep(NA, length(miss_m))
alt_names[2] <- 'Philantomba_monticola'
alt_names[1] <- 'Cercopithecus_lhoesti'
 #alt_names[4] <- 'Loxodonta_africana' # we actually do want to keep this na, but will add info back after processing
alt_names[3] <- 'Procolobus_badius'

cbind(miss_m, alt_names)

alt_names[!is.na(alt_names)] %in% x_names[, 1]  # The new names exist in phylogeny
alt_names[!is.na(alt_names)] %in% m.names  # Check for overlap with other names

# Loxodonta cyclotis is that has only very recently been recognized as a separate species 
# form Loxodonta africana, so it's very genetically similar Since Loxodonta africana already
# exists in the data, we're going to just make our correlation matrix without Loxondonta
# cyclotis and then copy that row/column into the appropriate location



# --------------------
# PART 3: Getting a phylogenetic correlation matrix for each phylogenetic tree
# --------------------

# Use all 1000 samples from the phylogenetic tree
use_post_samples <- 1 : length(x)

num_species <- length(m.names) - sum(is.na(alt_names))
s <- array(NA, dim = c(num_species, num_species, length(use_post_samples)))
dim(s)

# For each phylogenetic tree from the 1000 posterior samples, calculate the
# phylogenetic correlation matrix using the ape R package.
t1 <- Sys.time()
for (ii in 1 : length(use_post_samples)) {
  if (ii %% 10 == 0) print(ii)
  this_sample <- use_post_samples[ii]
  # Getting the correlation matrix:
  this_cov <- ape::vcv(x[[this_sample]], corr = TRUE)
  # Re-ordering to get species on the same order:
  this_cov <- this_cov[order(rownames(this_cov)), order(rownames(this_cov))]
  
  s[, , ii] <- this_cov
  if (ii == 1) {
    dimnames(s)[[1]] <- rownames(this_cov)
    dimnames(s)[[2]] <- colnames(this_cov)
  }
}
t2 <- Sys.time()
t2 - t1


# --------------------
# PART 5: Combining the correlation matrix in one phylogenetic matrix
# --------------------

# Take the mean correlation matrix across posterior samples and use this as the
# phylogenetic correlation matrix.
Cu_phylo <- apply(s, c(1, 2), mean)


# Use values from Loxodonta africana Loxodonta cyclotis
# And then reorder rows and columns to be consistent with original mammal ordering
w.L.africana <- which(row.names(Cu_phylo) == "Loxodonta_africana")
Cu.L.cyclotis <- Cu_phylo[, w.L.africana]

Cu_phylo <- rbind(Cu.L.cyclotis, Cu_phylo)
Cu_phylo <- cbind(c(1, Cu.L.cyclotis), Cu_phylo)

rownames(Cu_phylo)[1] <- colnames(Cu_phylo)[1] <- "Loxodonta_cyclotis"

#Create a random phylogeny to test later
v_rand <- rtree(29)
Cu_phylo2 <- round(ape::vcv(v_rand, corr = TRUE), 10)
colnames(Cu_phylo2) <- colnames(Cu_phylo)
rownames(Cu_phylo2) <- rownames(Cu_phylo)

gplots::heatmap.2(Cu_phylo, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE)


# For the alternative names we used, put in the original names:
for (ii in 1 : length(alt_names)) {
  if (!is.na(alt_names[ii])) {
    
    old_name <- miss_m[ii]
    wh_row <- which(rownames(Cu_phylo) == alt_names[ii])
    rownames(Cu_phylo)[wh_row] <- old_name
    colnames(Cu_phylo)[wh_row] <- old_name
    
  }
}

# Re-order the rows and columns to match the order of species in our data:
use_order <- sapply(m.names, function(r) which(rownames(Cu_phylo) == r))
use_order <- use_order[sapply(use_order, length) > 0]
use_order <- as.numeric(use_order)
Cu_phylo <- Cu_phylo[use_order, use_order]
gplots::heatmap.2(Cu_phylo, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE)

# Abbreviate names for plotting
plot.names <- gsub("_", " ", colnames(Cu_phylo))
g <- substr(gsub(" .*", "", plot.names), 1,3)
sp <- gsub(".* " , "", plot.names)
plot.names <- paste0(g, ". ", sp)
plot.cu <- Cu_phylo
colnames(plot.cu) <- rownames(plot.cu) <- plot.names

png(filename = "EDA/Phylo_mammals.png", width = 1000, height =1000)
superheat(X = plot.cu,
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 4,
          bottom.label.text.size = 4,
          bottom.label.size = 0.2, left.label.size = 0.12,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          heat.col.scheme = "grey", heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05), 
          title = "Phylogenetic correlation of mammals", 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 15)
dev.off()


if(save_files){
  save(Cu_phylo, file = paste0(save_path, 'Cu_phylo.dat'))
}

