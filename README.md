# African Frugivory interactions

This repo extends the model and sampler developed in Papadogeorgou (2023) to account for the effects of extreme taxonomic bias. In this project, we aim to impute missing mammal-plant interactions from a multi-study data set of mammalian frugivory in forested Sub-Saharan Africa.

## Data set

The key raw data sources are:

-   Frugivory only data: Frugivory.csv

-   Frugivory and folivory meta-anlaysis: Frugivory_folivory.csv

-   Site metadata: Site_metadata.csv

-   Study locations: Study_locations.xlsx

Additional information on the construction of the meta-analysis is included in the Background folder.

Plant phylogenies are acquired using the V.PhyloMaker R package, and the code to do so is available in the Analysis/ folder under name 1c_phylo_plants.R. A phylogenetic correlation matrix for vertebrates is computed via the the package ape using a consensus tree obtained from VertLife and provided in this repo in the file mammals2.nex; the code to perform this analysis is in the Analysis/ folder under the name 1b_phylo_verts.R.

## Code

The folder HelperScriptsNew/ includes functions (with substantial modifications from <https://github.com/gpapadog/BiasedNetwork>) that are used in the analysis code while the folder HelperScriptsOld/ includes more minor modifications for occurrence estimation within the original model for comparison purposes.

The code for the analysis is in the folder Analysis/. The numbers in the beginning of the file names represent the order with which the files should be used/run. In brief the content of each analysis file is as follows:

-   0a_network_cleaning.R: This code MUST be run before any subsequent analysis. Assembles key network and meta-data sources.

-   0b_species_cleaning.R: This code MUST be run before any subsequent analysis. Cleans species name and trait data.

-   1a_subset_data.R: This code MUST be run before any subsequent analysis. This file subsets data to include only mammals and plants.

-   1b_phylo_verts.R: This code MUST be run before any subsequent analysis. This file generates the vertebrate phylogenetic correlation matrix.

-   1c_phylo_plants.R: This code MUST be run before any subsequent analysis. This file generate the plant phylogenetic correlation matrix.

-   2a_subsetting.R: This code MUST be run before any subsequent analysis. This file further subsets the data to exclude isolated species.

-   2b_analysis_new_defaultprior.R: perform the analysis using default prior occurrence probabilities under the new sampler.

-   2c_analysis_new_expertprior.R: perform the analysis using expert defined prior occurrence probabilities under the new sampler

-   3a_cv_new_defaultprior.R: perform cross validation using default prior occurrence probabilities under the new sampler.

-   3b_cv_new_expertprior.R: perform cross validation using expert prior occurrence probabilities under the new sampler.

-   4a_plot_results_new.R: plot the posterior heatmap and cross validation performance metrics using the new sampler.

### Note

For replicating the analysis results in this code, you will need to specify a directory where processed data and results can be saved. We recommend you create a folder at the same level as the Analysis/, Data/, and HelperScripts/ folders that is named Results/ and ProcessedData/.
