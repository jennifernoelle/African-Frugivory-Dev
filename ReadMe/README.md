---
editor_options: 
  markdown: 
    wrap: 72
---

# African Frugivory interactions

This repo contains extensions to the code and model developed in
Papadogeorgou (2023) to account for the effects of extreme taxonomic
bias. In this project, we aim to impute missing mammal-plant
interactions from a multi-study data set.

## Data set

The key raw data sources are:

-   Frugivory only data: "Frugivory.csv"

-   Frugivory and folivory meta-anlaysis: "Frugivory_folivory.csv"

-   Site metadata: "Site_metadata.csv"

-   Study locations: "Study_locations.xlsx"

Additional information on the construction of the meta-analysis is
included in this folder.

Plant phylogenies are acquired using the V.PhyloMaker R package, and the
code to do so is available in the Analysis/ folder under name
1c_phylo_plants.R. A hylogenetic correlation matrix for vertebrates is
computed via the the package ape using a consensus tree obtained from
VertLife and provided in this repo in the file mammals2.nex; the code to
perform this analysis is in the Analysis/ folder under the name
1b_phylo_verts.R.

## Code

You must create the folder ProcessedData under the main directory for
the project

The folder HelperScriptsJK/ includes functions (modified from
<https://github.com/gpapadog/BiasedNetwork>) that are used in the
analysis code.

The code for the analysis is in the folder Analysis/. The numbers in the
beginning of the file names represent the order with which the files
should be used/run.

-   0a_data_prep.R

-   0b_more_cleaning.R

-   1a_subset_data.R: subset data to include only mammals and plants

-   1b_phylo_verts.R: generate vertebrate phylogenetic correlation
    matrix

-   1c_phylo_plants.R: generate plant phylogenetic correlation matrix

-   2a_subsetting.R: further subset the data to exclude isolated species

-   2b_analysis_defaultprior.R: perform the analysis using default prior
    occurrence probabilities

-   2b_analysis_experprior.R: perform the analysis using expert defined
    prior occurrence probabilities

-   0_Bias_visualization_cleaning.R: This code MUST be run before any
    subsequent analysis. Basic cleaning to ensure consistency among
    naming conventions between adjacency matrices, Focus data, traits
    data, and to remove duplicates. Transforms a few variables in the
    traits data. Plots taxonomic bias of the study.

-   1a_Mammals_subset_data.R: This code MUST be run before any
    subsequent analysis. This code loads the data and processes them
    into matrices and data frames that can be used in subsequent
    analysis, turning the partial observation matrices into full
    observation matrices by adding zeros appropriately. Then the file
    subsets out mammals and the plants relevant to mammals (i.e., any
    plant that any mammal is observed to eat somewhere in the data). The
    corresponding Obs, trait, F files are created. The processed data
    need to be saved to a local directory. That directory can be
    specified at the beginning of the code.

-   1b_phylo_verts.R This code MUST be run. It uses downloaded
    phylogenetic trees to acquire an estimate of the phylogenetic
    correlation matrix for the vertebrate species, then subsets mammals.

-   1b_phylo_plants.R This code MUST be run. It uses an existing R
    package to acquire an estimate of the phylogenetic correlation
    matrix for the plant species.

-   2a_analysis.R: The main code for performing the analysis using the
    proposed method. We ran four MCMC chains in parallel. You will need
    to create a folder where analysis results should be saved.

-   3a_cross_validation.R: Performing cross validation by helding out
    recorded interactions. Code for cross-validation based on our
    method.

-   4_trait_matching.R: Performing the trait matching algorithm of the
    manuscript

-   5_plot_results.R: Code that plots all the results that are shown in
    the manuscript

### Note

For replicating the analysis results in this code, you will need to
specify a directory where results from the individual studies can be
saved. We recommend you create a folder at the same level as the
Analysis/, Data/, and HelperScripts/ folders that is named Results/.

## References

Bello, C., Galetti, M., Montan, D., Pizo, M. A., Mariguela, T. C.,
Culot, L., Bufalo, F., Labecca, F., Pedrosa, F., Constantini, R., Emer,
C., Silva, W. R., da Silva, F. R., Ovaskainen, O., & Jordano, P. (2017).
Atlantic frugivory: A plant-frugivore interaction data set for the
Atlantic Forest. Ecology, 98(6), 1729.
<https://doi.org/10.1002/ecy.1818>

## Acknowledgements

This project has received funding from the European Research Council
(ERC) under the European Union's Horizon 2020 research and innovation
programme (grant agreement No 856506; ERC-synergy project LIFEPLAN)
