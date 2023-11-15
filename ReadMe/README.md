# African Frugivory interactions

In this project, we aim to learn the underlying bird-plant interaction network from a multi-study
data set. Individual studies have inherent taxonomic or geographical biases due to focusing
on a specific set of species or geographical area. We correct for these biases by extending network
modeling to account for truly possible but unrecorded interactions. We also accommodate species
covariates to improve efficiency. Finally, we employ an approach to trait matching which allows us
to identify the traits that are most inflential in forming and detecting species interactions.

## Data set


Plant phylogenies are acquired using the V.PhyloMaker R package, and the code to do so is available
in the Analysis/ folder under name 1c_phylo_plants.R.

## Code
You must create the folder ProcessedData under the main directory for the project

The folder HelperScriptsPapadogeorgou/ includes functions that are used in the analysis code.

The code for the analysis is in the folder Analysis/. The numbers in the beginning of the file 
names represent the order with which the files should be used/ran. 
A short description of each file is included here. Each file is commented heavily.

- 0_Bias_visualization_cleaning.R: This code MUST be run before any subsequent analysis.
Basic cleaning to ensure consistency among naming conventions
between adjacency matrices, Focus data, traits data, and to remove duplicates. Transforms a few variables in 
the traits data. Plots taxonomic bias of the study. 

- 1a_Mammals_subset_data.R: This code MUST be run before any subsequent analysis. This code loads
the data and processes them into matrices and data frames that can be used in subsequent analysis, 
turning the partial observation matrices into full observation matrices by adding zeros appropriately. 
Then the file subsets out mammals and the plants relevant to mammals (i.e., any plant that any mammal
is observed to eat somewhere in the data). The corresponding Obs, trait, F files are created. 
The processed data need to be saved to a local directory. That directory can be specified at the beginning of the code.

- 1b_phylo_verts.R This code MUST be run. It uses downloaded phylogenetic trees to acquire an
estimate of the phylogenetic correlation matrix for the vertebrate species, then subsets 
mammals.

- 1b_phylo_plants.R This code MUST be run. It uses an existing R package to acquire an estimate
of the phylogenetic correlation matrix for the plant species.

- 2a_analysis.R: The main code for performing the analysis using the proposed method. We ran four
MCMC chains in parallel. You will need to create a folder where analysis results should be saved.

- 3a_cross_validation.R: Performing cross validation by helding out recorded interactions.
Code for cross-validation based on our method.

- 4_trait_matching.R: Performing the trait matching algorithm of the manuscript 

- 5_plot_results.R: Code that plots all the results that are shown in the manuscript

### Note

For replicating the analysis results in this code, you will need to specify a directory where
results from the individual studies can be saved. We recommend you create a folder at the same
level as the Analysis/, Data/, and HelperScripts/ folders that is named Results/.

## References

Bello, C., Galetti, M., Montan, D., Pizo, M. A., Mariguela, T. C., Culot, L.,
Bufalo, F., Labecca, F., Pedrosa, F., Constantini, R., Emer, C., Silva, W.
R., da Silva, F. R., Ovaskainen, O., & Jordano, P. (2017). Atlantic frugivory:
A plant-frugivore interaction data set for the Atlantic Forest.
Ecology, 98(6), 1729. https://doi.org/10.1002/ecy.1818

## Acknowledgements

This project has received funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (grant agreement No 856506; ERC-synergy project LIFEPLAN)
