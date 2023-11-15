# African-Frugivory
- Plant trait and phylogeny much better known in Madagascar
- Plant ocucrrence (plot maps) and elephant density are available for Gabon

## Application 1: Predicting unobserved interactions
- Predict non-observed interactions for the purpose of estimating the effects of defaunation on above ground biomass
- Existing estimates fail to account for unobserved interactions and hence underestimate deforestation
- As a first step, ignore the spatial component and environmental covariates
- Phylogenetic distance is available for most animals
- Phylogenetic distance is less available for plants, likely only available at higher levels. Consider allowing rho to vary by genus, order, etc. 

Open questions/notes
- Do we want to increase the number of traits even if there may be a lot of missing values
- Some animal trates are at the genus level
- Can we incorporate plant IUCN status/threat level as a proxy for density. Is this circular? Can we put it into the detection submodel?
- Can we infer non-frugivory interactions if we observed an herbivory interaction on the same tree?

Special opportunity with Madagascar data: more detailed interaction networks will be created, combining field observations, ecological knowledge, surveys of local communities. We could use this data to validate the literature-based predicted network. 


## Application 2: Spatial Extensions
- Incorporate spatial information, predict networks at new locations
- Add home range size

## Application 3: Predicting interactions with extinct species
- Motivation: predict lost ecosystem services, e.g. carbon storage, how did extinct lemurs play a role in disbursal of ethnobotanically important plants
- For extinct lemurs, we have limited trait data, e.g., body size


## Application 4: Incorporate occurrence data (Gabon)
- Preliminary evidence suggests that elephant density is associated with greater above ground biomass
- Effects of trampling versus seed dispersal
