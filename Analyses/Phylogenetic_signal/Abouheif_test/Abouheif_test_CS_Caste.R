########## Abouheif test for phylogenetic signal in Caste and colony size ############
##Louis Bell-Roberts
#23/10/2023

library(adephylo)
library(phylobase)
library(geiger)
library(phytools)
library(ade4)
library(ape)
library(dplyr)

#Read in data
ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/AntClm_CSize_all_230710.csv") #After species with unusual biology removed = 474 species
ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/Tropical_ant_data.csv") #After species with unusual biology removed = 474 species
ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/Temperate_ant_data_230804.csv") #After species with unusual biology removed = 474 species
ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/Both_ant_data_230804.csv") #After species with unusual biology removed = 474 species

ant_data2 <- dplyr::select(ant_data, colony.size.y, poly_id)
ant_data2$colony.size.y <- log10(ant_data2$colony.size.y)

#Read in trees
tree_set<- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/400_trees/400trees_pruned_CSize_230712.tre")
tree_set<- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/Tropical_vs_temperate/400_trees/Tropical_ant_trees_pruned.tre")
tree_set<- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/Tropical_vs_temperate/400_trees/Temperate_ant_trees_pruned.tre")
tree_set<- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/Tropical_vs_temperate/400_trees/Both_ant_trees_pruned.tre")

#Combine data and tree
phylotraits <- phylo4d(tree_set[[1]], ant_data2)

#Run test
abouheif.test <- abouheif.moran(phylotraits, method="oriAbouheif")
abouheif.test

plot(abouheif.test)

############## 
# library(phyr)
# library(rr2)
# rownames(ant_data) <- ant_data$valid_species_name
# ref_mod <- lm(log10(colony.size.y) ~ poly_id, data = ant_data)
# CS_Caste <- pglmm_compare(log10(colony.size.y) ~ poly_id, family = "gaussian", phy = tree_set[[1]], data = ant_data, REML = T)
# summary(CS_Caste)
# Int_mod <- pglmm_compare(log10(colony.size.y) ~ 1, family = "gaussian", phy = tree_set[[1]], data = ant_data, REML = T)
# 
# R2(CS_Caste) #Total variance explained
# R2(mod = CS_Caste, mod.r = ref_mod) #Explained by phy
# R2(mod = CS_Caste, mod.r = Int_mod) #Explained by fixed effects





