############################### BPMMs: MCMCglmm analyses of variation in worker size (monomorphic/polymorphic) and colony size ########################################
####To quantifying phylogenetic uncertainty in parameter estimates, analyses run over a sample of 400 phylogenetic trees
###Written for West group server
##Louis Bell-Roberts
#18/03/2024


.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))

#############             #############               #############             

library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
library(MCMCglmm)
library(arm)
library(lme4)
library(doParallel)
library(coda)

registerDoParallel(50)


#Read in ant data
##Global
ant_data_global <- read.csv("/drives/4tb/Louis/AntEnv/Ant_data/CS/New_transformations/AntClm_CSize_all_230710.csv") #474 species
ant_data_global <- ant_data_global %>% rename(animal = valid_species_name)
##Tropical
ant_data_tropical <- read.csv("/drives/4tb/Louis/AntEnv/Ant_data/CS/New_transformations/Tropical_ant_data.csv") #116 species
ant_data_tropical <- ant_data_tropical %>% rename(animal = valid_species_name)
##Temperate
ant_data_temperate <- read.csv("/drives/4tb/Louis/AntEnv/Ant_data/CS/New_transformations/Temperate_ant_data_230804.csv", ) #190 species
ant_data_temperate <- ant_data_temperate %>% rename(animal = valid_species_name)
#Species existing in both tropical and temperate regions
ant_data_both <- read.csv("/drives/4tb/Louis/AntEnv/Ant_data/CS/New_transformations/Both_ant_data_230804.csv") #168 species
ant_data_both <- ant_data_both %>% rename(animal = valid_species_name)

#Read in phylogenies
##Global
tree_set_global<- read.tree("/drives/4tb/Louis/AntEnv/Ant_trees/CS/400trees_pruned_CSize_230712.tre")
##Tropical
tree_set_tropical<- read.tree("/drives/4tb/Louis/AntEnv/Ant_trees/CS/Tropical_ant_trees_pruned.tre")
#Temperate
tree_set_temperate<- read.tree("/drives/4tb/Louis/AntEnv/Ant_trees/CS/Temperate_ant_trees_pruned.tre")
#Species existing in both tropical and temperate regions
tree_set_both<- read.tree("/drives/4tb/Louis/AntEnv/Ant_trees/CS/Both_ant_trees_pruned.tre")

##Set prior
# R is for residual variance
# G is the phylogenetic/additive genetic variance
prior1 <- list(G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)), #Chi-squared prior
               R = list(V = 1, fix = 1)) #Residual variance fixed to 1


##########
#Global
##########

#Set working directory for location of models to be saved
setwd('/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS_Caste/Outputs/Global/')

#Parallel models where the outputs are save as .rds files
##

#Loop through each element in tree_set_global and perform parallel computations
foreach(i = seq_along(tree_set_global)) %dopar% {
  # Create the MCMCglmm model here
  model1 <- MCMCglmm(poly_id ~ log10(colony.size.y), random = ~animal, family = "categorical", prior = prior1, pedigree = tree_set_global[[i]], data = ant_data_global, nitt = 11000000, burnin = 1000000, thin = 10000, verbose = F)
  model2 <- MCMCglmm(poly_id ~ log10(colony.size.y), random = ~animal, family = "categorical", prior = prior1, pedigree = tree_set_global[[i]], data = ant_data_global, nitt = 11000000, burnin = 1000000, thin = 10000, verbose = F)
  
  # Generate a unique filename for the model
  filename1 <- paste0("1st_run/CS_Caste_global_11M_1M_10k_", i, "_1stRun.rds")
  filename2 <- paste0("2nd_run/CS_Caste_global_11M_1M_10k_", i, "_2ndRun.rds")
  
  # Save the model as an .rds file
  saveRDS(model1, file = filename1)
  saveRDS(model2, file = filename2)
  
}



##########
#Tropical
##########

#Set working directory for location of models to be saved
setwd('/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS_Caste/Outputs/Tropical/')

#Parallel models where the outputs are save as .rds files
##

#Loop through each element in tree_set_tropical and perform parallel computations
foreach(i = seq_along(tree_set_tropical)) %dopar% {
  # Create the MCMCglmm model here
  model1 <- MCMCglmm(poly_id ~ log10(colony.size.y), random = ~animal, family = "categorical", prior = prior1, pedigree = tree_set_tropical[[i]], data = ant_data_tropical, nitt = 11000000, burnin = 1000000, thin = 10000, verbose = F)
  model2 <- MCMCglmm(poly_id ~ log10(colony.size.y), random = ~animal, family = "categorical", prior = prior1, pedigree = tree_set_tropical[[i]], data = ant_data_tropical, nitt = 11000000, burnin = 1000000, thin = 10000, verbose = F)
  
  # Generate a unique filename for the model
  filename1 <- paste0("1st_run/CS_Caste_tropical_11M_1M_10k_", i, "_1stRun.rds")
  filename2 <- paste0("2nd_run/CS_Caste_tropical_11M_1M_10k_", i, "_2ndRun.rds")
  
  # Save the model as an .rds file
  saveRDS(model1, file = filename1)
  saveRDS(model2, file = filename2)
  
}



##########
#Temperate
##########

#Set working directory for location of models to be saved
setwd('/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS_Caste/Outputs/Temperate/')

#Parallel models where the outputs are save as .rds files
##

# Loop through each element in tree_set_temperate and perform parallel computations
foreach(i = seq_along(tree_set_temperate)) %dopar% {
  # Create the MCMCglmm model here
  model1 <- MCMCglmm(poly_id ~ log10(colony.size.y), random = ~animal, family = "categorical", prior = prior1, pedigree = tree_set_temperate[[i]], data = ant_data_temperate, nitt = 11000000, burnin = 1000000, thin = 10000, verbose = F)
  model2 <- MCMCglmm(poly_id ~ log10(colony.size.y), random = ~animal, family = "categorical", prior = prior1, pedigree = tree_set_temperate[[i]], data = ant_data_temperate, nitt = 11000000, burnin = 1000000, thin = 10000, verbose = F)
  
  # Generate a unique filename for the model
  filename1 <- paste0("1st_run/CS_Caste_temperate_11M_1M_10k_", i, "_1stRun.rds")
  filename2 <- paste0("2nd_run/CS_Caste_temperate_11M_1M_10k_", i, "_2ndRun.rds")
  
  # Save the model as an .rds file
  saveRDS(model1, file = filename1)
  saveRDS(model2, file = filename2)
  
}



##########
##Species existing in both tropical and temperate regions (Both)
##########

#Set working directory for location of models to be saved
setwd('/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS_Caste/Outputs/Both/')

#Parallel models where the outputs are save as .rds files
##

# Loop through each element in tree_set_both and perform parallel computations
foreach(i = seq_along(tree_set_both)) %dopar% {
  # Create the MCMCglmm model here
  model1 <- MCMCglmm(poly_id ~ log10(colony.size.y), random = ~animal, family = "categorical", prior = prior1, pedigree = tree_set_both[[i]], data = ant_data_both, nitt = 11000000, burnin = 1000000, thin = 10000, verbose = F)
  model2 <- MCMCglmm(poly_id ~ log10(colony.size.y), random = ~animal, family = "categorical", prior = prior1, pedigree = tree_set_both[[i]], data = ant_data_both, nitt = 11000000, burnin = 1000000, thin = 10000, verbose = F)
  
  # Generate a unique filename for the model
  filename1 <- paste0("1st_run/CS_Caste_both_11M_1M_10k_", i, "_1stRun.rds")
  filename2 <- paste0("2nd_run/CS_Caste_both_11M_1M_10k_", i, "_2ndRun.rds")
  
  # Save the model as an .rds file
  saveRDS(model1, file = filename1)
  saveRDS(model2, file = filename2)
  
}


############ END ##########################


