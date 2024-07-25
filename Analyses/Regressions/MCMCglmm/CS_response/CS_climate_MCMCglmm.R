############################### BPMMs: MCMCglmm analyses of variation in colony size and climatic variables (TMP PRE DTR) ########################################
####To quantifying phylogenetic uncertainty in parameter estimates, analyses run over a sample of 400 phylogenetic trees
###Written for West group server
##Louis Bell-Roberts
#25/10/2023

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
#Attempt two priors - Cornwallis prior (with and without the intercept), the gelman prior (nations). See if the results differ at all. Ask for verification in an email to Anna or Charlie/Philip

# B is for the fixed effects
# R is for residual variance
# G is the phylogenetic/additive genetic variance
prior1 <- list(G = list(G1 = list(V = 1, nu = 0.002)),
               R = list(V = 1, nu = 0.002))

# prior2 <- list(B = list(mu = rep(0,3), V = diag(3)*(1+pi^2/3)),
#                R = list(V = 1, fix = 1),
#                G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))) #X2 prior
# 
# prior4 <- list(B = list(mu = rep(0,4), V = diag(4)*(1+pi^2/3)),
#                R = list(V = 1, fix = 1),
#                G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))) #Fisher prior



##########
#Global
##########

#Set working directory for location of models to be saved
setwd('/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS/Outputs/Global/1st_run_11M_1M_10k/')
registerDoParallel(50)

#Parallel models where the outputs are save as .rds files
##

# Initialize a list to store model filenames
global_filenames <- character(length(tree_set_global))
# 
# Loop through each element in tree_set_global and perform parallel computations
# foreach(i = seq_along(tree_set_global)) %dopar% {
#   # Create the MCMCglmm model here
#   model <- MCMCglmm(log10(colony.size.y) ~ TMPavg_avg + PREavg_avg + DTRavg_avg, random = ~animal, family = "gaussian", prior = prior1, pedigree = tree_set_global[[i]], data = ant_data_global, nitt = 11000000, burnin = 1000000, thin = 10000, verbose = F)
# 
#   # Generate a unique filename for the model
#   filename <- paste0("CS_climate_global_11M_1M_10k_", i, "_2ndRun.rds")
# 
#   # Save the model as an .rds file
#   saveRDS(model, file = filename)
# 
#   # Store the filename in the list
#   global_filenames[i] <- filename
# }



##########
#Tropical
##########

#Set working directory for location of models to be saved
setwd('/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS/Outputs/Tropical/2nd_run_11M_1M_10k/')
registerDoParallel(50)

#Parallel models where the outputs are save as .rds files
##

# Initialize a list to store model filenames
tropical_filenames <- character(length(tree_set_tropical))

# Loop through each element in tree_set_tropical and perform parallel computations
# foreach(i = seq_along(tree_set_tropical)) %dopar% {
#   # Create the MCMCglmm model here
#   model <- MCMCglmm(log10(colony.size.y) ~ TMPavg_avg + PREavg_avg + DTRavg_avg, random = ~animal, family = "gaussian", prior = prior1, pedigree = tree_set_tropical[[i]], data = ant_data_tropical, nitt = 11000000, burnin = 1000000, thin = 10000, verbose = F)
#   
#   # Generate a unique filename for the model
#   filename <- paste0("CS_climate_tropical_11M_1M_10k_", i, "_2ndRun.rds")
#   
#   # Save the model as an .rds file
#   saveRDS(model, file = filename)
#   
#   # Store the filename in the list
#   tropical_filenames[i] <- filename
# }

#Read in the MCMCglmm objects

# # get list of all RDS files in the working directory
# model_files_CS_tropical <- list.files(pattern = "\\.rds$")
# # read in all models using lapply()
# mcmc_list_CS_tropical <- lapply(model_files_CS_tropical, readRDS)


##########
#Temperate
##########

#Set working directory for location of models to be saved
setwd('/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS/Outputs/Temperate/2nd_run_11M_1M_10k/')
registerDoParallel(50)

#Parallel models where the outputs are save as .rds files
##

# Initialize a list to store model filenames
temperate_filenames <- character(length(tree_set_temperate))

# Loop through each element in tree_set_temperate and perform parallel computations
# foreach(i = seq_along(tree_set_temperate)) %dopar% {
#   # Create the MCMCglmm model here
#   model <- MCMCglmm(log10(colony.size.y) ~ TMPavg_avg + PREavg_avg + DTRavg_avg, random = ~animal, family = "gaussian", prior = prior1, pedigree = tree_set_temperate[[i]], data = ant_data_temperate, nitt = 11000000, burnin = 1000000, thin = 10000, verbose = F)
#   
#   # Generate a unique filename for the model
#   filename <- paste0("CS_climate_temperate_11M_1M_10k_", i, "_2ndRun.rds")
#   
#   # Save the model as an .rds file
#   saveRDS(model, file = filename)
#   
#   # Store the filename in the list
#   temperate_filenames[i] <- filename
# }

#Read in the MCMCglmm objects

# # get list of all RDS files in the working directory
# model_files_CS_temperate <- list.files(pattern = "\\.rds$")
# # read in all models using lapply()
# mcmc_list_CS_temperate <- lapply(model_files_CS_temperate, readRDS)


##########
##Species existing in both tropical and temperate regions (Both)
##########

#Set working directory for location of models to be saved
setwd('/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS/Outputs/Both/2nd_run_11M_1M_10k/')
registerDoParallel(50)

#Parallel models where the outputs are save as .rds files
##

# Initialize a list to store model filenames
both_filenames <- character(length(tree_set_both))

# Loop through each element in tree_set_both and perform parallel computations
# foreach(i = seq_along(tree_set_both)) %dopar% {
#   # Create the MCMCglmm model here
#   model <- MCMCglmm(log10(colony.size.y) ~ TMPavg_avg + PREavg_avg + DTRavg_avg, random = ~animal, family = "gaussian", prior = prior1, pedigree = tree_set_both[[i]], data = ant_data_both, nitt = 11000000, burnin = 1000000, thin = 10000, verbose = F)
#   
#   # Generate a unique filename for the model
#   filename <- paste0("CS_climate_both_11M_1M_10k_", i, "_2ndRun.rds")
#   
#   # Save the model as an .rds file
#   saveRDS(model, file = filename)
#   
#   # Store the filename in the list
#   both_filenames[i] <- filename
# }

#Read in the MCMCglmm objects

# # get list of all RDS files in the working directory
# model_files_CS_both <- list.files(pattern = "\\.rds$")
# # read in all models using lapply()
# mcmc_list_CS_both <- lapply(model_files_CS_both, readRDS)

############ END ##########################


