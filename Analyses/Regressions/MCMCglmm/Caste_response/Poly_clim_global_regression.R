############################### BPMMs: MCMCglmm analyses of variation in worker size (monomorphic/polymorphic) and climatic variables (TMP PRE DTR) ########################################
####To quantifying phylogenetic uncertainty in parameter estimates, analyses run over a sample of 400 phylogenetic trees
###Written for West group server
##Louis Bell-Roberts
#25/10/2023


#############             #############               #############             

library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
library(geiger)
library(MuMIn)
library(MCMCglmm)
library(arm)
library(lme4)
library(purrr)


#Read in ant data
##Global
ant_data<- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/AntClm_CSize_all_230710.csv")
ant_data <- ant_data %>% dplyr::rename(animal = valid_species_name)

##Tropical
##Temperate
#Species existing in both tropical and temperate regions

#Read in phylogenies
##Global
tree_set<- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/400_trees/400trees_pruned_CSize_230712.tre") 
##Tropical
#Temperate
#Species existing in both tropical and temperate regions

##Set prior
#Attempt two priors - Cornwallis prior (with and without the intercept), the gelman prior (nations). See if the results differ at all. Ask for verification in an email to Anna or Charlie/Philip

# B is for the fixed effects
# R is for residual variance
# G is the phylogenetic/additive genetic variance
prior1 <- list(G = list(G1 = list(V = 1, nu = 0.002)),
               R = list(V = 1, fix = 1))

prior1.1 <- list(G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)),
               R = list(V = 1, fix = 1))

prior1.2 <- list(B = list(mu = rep(0,4), V = diag(4)*(1)),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)),
                 R = list(V = 1, fix = 1))

prior2 <- list(B = list(mu = rep(0,4), V = diag(4)*(1+pi^2/3)),
               R = list(V = 1, fix = 1),
               G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))) #X2 prior

prior3 <- list(B = list(mu = rep(0,4), V = gelman.prior(~ TMPavg_avg + PREavg_avg + DTRavg_avg, data = ant_data, scale = sqrt(1+pi^2/3))),
               R = list(V = 1, fix = 1),
               G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

prior3.1 <- list(B = list(mu = rep(0,4), V = gelman.prior(~ TMPavg_avg + PREavg_avg + DTRavg_avg, data = ant_data, scale = sqrt(1+pi^2/3))),
               R = list(V = 1, fix = 1),
               G = list(G1 = list(V = diag(1)*0.1, nu = 1)))

prior <- list(R = list(V = 1, fix = 1),
              G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))) #Fisher prior


##########
#Global
##########

#Set working directory for location of models to be saved
# setwd("/drives/4tb/Louis/...")

Basic_model <- MCMCglmm(poly_id ~ TMPavg_avg + PREavg_avg + DTRavg_avg, random = ~animal, family = "categorical", prior = prior1, pedigree = tree_set[[1]], data = ant_data, nitt = 1100000, burnin = 100000, thin = 1000, verbose = F) 

Basic_X2_model <- MCMCglmm(poly_id ~ TMPavg_avg + PREavg_avg + DTRavg_avg, random = ~animal, family = "categorical", prior = prior1.1, pedigree = tree_set[[1]], data = ant_data, nitt = 1100000, burnin = 100000, thin = 1000, verbose = F) 

Low_var_model <- MCMCglmm(poly_id ~ TMPavg_avg + PREavg_avg + DTRavg_avg, random = ~animal, family = "categorical", prior = prior1.2, pedigree = tree_set[[1]], data = ant_data, nitt = 1100000, burnin = 100000, thin = 1000, verbose = F) 

Charlie_model <- MCMCglmm(poly_id ~ TMPavg_avg + PREavg_avg + DTRavg_avg, random = ~animal, family = "categorical", prior = prior2, pedigree = tree_set[[1]], data = ant_data, nitt = 1100000, burnin = 100000, thin = 1000, verbose = F) 

Gelman_model <- MCMCglmm(poly_id ~ TMPavg_avg + PREavg_avg + DTRavg_avg, random = ~animal, family = "categorical", prior = prior3, pedigree = tree_set[[1]], data = ant_data, nitt = 1100000, burnin = 100000, thin = 1000, verbose = F) 

Jarrod_model <- MCMCglmm(poly_id ~ TMPavg_avg + PREavg_avg + DTRavg_avg, random = ~animal, family = "categorical", prior = prior3.1, pedigree = tree_set[[1]], data = ant_data, nitt = 1100000, burnin = 100000, thin = 1000, verbose = F) 

summary(Basic_model)
summary(Basic_X2_model)
summary(Low_var_model)
summary(Charlie_model)
summary(Gelman_model)
summary(Jarrod_model)

autocorr.diag(Basic_model$Sol)
autocorr.diag(Basic_model$VCV)
autocorr.diag(Basic_X2_model$Sol)
autocorr.diag(Basic_X2_model$VCV)
autocorr.diag(Charlie_model$Sol)
autocorr.diag(Charlie_model$VCV)
autocorr.diag(Gelman_model$Sol)
autocorr.diag(Gelman_model$VCV)
autocorr.diag(Jarrod_model$Sol)
autocorr.diag(Jarrod_model$VCV)

plot(Basic_model$VCV)
plot(Basic_X2_model$VCV)
plot(Charlie_model$VCV)
plot(Gelman_model$VCV)
plot(Jarrod_model$VCV)

plot(Basic_model$Sol)
plot(Basic_X2_model$Sol)
plot(Charlie_model$Sol)
plot(Gelman_model$Sol)
plot(Jarrod_model$Sol)

#Parallel models where the outputs are save as .rds files
##

# Initialize a list to store model filenames
global_filenames <- character(length(tree_set))

# Loop through each element in tree_set and perform parallel computations
foreach(i = seq_along(tree_set)) %dopar% {
  # Create the MCMCglmm model here
  model <- MCMCglmm(poly_id ~ PREDICTORS, random = ~valid_species_name, family = "categorical", prior = prior2, pedigree = tree_set[[i]], data = ant_data, nitt = 11000, burnin = 1000, thin = 10, verbose = FALSE)
  
  # Generate a unique filename for the model
  filename <- paste0("MCMC_phylo_sig_model_", i, ".rds")
  
  # Save the model as an .rds file
  saveRDS(model, file = filename)
  
  # Store the filename in the list
  global_filenames[i] <- filename
}




##########
#Tropical
##########
mulTree_data_CS_poly <- as.mulTree(data = data_CS_poly, tree = PT_data_CS_poly,
                                   taxa = "animal")

CS_poly <- worker_CoV ~ colony.size
mulTree(mulTree.data = mulTree_data_CS_poly, formula = CS_poly, priors = prior1,
        parameters = mul_parameters1M, output = "CS_poly_1M", ESS = 1000,
        chains = 2, family = "gaussian")


##########
#Temperate
##########
mulTree_data_CS_poly <- as.mulTree(data = data_CS_poly, tree = PT_data_CS_poly,
                                   taxa = "animal")

CS_poly <- worker_CoV ~ colony.size
mulTree(mulTree.data = mulTree_data_CS_poly, formula = CS_poly, priors = prior1,
        parameters = mul_parameters1M, output = "CS_poly_1M", ESS = 1000,
        chains = 2, family = "gaussian")


##########
#Both regions
##########
mulTree_data_CS_poly <- as.mulTree(data = data_CS_poly, tree = PT_data_CS_poly,
                                   taxa = "animal")

CS_poly <- worker_CoV ~ colony.size
mulTree(mulTree.data = mulTree_data_CS_poly, formula = CS_poly, priors = prior1,
        parameters = mul_parameters1M, output = "CS_poly_1M", ESS = 1000,
        chains = 2, family = "gaussian")
############ END ##########################


