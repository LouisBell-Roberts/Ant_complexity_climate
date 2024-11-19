########################Calculate measure of phylogenetic signal (heritability) - Climatic predictors and ant traits######################
###Calculate over 400 trees using MCMCglmm
##Run MCMCglmm models in parallel on the West server
#Louis Bell-Roberts
#07/08/2023


.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))

#Load packages
library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
library(geiger)
library(doParallel)
library(MCMCglmm)
library(coda)
library(gtools)

registerDoParallel(40)


#Set variables so that they're in the correct structure and remove species that are not on Antcat. Assign species with NA in Antcat column with 0 (registered)
pre_data <- read.csv("/drives/4tb/Louis/AntEnv/Ant_data/CS/New_transformations/AntClm_CSize_all_230710.csv") #After species with unusual biology removed = 474 species

#Read in ant trees
ant_trees <-read.tree("/drives/4tb/Louis/AntEnv/Ant_trees/CS/400trees_pruned_CSize_230712.tre") 

pre_data$valid_species_name <- gsub("_", ".", pre_data$valid_species_name, fixed = T)
is.ultrametric(ant_trees) # true
is.binary(ant_trees)

#Rename valid_species_name column
pre_data <- pre_data %>% 
  rename(animal=valid_species_name,
         PREavg=PREavg_avg,
         TMPsd=TMPsd_avg)
############## End of generic code for all analyses #############



########################################## 

###################### 
#Signal in TMPsd
###################### 

########################################## 


#Select data
TMP_data_select <- dplyr::select(pre_data, animal, TMPsd)

#########
#MCMCglmm models
#########

##Set prior
# B is for the fixed effect to help with mixing
# R is for residual variance
# G is the phylogenetic/additive genetic variance
prior1 <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)


##############################
#Running models in parallel
##############################

#Choose a location to save MCMCglmm models to

setwd("/drives/4tb/Louis/AntEnv/Analyses/Phylogenetic_signal/Outputs/TMPsd/1M/1st_run/") #Create a new folder to receive the model runs for each variable

#Parallel models where the outputs are save as .rds files

# Initialize a list to store model filenames
model_filenames_TMPsd <- character(length(ant_trees))

# Loop through each element in ant_trees and perform parallel computations
foreach(i = seq_along(ant_trees)) %dopar% {
  # Create the MCMCglmm model here
  model <- MCMCglmm(sqrt(TMPsd) ~ 1, random = ~animal, family = "gaussian", prior = prior1, pedigree = ant_trees[[i]], data = TMP_data_select, nitt = 1100000, burnin = 100000, thin = 1000, verbose = FALSE)

  # Generate a unique filename for the model
  filename <- paste0("MCMC_phylo_sig_model_TMPsd_1M_", i, ".rds")

  # Save the model as an .rds file
  saveRDS(model, file = filename)

  # Store the filename in the list
  model_filenames_TMPsd[i] <- filename
}

setwd("/drives/4tb/Louis/AntEnv/Analyses/Phylogenetic_signal/Outputs/TMPsd/1M/2nd_run/") #Create a new folder to receive the model runs for each variable

foreach(i = seq_along(ant_trees)) %dopar% {
  # Create the MCMCglmm model here
  model <- MCMCglmm(sqrt(TMPsd) ~ 1, random = ~animal, family = "gaussian", prior = prior1, pedigree = ant_trees[[i]], data = TMP_data_select, nitt = 1100000, burnin = 100000, thin = 1000, verbose = FALSE)
  
  # Generate a unique filename for the model
  filename <- paste0("MCMC_phylo_sig_model_TMPsd_1M_", i, ".rds")
  
  # Save the model as an .rds file
  saveRDS(model, file = filename)
  
  # Store the filename in the list
  model_filenames_TMPsd[i] <- filename
}

#Read in the MCMCglmm objects

# get list of all RDS files in the working directory
model_files_TMPsd <- list.files(pattern = "\\.rds$")
# read in all models using lapply()
mcmc_list_TMPsd <- lapply(model_files_TMPsd, readRDS)

#Check effective size for each model seperately
effectiveSize(mcmc_list_TMPsd[[4]]$VCV)


##############################
#Combine posterior distribution of the random effects - VCV
##############################

# extract the VCV objects from each model output in the list and create a list of the mcmc objects
vcv_matrices_TMPsd <- lapply(mcmc_list_TMPsd, function(model) model$VCV)
#Combine the mcmc objects together
combined_vcv_TMPsd <- do.call(rbind, vcv_matrices_TMPsd)
#Assign as an mcmc object
MCMC_combined_TMPsd <- as.mcmc(combined_vcv_TMPsd)


#Calculate phylogenetic signal (heritability)
herit_TMPsd <- MCMC_combined_TMPsd[ , "animal"] / (MCMC_combined_TMPsd[ , "animal"] + MCMC_combined_TMPsd[ , "units"])

effectiveSize(herit_TMPsd)

posterior.mode(herit_TMPsd) #0.7433347; CI = 0.622812 to 0.8593536

#Compute 95% credible interval
HPDinterval(herit_TMPsd)

#Posterior mode of additive genetic variance
posterior.mode(MCMC_combined_TMPsd[, "animal"]) #7.470811; CI = 5.020869 to 11.09869
#95% credible interval of additive genetic variance
HPDinterval(MCMC_combined_TMPsd[, "animal"])

#Posterior mode of additive genetic variance
posterior.mode(MCMC_combined_TMPsd[, "units"]) #2.625827; CI = 1.794206 to 3.421187
#95% credible interval of additive genetic variance
HPDinterval(MCMC_combined_TMPsd[, "units"])



########################################## 

###################### 
#Signal in PREavg
###################### 

########################################## 


#Select data
PRE_data_select <- dplyr::select(pre_data, animal, PREavg)

#########
#MCMCglmm models
#########

##Set prior
# B is for the fixed effect to help with mixing
# R is for residual variance
# G is the phylogenetic/additive genetic variance
prior1 <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)


##############################
#Running models in parallel
##############################

#Choose a location to save MCMCglmm models to

setwd("/drives/4tb/Louis/AntEnv/Analyses/Phylogenetic_signal/Outputs/PREavg/1M/1st_run/")

#Parallel models where the outputs are save as .rds files

# Initialize a list to store model filenames
model_filenames_PREavg <- character(length(ant_trees))

# Loop through each element in phylo_test_list_CS and perform parallel computations
foreach(i = seq_along(ant_trees)) %dopar% {
  # Create the MCMCglmm model here
  model <- MCMCglmm(sqrt(PREavg) ~ 1, random = ~animal, family = "gaussian", prior = prior1, pedigree = ant_trees[[i]], data = PRE_data_select, nitt = 1100000, burnin = 100000, thin = 1000, verbose = FALSE)

  # Generate a unique filename for the model
  filename <- paste0("MCMC_phylo_sig_model_PREavg_1M_", i, ".rds")

  # Save the model as an .rds file
  saveRDS(model, file = filename)

  # Store the filename in the list
  model_filenames_PREavg[i] <- filename
}

setwd("/drives/4tb/Louis/AntEnv/Analyses/Phylogenetic_signal/Outputs/PREavg/1M/2nd_run/")

foreach(i = seq_along(ant_trees)) %dopar% {
  # Create the MCMCglmm model here
  model <- MCMCglmm(sqrt(PREavg) ~ 1, random = ~animal, family = "gaussian", prior = prior1, pedigree = ant_trees[[i]], data = PRE_data_select, nitt = 1100000, burnin = 100000, thin = 1000, verbose = FALSE)
  
  # Generate a unique filename for the model
  filename <- paste0("MCMC_phylo_sig_model_PREavg_1M_", i, ".rds")
  
  # Save the model as an .rds file
  saveRDS(model, file = filename)
  
  # Store the filename in the list
  model_filenames_PREavg[i] <- filename
}

#Read in the MCMCglmm objects

# get list of all RDS files in the working directory
model_files_PREavg <- list.files(pattern = "\\.rds$")
# read in all models using lapply()
mcmc_list_PREavg <- lapply(model_files_PREavg, readRDS)

#Check effective size for each model seperately
effectiveSize(mcmc_list_PREavg[[1]]$VCV)


##############################
#Combine posterior distribution of the random effects - VCV
##############################

# extract the VCV objects from each model output in the list and create a list of the mcmc objects
vcv_matrices_PREavg <- lapply(mcmc_list_PREavg, function(model) model$VCV)
#Combine the mcmc objects together
combined_vcv_PREavg <- do.call(rbind, vcv_matrices_PREavg)
#Assign as an mcmc object
MCMC_combined_PREavg <- as.mcmc(combined_vcv_PREavg)


#Calculate phylogenetic signal (heritability)
herit_PREavg <- MCMC_combined_PREavg[ , "animal"] / (MCMC_combined_PREavg[ , "animal"] + MCMC_combined_PREavg[ , "units"])

effectiveSize(herit_PREavg)

posterior.mode(herit_PREavg) #0.7103726; CI = 0.5068401 to 0.9045107

#Compute 95% credible interval
HPDinterval(herit_PREavg)

#Posterior mode of additive genetic variance
posterior.mode(MCMC_combined_PREavg[, "animal"]) #4319.562; CI = 2298.54 to 8708.851
#95% credible interval of additive genetic variance
HPDinterval(MCMC_combined_PREavg[, "animal"])

#Posterior mode of additive genetic variance
posterior.mode(MCMC_combined_PREavg[, "units"]) #2040.77; CI = 1097.923 to 2863.637
#95% credible interval of additive genetic variance
HPDinterval(MCMC_combined_PREavg[, "units"])







########################################## 

###################### 
#Signal in Colony size
###################### 

########################################## 


#Select data
CS_data_select <- dplyr::select(pre_data, animal, colony.size.y)

#########
#MCMCglmm models
#########

##Set prior
# B is for the fixed effect to help with mixing
# R is for residual variance
# G is the phylogenetic/additive genetic variance
prior1 <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)


##############################
#Running models in parallel
##############################

#Choose a location to save MCMCglmm models to

setwd("/drives/4tb/Louis/AntEnv/Analyses/Phylogenetic_signal/Outputs/Colony_size/1M/1st_run/")

#Parallel models where the outputs are save as .rds files

# Initialize a list to store model filenames
model_filenames_CS <- character(length(ant_trees))

# Loop through each element in phylo_test_list_CS and perform parallel computations
foreach(i = seq_along(ant_trees)) %dopar% {
  # Create the MCMCglmm model here
  model <- MCMCglmm(log10(colony.size.y) ~ 1, random = ~animal, family = "gaussian", prior = prior1, pedigree = ant_trees[[i]], data = CS_data_select, nitt = 1100000, burnin = 100000, thin = 1000, verbose = FALSE)
  
  # Generate a unique filename for the model
  filename <- paste0("MCMC_phylo_sig_model_CS_1M_", i, ".rds")
  
  # Save the model as an .rds file
  saveRDS(model, file = filename)
  
  # Store the filename in the list
  model_filenames_CS[i] <- filename
}

setwd("/drives/4tb/Louis/AntEnv/Analyses/Phylogenetic_signal/Outputs/Colony_size/1M/2nd_run/")

foreach(i = seq_along(ant_trees)) %dopar% {
  # Create the MCMCglmm model here
  model <- MCMCglmm(log10(colony.size.y) ~ 1, random = ~animal, family = "gaussian", prior = prior1, pedigree = ant_trees[[i]], data = CS_data_select, nitt = 1100000, burnin = 100000, thin = 1000, verbose = FALSE)
  
  # Generate a unique filename for the model
  filename <- paste0("MCMC_phylo_sig_model_CS_1M_", i, ".rds")
  
  # Save the model as an .rds file
  saveRDS(model, file = filename)
  
  # Store the filename in the list
  model_filenames_CS[i] <- filename
}


#Read in the MCMCglmm objects
setwd("/drives/4tb/Louis/AntEnv/Analyses/Phylogenetic_signal/Outputs/Colony_size/1M/1st_run/")

# get list of all RDS files in the working directory
model_files_CS <- list.files(pattern = "\\.rds$")
# Sort the files by name
sorted_model_files_CS <- mixedsort(model_files_CS)

# read in all models using lapply()
mcmc_list_CS <- lapply(sorted_model_files_CS, readRDS)

#Check effective size for each model seperately
effectiveSize(mcmc_list_CS[[1]]$VCV)


##############################
#Combine posterior distribution of the random effects - VCV
##############################

# extract the VCV objects from each model output in the list and create a list of the mcmc objects
vcv_matrices_CS <- lapply(mcmc_list_CS, function(model) model$VCV)
#Combine the mcmc objects together
combined_vcv_CS <- do.call(rbind, vcv_matrices_CS)
#Assign as an mcmc object
MCMC_combined_CS <- as.mcmc(combined_vcv_CS)


#Calculate phylogenetic signal (heritability)
herit_CS <- MCMC_combined_CS[ , "animal"] / (MCMC_combined_CS[ , "animal"] + MCMC_combined_CS[ , "units"])

effectiveSize(herit_CS)

posterior.mode(herit_CS) #0.7225756; CI = 0.5812703 to 0.8388493

#Compute 95% credible interval
HPDinterval(herit_CS)

#Posterior mode of additive genetic variance
posterior.mode(MCMC_combined_CS[, "animal"]) #1.191457; CI = 0.7809104 to 1.837078
#95% credible interval of additive genetic variance
HPDinterval(MCMC_combined_CS[, "animal"])

#Posterior mode of additive genetic variance
posterior.mode(MCMC_combined_CS[, "units"]) #0.4935909; CI = 0.3504974 to 0.6430719
#95% credible interval of additive genetic variance
HPDinterval(MCMC_combined_CS[, "units"])








########################################## 

###################### 
#Signal in Worker polymorphism
######################

#Since worker polymorphism is a binary trait, the residual variance is not identifiable (fix = 1) and so it's set to 1. We estimate the amount of variation in worker polymorphism explained by shared ancestry between species calculated as the intraclass correlation coefficient (phylogenetic variance/ (phylogenetic variance + 1) + π^(2/3)) for worker polymorphism — an analogous measure of the amount of variation explained by phylogenetic history appropriate for binary traits (Agreed by Villemereuil in his tutorial).

########################################## 


#Select data
WP_data_select <- dplyr::select(pre_data, animal, poly_id)

#########
#MCMCglmm models
#########

##Set prior
# B is for the fixed effect to help with mixing
# R is for residual variance
# G is the phylogenetic/additive genetic variance
prior_X2 <- list(G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)),
               R = list(V = 1, fix = 1)) #X2 random effect


##############################
#Running models in parallel
##############################

#Choose a location to save MCMCglmm models to

setwd("/drives/4tb/Louis/AntEnv/Analyses/Phylogenetic_signal/Outputs/Worker_polymorphism/1M/1st_run/")

#Parallel models where the outputs are save as .rds files

# Initialize a list to store model filenames
model_filenames_WP <- character(length(ant_trees))

# Loop through each element in phylo_test_list_WP and perform parallel computations
foreach(i = seq_along(ant_trees)) %dopar% {
  # Create the MCMCglmm model here
  model <- MCMCglmm(poly_id ~ 1, random = ~animal, family = "categorical", prior = prior_X2, pedigree = ant_trees[[i]], data = WP_data_select, nitt = 1100000, burnin = 100000, thin = 1000, verbose = FALSE)
  
  # Generate a unique filename for the model
  filename <- paste0("MCMC_phylo_sig_model_WP_1M_", i, ".rds")
  
  # Save the model as an .rds file
  saveRDS(model, file = filename)
  
  # Store the filename in the list
  model_filenames_WP[i] <- filename
}

setwd("/drives/4tb/Louis/AntEnv/Analyses/Phylogenetic_signal/Outputs/Worker_polymorphism/1M/2nd_run/")

foreach(i = seq_along(ant_trees)) %dopar% {
  # Create the MCMCglmm model here
  model <- MCMCglmm(poly_id ~ 1, random = ~animal, family = "categorical", prior = prior_X2, pedigree = ant_trees[[i]], data = WP_data_select, nitt = 1100000, burnin = 100000, thin = 1000, verbose = FALSE)
  
  # Generate a unique filename for the model
  filename <- paste0("MCMC_phylo_sig_model_WP_1M_", i, ".rds")
  
  # Save the model as an .rds file
  saveRDS(model, file = filename)
  
  # Store the filename in the list
  model_filenames_WP[i] <- filename
}


#Read in the MCMCglmm objects
setwd("/drives/4tb/Louis/AntEnv/Analyses/Phylogenetic_signal/Outputs/Worker_polymorphism/1M/1st_run/")

# get list of all RDS files in the working directory
model_files_WP <- list.files(pattern = "\\.rds$")

# Sort the files by name
sorted_model_files_WP <- mixedsort(model_files_WP)

# read in all models using lapply()
mcmc_list_WP <- lapply(sorted_model_files_WP, readRDS)

#Check effective size for each model seperately
effectiveSize(mcmc_list_WP[[1]]$VCV)


##############################
#Combine posterior distribution of the random effects - VCV
##############################

# extract the VCV objects from each model output in the list and create a list of the mcmc objects
vcv_matrices_WP <- lapply(mcmc_list_WP, function(model) model$VCV)
#Combine the mcmc objects together
combined_vcv_WP <- do.call(rbind, vcv_matrices_WP)
#Assign as an mcmc object
MCMC_combined_WP <- as.mcmc(combined_vcv_WP)


#Calculate phylogenetic signal (heritability)
herit_WP <- MCMC_combined_WP[ , "animal"] / (MCMC_combined_WP[ , "animal"] + 1 + (pi^2)/3)

effectiveSize(herit_WP)

posterior.mode(herit_WP) #0.7701589; CI = 0.6476877 to 0.8543225

#Compute 95% credible interval
HPDinterval(herit_WP)

#Posterior mode of additive genetic variance
posterior.mode(MCMC_combined_WP[, "animal"]) #12.87283; CI = 6.74413 to 22.57806
#95% credible interval of additive genetic variance
HPDinterval(MCMC_combined_WP[, "animal"])

#Posterior mode of additive genetic variance
posterior.mode(MCMC_combined_WP[, "units"]) #1; CI = 1 to 1
#95% credible interval of additive genetic variance
HPDinterval(MCMC_combined_WP[, "units"])
