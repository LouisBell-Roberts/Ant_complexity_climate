############################### BPMMs: Parameter estimation for MCMCglmm analyses of variation in colony size and climatic variables (TMP PRE DTR) ########################################
####To quantifying phylogenetic uncertainty in parameter estimates, analyses run over a sample of 400 phylogenetic trees
###Written for Louis' laptop
##Louis Bell-Roberts
#30/10/2023


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


## Function for reading in MCMCglmm models and combining model outputs for subsequent use in parameter estimation

combine_model_outputs <- function(directory = getwd(), file_pattern = "\\.rds") {
  # Get list of all RDS files in the specified directory
  model_files <- list.files(path = directory, pattern = file_pattern, full.names = TRUE)
  
  # Read in all models using lapply()
  mcmc_list <- lapply(model_files, readRDS)
  
  # Extract the Sol objects from each model output and create a list of the mcmc objects
  sol_matrices <- lapply(mcmc_list, function(model) model$Sol)
  
  # Combine the mcmc objects together
  combined_sol <- do.call(rbind, sol_matrices)
  
  # Assign as an mcmc object
  MCMC_combined <- as.mcmc(combined_sol)
  
  return(MCMC_combined)
}


##################################################


##########
#Global
##########

#When using the combine_model_outputs function - choose the directory folder containing MCMCglmm objects
MCMC_combined_CS_global <- combine_model_outputs(directory = "/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/CS_response/Outputs/Global/1st_run_11M_1M_10k/", file_pattern = "\\.rds")
length(MCMC_combined_CS_global[,1]) #Should be 400000

#Summary outputs
effectiveSize(MCMC_combined_CS_global)
posterior.mode(MCMC_combined_CS_global) #2.940924640; 0.009696118; -0.002150927; -0.027389980

#Compute 95% credible interval
HPDinterval(MCMC_combined_CS_global) #Only PRE significant: 1.968556000 to 4.0821803418; -0.013402148 to 0.0331984928; -0.004071438 to -0.0004697124; -0.081838881 to 0.0233957048


###
#This code was turned into a function - therefore, now redundant
###

# #Read in the MCMCglmm objects
# 
# # get list of all RDS files in the working directory
# model_files_CS_global <- list.files(pattern = "\\.rds$")
# # read in all models using lapply()
# mcmc_list_CS_global <- lapply(model_files_CS_global, readRDS)
# effectiveSize(mcmc_list_CS_global[[201]]$Sol) #Check the effective sample size for each of the fixed effects
# 
# #Combine posterior distribution of the fixed effects - Sol
# 
# # extract the Sol objects from each model output in the list and create a list of the mcmc objects
# sol_matrices_CS_global <- lapply(mcmc_list_CS_global, function(model) model$Sol)
# #Combine the mcmc objects together
# combined_sol_CS_global <- do.call(rbind, sol_matrices_CS_global)
# #Assign as an mcmc object
# MCMC_combined_CS_global <- as.mcmc(combined_sol_CS_global)
# 



##########
#Tropical
##########
#When using the combine_model_outputs function - choose the directory folder containing MCMCglmm objects
MCMC_combined_CS_tropical <- combine_model_outputs(directory = "/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/CS_response/Outputs/Tropical/1st_run_11M_1M_10k/", file_pattern = "\\.rds")
length(MCMC_combined_CS_tropical[,1]) #Should be 400000

#Summary outputs
effectiveSize(MCMC_combined_CS_tropical)
posterior.mode(MCMC_combined_CS_tropical) #2.086836221; 0.027301367; -0.001622869; -0.067832378

#Compute 95% credible interval
HPDinterval(MCMC_combined_CS_tropical) #None significant: -3.758841945 to 7.943346925; -0.173435932 to 0.216132191; -0.004512138 to 0.001290538; -0.103794571 to 0.212674466




##########
#Temperate
##########
MCMC_combined_CS_temperate <- combine_model_outputs(directory = "/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/CS_response/Outputs/Temperate/1st_run_11M_1M_10k/", file_pattern = "\\.rds")
length(MCMC_combined_CS_temperate[,1]) #Should be 400000

#Summary outputs
effectiveSize(MCMC_combined_CS_temperate)
posterior.mode(MCMC_combined_CS_temperate) #3.455701022; -0.009198375; -0.006663124; -0.016604765

#Compute 95% credible interval
HPDinterval(MCMC_combined_CS_temperate) #Only PRE significant: 2.27034509 to 4.517006365; -0.03838252 to 0.021080642; -0.01121559 to -0.002254338; -0.07910841 to 0.048575952




##########
##Species existing in both tropical and temperate regions (Both)
##########
MCMC_combined_CS_both <- combine_model_outputs(directory = "/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/CS_response/Outputs/Both/1st_run_11M_1M_10k/", file_pattern = "\\.rds")
length(MCMC_combined_CS_both[,1]) #Should be 400000

#Summary outputs
effectiveSize(MCMC_combined_CS_both)
posterior.mode(MCMC_combined_CS_both) #3.757708779; 0.034748021; -0.006377403; -0.082747549

#Compute 95% credible interval
HPDinterval(MCMC_combined_CS_both) #Only PRE significant: 1.12117194 to 6.2567530316; -0.01390962 to 0.0920878895; -0.01279867 to -0.0001097685; -0.25392248 to 0.0683470469





