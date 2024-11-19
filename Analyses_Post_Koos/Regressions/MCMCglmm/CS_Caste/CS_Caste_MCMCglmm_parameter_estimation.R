############################### BPMMs: Parameter estimation for MCMCglmm analyses of variation in worker size (monomorphic/polymorphic) and colony size ########################################
####To quantifying phylogenetic uncertainty in parameter estimates, analyses run over a sample of 400 phylogenetic trees
###Written for West server
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
MCMC_combined_Caste_global <- combine_model_outputs(directory = "/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS_Caste/Outputs/Global/1st_run/", file_pattern = "\\.rds")
length(MCMC_combined_Caste_global[,1]) #Should be 400000

#Summary outputs
effectiveSize(MCMC_combined_Caste_global)
posterior.mode(MCMC_combined_Caste_global) #-7.441706; 1.843471 

#Compute 95% credible interval
HPDinterval(MCMC_combined_Caste_global) #-10.628002 -4.859915; 1.327901  2.472161



##########
#Tropical
##########
#When using the combine_model_outputs function - choose the directory folder containing MCMCglmm objects
MCMC_combined_Caste_tropical <- combine_model_outputs(directory = "/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS_Caste/Outputs/Tropical/1st_run/", file_pattern = "\\.rds")
length(MCMC_combined_Caste_tropical[,1]) #Should be 400000

#Summary outputs
effectiveSize(MCMC_combined_Caste_tropical)
posterior.mode(MCMC_combined_Caste_tropical) #-5.997881; 1.382855

#Compute 95% credible interval
HPDinterval(MCMC_combined_Caste_tropical) #-9.6470026 -3.411486; 0.7356547  2.290619




##########
#Temperate
##########
MCMC_combined_Caste_temperate <- combine_model_outputs(directory = "/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS_Caste/Outputs/Temperate/1st_run/", file_pattern = "\\.rds")
length(MCMC_combined_Caste_temperate[,1]) #Should be 400000

#Summary outputs
effectiveSize(MCMC_combined_Caste_temperate)
posterior.mode(MCMC_combined_Caste_temperate) #-8.315201; 1.837393

#Compute 95% credible interval
HPDinterval(MCMC_combined_Caste_temperate) #-12.514578 -4.670468; 1.051748  2.902802




##########
##Species existing in both tropical and temperate regions (Both)
##########
MCMC_combined_Caste_both <- combine_model_outputs(directory = "/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS_Caste/Outputs/Both/1st_run/", file_pattern = "\\.rds")
length(MCMC_combined_Caste_both[,1]) #Should be 400000

#Summary outputs
effectiveSize(MCMC_combined_Caste_both)
posterior.mode(MCMC_combined_Caste_both) #-5.673313; 1.670599

#Compute 95% credible interval
HPDinterval(MCMC_combined_Caste_both) #-8.861151 -3.566728; 1.057250  2.453883





