############################### BPMMs: Parameter estimation for MCMCglmm analyses of variation in worker size and climatic variables (TMP PRE DTR) ########################################
####To quantifying phylogenetic uncertainty in parameter estimates, analyses run over a sample of 400 phylogenetic trees
###Written for Louis' laptop
##Louis Bell-Roberts
#29/10/2023


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
MCMC_combined_Caste_global <- combine_model_outputs(directory = "/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/Caste_response/Outputs/Global/No_fixed_effect_X2_rand/1st_run_11M_1M_10k/", file_pattern = "\\.rds")
length(MCMC_combined_Caste_global[,1]) #Should be 400000

#Summary outputs
effectiveSize(MCMC_combined_Caste_global)
posterior.mode(MCMC_combined_Caste_global) #-1.09297134; 0.16091680; -0.01633298; -0.22558481

#Compute 95% credible interval
HPDinterval(MCMC_combined_Caste_global) #TMP and PRE significant: -5.92738176 to 4.167174335; 0.04989480 to 0.284414644; -0.02923016 to -0.005347469; -0.56165504 to 0.031634083



##########
#Tropical
##########
#When using the combine_model_outputs function - choose the directory folder containing MCMCglmm objects
MCMC_combined_Caste_tropical <- combine_model_outputs(directory = "/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/Caste_response/Outputs/Tropical/No_fixed_effect_X2_rand/1st_run_11M_1M_10k/", file_pattern = "\\.rds")
length(MCMC_combined_Caste_tropical[,1]) #Should be 400000

#Summary outputs
effectiveSize(MCMC_combined_Caste_tropical)
posterior.mode(MCMC_combined_Caste_tropical) #-7.979796502; 0.372233624; -0.003548668; -0.308833657

#Compute 95% credible interval
HPDinterval(MCMC_combined_Caste_tropical) #None significant: -31.51510275 to 14.436683015; -0.33850812 to 1.191552594; -0.01915795 to 0.008875661; -1.05024785 to 0.270761176




##########
#Temperate
##########
MCMC_combined_Caste_temperate <- combine_model_outputs(directory = "/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/Caste_response/Outputs/Temperate/No_fixed_effect_X2_rand/1st_run_11M_1M_10k/", file_pattern = "\\.rds")
length(MCMC_combined_Caste_temperate[,1]) #Should be 400000

#Summary outputs
effectiveSize(MCMC_combined_Caste_temperate)
posterior.mode(MCMC_combined_Caste_temperate) #-1.62089575; 0.14247321; -0.01864425; -0.14820654

#Compute 95% credible interval
HPDinterval(MCMC_combined_Caste_temperate) #None significant: -7.71055751 to 4.556643008; -0.01474935 to 0.323105801; -0.04829187 to 0.003919135; -0.57944455 to 0.214010752




##########
##Species existing in both tropical and temperate regions (Both)
##########
MCMC_combined_Caste_both <- combine_model_outputs(directory = "/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/Caste_response/Outputs/Both/No_fixed_effect_X2_rand/1st_run_11M_1M_10k/", file_pattern = "\\.rds")
length(MCMC_combined_Caste_both[,1]) #Should be 400000

#Summary outputs
effectiveSize(MCMC_combined_Caste_both)
posterior.mode(MCMC_combined_Caste_both) #-2.49458532; 0.28510567; -0.02997316; -0.07837542

#Compute 95% credible interval
HPDinterval(MCMC_combined_Caste_both) #TMP and PRE significant: -13.36327173 to 6.665367993; 0.08645464 to 0.539172942; -0.05563978 to -0.005704148; -0.68850595  0.531275125





