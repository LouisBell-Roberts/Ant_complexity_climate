########## R2 calculations for colony size using MCMCglmm - Shinichi method ###########
####Across 400 trees
###Continuous data
#29/10/2023
#Louis Bell-Roberts


#Loading essential packages ####
library(data.table)
library(ape)
library(phytools)
library(MCMCglmm)
library(tidyverse)
library(gtools)
library(arm)
library(purrr)
library(gtools)
#

#########################################################

##
#####Custom functions
##


###
#
#Write a function that reads in all of the individual models
model_reader_function <- function() {
  
  ## get list of all rda files in the working directory
  model_files <- list.files(pattern = ".rds")
  ##sort the file names alphanumerically
  sorted_files <- mixedsort(model_files)
  ## read in all models
  # Create an empty list to store the loaded data
  mcmc_list <- list()
  
  #Read in each MCMCglmm model
  mcmc_list <- lapply(sorted_files, function(file) {
    readRDS(file)
  }) #mcmc_list is the list of all of the models read in
  return(mcmc_list)
}



###
#
#Function for calculating R2_marginal along with associated CI for gaussian response variables. Works on models with 3 fixed effects and 1 random effect (=phylogeny)
MCMCglmm_R2_marginal_with_CI_gaussian <- function(mcmc_object) {
  mmF_list <- split(mcmc_object$Sol, seq(nrow(mcmc_object$Sol)))
  vmVarF <- purrr::map(mmF_list, ~.x %*% t(mcmc_object$X)) %>% map_dbl(~var(as.vector(.x)))
  
  # getting R2_marginal
  R2_marginal <- vmVarF / (vmVarF + mcmc_object$VCV[, 1] + mcmc_object$VCV[, 2])
  result <- list(
    mean = mean(R2_marginal),
    posterior_mode = posterior.mode(R2_marginal),
    credible_interval = HPDinterval(R2_marginal)
  )
  return(result)
}



# ###
# #
# #Function for calculating R2_marginal along with associated CI for binary response variables. Works on models with 3 fixed effects and 1 random effect (=phylogeny)
# 
# MCMCglmm_R2_marginal_with_CI_binary <- function(mcmc_object) {
#   c2 <- (16 * sqrt(3)/(15 * pi))^2
#   mcmc_object$Sol1 <- mcmc_object$Sol/sqrt(1+c2 * mcmc_object$VCV[, "units"]) #Create new Sol variable
#   mcmc_object$VCV1 <- mcmc_object$VCV/(1+c2 * mcmc_object$VCV[, "units"]) #Create new VCV variable
#   
#   mmF_list <- split(mcmc_object$Sol1, seq(nrow(mcmc_object$Sol1)))
#   vmVarF <- purrr::map(mmF_list, ~.x %*% t(mcmc_object$X)) %>% map_dbl(~var(as.vector(.x)))
#   
#   # getting R2_marginal
#   R2_marginal <- vmVarF / (vmVarF + mcmc_object$VCV1[, 1] + pi^2/3) #VCV for residual intentionally not used (replaced with pi^2/3)
#   result <- list(
#     mean = mean(R2_marginal),
#     posterior_mode = posterior.mode(R2_marginal),
#     credible_interval = HPDinterval(R2_marginal)
#   )
#   return(result)
# }



###
#
#Function for calculating R2_conditional along with associated CI for Gaussian response variables. Works on models with 3 fixed effects and 1 random effect (=phylogeny)
MCMCglmm_R2_conditional_with_CI_gaussian <- function(mcmc_object) {
  mmF_list <- split(mcmc_object$Sol, seq(nrow(mcmc_object$Sol)))
  vmVarF <- purrr::map(mmF_list, ~.x %*% t(mcmc_object$X)) %>% map_dbl(~var(as.vector(.x)))
  
  # getting R2_conditional
  R2_conditional <- (vmVarF + mcmc_object$VCV[, 1]) / (vmVarF + mcmc_object$VCV[, 1] + mcmc_object$VCV[, 2])
  result <- list(
    mean = mean(R2_conditional),
    posterior_mode = posterior.mode(R2_conditional),
    credible_interval = HPDinterval(R2_conditional)
  )
  return(result)
}



# ###
# #
# #Function for calculating R2_conditional along with associated CI for binary response variables. Works on models with 3 fixed effects and 1 random effect (=phylogeny)
# MCMCglmm_R2_conditional_with_CI_binary <- function(mcmc_object) {
#   c2 <- (16 * sqrt(3)/(15 * pi))^2
#   mcmc_object$Sol1 <- mcmc_object$Sol/sqrt(1+c2 * mcmc_object$VCV[, "units"]) #Create new Sol variable
#   mcmc_object$VCV1 <- mcmc_object$VCV/(1+c2 * mcmc_object$VCV[, "units"]) #Create new VCV variable
#   
#   mmF_list <- split(mcmc_object$Sol1, seq(nrow(mcmc_object$Sol1)))
#   vmVarF <- purrr::map(mmF_list, ~.x %*% t(mcmc_object$X)) %>% map_dbl(~var(as.vector(.x)))
#   
#   # getting R2_conditional
#   R2_conditional <- (vmVarF + mcmc_object$VCV1[, 1]) / (vmVarF + mcmc_object$VCV1[, 1] + pi^2/3)
#   result <- list(
#     mean = mean(R2_conditional),
#     posterior_mode = posterior.mode(R2_conditional),
#     credible_interval = HPDinterval(R2_conditional)
#   )
#   return(result)
# }




###
#
#Function that converts marginal or conditional R2 list result to a dataframe
convert_results_list_to_dataframe <- function(results_list) {
  # Use lapply to extract values from the results_list
  mean_R2 <- unlist(lapply(results_list, function(x) x$mean))
  post_mode_R2 <- unlist(lapply(results_list, function(x) x$posterior_mode))
  lower_CI <- unlist(lapply(results_list, function(x) x$credible_interval[1]))
  upper_CI <- unlist(lapply(results_list, function(x) x$credible_interval[2]))
  
  # Create a data frame from the extracted values
  R2_distribution <- data.frame(mean_R2, post_mode_R2, lower_CI, upper_CI)
  
  return(R2_distribution)
}

#########################################################




#
#########################################################
#Colony size
#########################################################
#




#

###########
#Global
###########

#

#Set working directory
setwd('/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/CS_response/Outputs/Global/1st_run_11M_1M_10k/')

#Read in all of the individual models
mcmc_list_CS_1st_global <- model_reader_function()
#


####
#Marginal
####

#Calculate R2_marginal values across the 400 trees along with their associated CI
R2_marginal_single_global <- MCMCglmm_R2_marginal_with_CI_gaussian(mcmc_list_CS_1st_global[[1]]) #Apply to a single MCMCglmm model
R2_marginal_results_list_global <- lapply(mcmc_list_CS_1st_global, MCMCglmm_R2_marginal_with_CI_gaussian) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_marginal_distribution_global <- convert_results_list_to_dataframe(R2_marginal_results_list_global)
quantile(R2_marginal_distribution_global$post_mode_R2) #0.007999837


####
#Conditional
####

#Calculate R2_conditional values across the 400 trees along with their associated CI
R2_conditional_single_globa <- MCMCglmm_R2_conditional_with_CI_gaussian(mcmc_list_CS_1st_global[[1]]) #Apply to a single MCMCglmm model
R2_conditional_results_list_global <- lapply(mcmc_list_CS_1st_global, MCMCglmm_R2_conditional_with_CI_gaussian) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_conditional_distribution_global <- convert_results_list_to_dataframe(R2_conditional_results_list_global)
quantile(R2_conditional_distribution_global$post_mode_R2) #0.7257148



#

###########
#Tropical
###########

#

#Set working directory
setwd('/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/CS_response/Outputs/Tropical/1st_run_11M_1M_10k/')

#Read in all of the individual models
mcmc_list_CS_1st_tropical <- model_reader_function()
#


####
#Marginal
####

#Calculate R2_marginal values across the 400 trees along with their associated CI
R2_marginal_single_tropical <- MCMCglmm_R2_marginal_with_CI_gaussian(mcmc_list_CS_1st_tropical[[1]]) #Apply to a single MCMCglmm model
R2_marginal_results_list_tropical <- lapply(mcmc_list_CS_1st_tropical, MCMCglmm_R2_marginal_with_CI_gaussian) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_marginal_distribution_tropical <- convert_results_list_to_dataframe(R2_marginal_results_list_tropical)
quantile(R2_marginal_distribution_tropical$post_mode_R2) #0.008515972


####
#Conditional
####

#Calculate R2_conditional values across the 400 trees along with their associated CI
R2_conditional_single_tropical <- MCMCglmm_R2_conditional_with_CI_gaussian(mcmc_list_CS_1st_tropical[[1]]) #Apply to a single MCMCglmm model
R2_conditional_results_list_tropical <- lapply(mcmc_list_CS_1st_tropical, MCMCglmm_R2_conditional_with_CI_gaussian) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_conditional_distribution_tropical <- convert_results_list_to_dataframe(R2_conditional_results_list_tropical)
quantile(R2_conditional_distribution_tropical$post_mode_R2) #0.8821487


#

###########
#Temperate
###########

#

#Set working directory
setwd('/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/CS_response/Outputs/Temperate/1st_run_11M_1M_10k/')

#Read in all of the individual models
mcmc_list_CS_1st_temperate <- model_reader_function()
#


####
#Marginal
####

#Calculate R2_marginal values across the 400 trees along with their associated CI
R2_marginal_single_temperate <- MCMCglmm_R2_marginal_with_CI_gaussian(mcmc_list_CS_1st_temperate[[1]]) #Apply to a single MCMCglmm model
R2_marginal_results_list_temperate <- lapply(mcmc_list_CS_1st_temperate, MCMCglmm_R2_marginal_with_CI_gaussian) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_marginal_distribution_temperate <- convert_results_list_to_dataframe(R2_marginal_results_list_temperate)
quantile(R2_marginal_distribution_temperate$post_mode_R2) #0.04081920


####
#Conditional
####

#Calculate R2_conditional values across the 400 trees along with their associated CI
R2_conditional_single_temperate <- MCMCglmm_R2_conditional_with_CI_gaussian(mcmc_list_CS_1st_temperate[[1]]) #Apply to a single MCMCglmm model
R2_conditional_results_list_temperate <- lapply(mcmc_list_CS_1st_temperate, MCMCglmm_R2_conditional_with_CI_gaussian) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_conditional_distribution_temperate <- convert_results_list_to_dataframe(R2_conditional_results_list_temperate)
quantile(R2_conditional_distribution_temperate$post_mode_R2) #0.7586201


#

###########
#Both
###########

#

#Set working directory
setwd('/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/CS_response/Outputs/Both/1st_run_11M_1M_10k/')

#Read in all of the individual models
mcmc_list_CS_1st_both <- model_reader_function()
#


####
#Marginal
####

#Calculate R2_marginal values across the 400 trees along with their associated CI
R2_marginal_single_both <- MCMCglmm_R2_marginal_with_CI_gaussian(mcmc_list_CS_1st_both[[1]]) #Apply to a single MCMCglmm model
R2_marginal_results_list_both <- lapply(mcmc_list_CS_1st_both, MCMCglmm_R2_marginal_with_CI_gaussian) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_marginal_distribution_both <- convert_results_list_to_dataframe(R2_marginal_results_list_both)
quantile(R2_marginal_distribution_both$post_mode_R2) #0.01865943


####
#Conditional
####

#Calculate R2_conditional values across the 400 trees along with their associated CI
R2_conditional_single_both <- MCMCglmm_R2_conditional_with_CI_gaussian(mcmc_list_CS_1st_both[[1]]) #Apply to a single MCMCglmm model
R2_conditional_results_list_both <- lapply(mcmc_list_CS_1st_both, MCMCglmm_R2_conditional_with_CI_gaussian) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_conditional_distribution_both <- convert_results_list_to_dataframe(R2_conditional_results_list_both)
quantile(R2_conditional_distribution_both$post_mode_R2) #0.7577073












