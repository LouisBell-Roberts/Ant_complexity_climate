######################## AntEnv: BPMM analysis diagnostic checks over 400 models - Both dataset ######################
###For both colony size and variation in worker size response variables
##Run diagnostic checks: ESS, Gelman-Rubin scale reduction factor test (SRF), autocorrelation
#Louis Bell-Roberts
#25/10/2023


library(tidyverse)
library(MCMCglmm)
library(ggplot2)
library(mulTree)
library(gtools)

##
#####Custom functions
##

###
#Calculates SRF (Gelman-Rubin) estimate
srf_function <- function(multree_conv) {
  
  # Initialize an empty vector to store the SRFs
  list_psrf <- c()
  
  # Loop through the list of gelman.diag objects
  for (i in 1:length(multree_conv)) {
    # Get the first column of the ith object and calculate its range
    ith_psrf <- multree_conv[[i]]$psrf #Need to look at the psrf for the upper CI
    
    list_psrf <- c(list_psrf, ith_psrf)
  }
  #PSR statistics
  return(list_psrf)
}


###
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
#Write a function that tests levels of autocorrelation for the fixed effects for each model - provides the values in the second row of the autocorr.diag function. These should be <0.1
auto_cor_function_fix <- function(mcmc_list){
  
  # Apply the autocorr.diag function to each element of the list
  autocorr_list <- lapply(mcmc_list, function(x) {
    autocorr.diag(x$Sol)
  })
  
  # Use lapply to extract the second row of each matrix
  second_row_list <- lapply(autocorr_list, function(x) {
    x[2,]
  })
  
  # Convert the list to a vector
  second_row_vec <- unlist(second_row_list)
  return(second_row_vec)
}



###
#Write a function that tests levels of autocorrelation for the random effects for each model - provides the values in the second row of the autocorr.diag function. These should be <0.1
auto_cor_function_rand <- function(mcmc_list){
  
  # Apply the autocorr.diag function to each element of the list
  autocorr_list <- lapply(mcmc_list, function(x) {
    autocorr.diag(x$VCV)
  })
  
  # Use lapply to extract the second row of each matrix
  second_row_list <- lapply(autocorr_list, function(x) {
    x[2,]
  })
  
  # Convert the list to a vector
  second_row_vec <- unlist(second_row_list)
  return(second_row_vec)
}


###
#Two functions: one that calculates the ESS values for the fixed effects for each model across the full sample of models and one that calculates the ESS values for the random effects for each model across the full sample of models

##Fixed
ESS_function_fix <- function(mcmc_list) {
  
  ESS_list <- lapply(mcmc_list, function(x) {
    effectiveSize(x$Sol)
  })
  
  # Convert the list to a vector
  ESS_vec <- unlist(ESS_list)
  
  return(ESS_vec)
}
##Random
ESS_function_rand <- function(mcmc_list) {
  
  ESS_list <- lapply(mcmc_list, function(x) {
    effectiveSize(x$VCV)
  })
  
  # Convert the list to a vector
  ESS_vec <- unlist(ESS_list)
  
  return(ESS_vec)
}




##################

####################################################
#Colony size
####################################################

##################

#Set working directory
setwd("/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/CS_response/Outputs/Both/1st_run_11M_1M_10k/")

#Read in all of the individual models
mcmc_list_CS_1st <- model_reader_function()

#Set working directory again
setwd("/Volumes/ADATA SE800/DOL_worker_castes/AntEnv/Regressions/MCMCglmm/CS_response/Outputs/Both/2nd_run_11M_1M_10k/")

mcmc_list_CS_2nd <- model_reader_function()


#######
#Calculate whether scale-reduction factor (Gelman-Rubin) is <1.1 for all 400 models
#######

# Create a list to store the paired lists
paired_lists_CS_Sol <- list()

# Iterate over the indices of the lists
for (i in 1:length(mcmc_list_CS_1st)) {
  # Select the "$Sol" element from each object in list_a and list_b
  CS_pair_Sol <- list(mcmc_list_CS_1st[[i]]$Sol, mcmc_list_CS_2nd[[i]]$Sol)
  
  # Add the pair to the list of paired lists
  paired_lists_CS_Sol[[i]] <- CS_pair_Sol
}

#Run gelman.diag function over the list of paired chain runs
gelman_results_CS_Sol <- lapply(paired_lists_CS_Sol, gelman.diag)

#Quantify whether each Sol passes srf test
##Create a single large data frame which holds all of the srf scores, and can then check whether they exceed 1.1

# Initialize an empty data frame
CS_Sol_df <- data.frame()
# Combine matrices into a single data frame
for (i in 1:length(gelman_results_CS_Sol)) {
  # Convert the matrix to a data frame
  CS_Sol_temporary_df <- as.data.frame(gelman_results_CS_Sol[[i]]$psrf)
  
  # Add an identifier column to keep track of the original matrix
  CS_Sol_temporary_df$matrix_id <- i
  
  # Append the temporary data frame to the combined data frame
  CS_Sol_df <- rbind(CS_Sol_df, CS_Sol_temporary_df)
} #Passes test
#All Sols <1.1 NEED TO LOOK AT UPPER CI - all trees except one have psrf <1.1
#Within the data frame produced, use 'matrix_id' to identify which tree the PSRF score came from 

#Summary of results
quantile(unlist(gelman_results_CS_Sol)) #All values <1.1
gelman_results_CS_Sol[[103]] #View for a single model


#######
##Repeat process for the VCV objects
#######

#Create a list of paired mcmc objects

paired_lists_CS_VCV <- list()

# Iterate over the indices of the lists
for (i in 1:length(mcmc_list_CS_1st)) {
  # Select the "$Sol" element from each object in list_a and list_b
  CS_pair_VCV <- list(mcmc_list_CS_1st[[i]]$VCV, mcmc_list_CS_2nd[[i]]$VCV)
  
  # Add the pair to the list of paired lists
  paired_lists_CS_VCV[[i]] <- CS_pair_VCV
}

#Calculate srf scores for each mcmc object pair in the list
CS_gelman_results_VCV <- lapply(paired_lists_CS_VCV, gelman.diag)

#Quantify whether each VCV passes srf test
##Create a single large data frame which holds all of the srf scores, and can then check whether they exceed 1.1

# Initialize an empty data frame
CS_VCV_df <- data.frame()
# Combine matrices into a single data frame
for (i in 1:length(CS_gelman_results_VCV)) {
  # Convert the matrix to a data frame
  CS_VCV_temporary_df <- as.data.frame(CS_gelman_results_VCV[[i]]$psrf)
  
  # Add an identifier column to keep track of the original matrix
  CS_VCV_temporary_df$matrix_id <- i
  
  # Append the temporary data frame to the combined data frame
  CS_VCV_df <- rbind(CS_VCV_df, CS_VCV_temporary_df)
} #All VCVs <1.1 NEED TO LOOK AT UPPER CI
# Good enough - only one value over 1.1 (1.118662)

#Summary of results
quantile(unlist(CS_gelman_results_VCV)) #All values <1.1
CS_gelman_results_VCV[[103]] #View for a single model




########
#Test levels of autocorrelation for the fixed and random effects for each model 
########

##Fixed effects

# Apply the autocorr.diag function to each element of the list
autocorr_list_CS_fix <- lapply(mcmc_list_CS_1st, function(x) {
  autocorr.diag(x$Sol)
})

# Use lapply to extract the second row of each matrix
second_row_list_fix_CS <- lapply(autocorr_list_CS_fix, function(x) {
  x[2,]
})

# Convert the list to a vector
second_row_vec_fix_CS <- unlist(second_row_list_fix_CS)
quantile(second_row_vec_fix_CS) # -0.0968987951 to 0.0974697560; All values <0.1
# View(as.data.frame(second_row_vec_fix_CS))

##Random effects

# Apply the autocorr.diag function to each element of the list
autocorr_list_CS_rand <- lapply(mcmc_list_CS_1st, function(x) {
  autocorr.diag(x$VCV)
})

# Use lapply to extract the second row of each matrix
second_row_list_rand_CS <- lapply(autocorr_list_CS_rand, function(x) {
  x[2,]
})

# Convert the list to a vector
second_row_vec_rand_CS <- unlist(second_row_list_rand_CS)
quantile(second_row_vec_rand_CS) # -0.094448644 to 0.099429489; All values <0.1
# View(as.data.frame(second_row_vec_rand_CS))





#######
#Calculate effective sample size for both the fixed and random effects for each model
#######

##Fixed effects

#Test ESS levels for the fixed effects for each model 
# Apply the autocorr.diag function to each element of the list
ESS_list_CS_fix <- lapply(mcmc_list_CS_1st, function(x) {
  effectiveSize(x$Sol)
})

# Convert the list to a vector
ESS_vec_CS_fix <- unlist(ESS_list_CS_fix)
quantile(ESS_vec_CS_fix) #lowest ESS = 492.1955; generally values are around 1000
# View(as.data.frame(ESS_vec_CS_fix))


##Random effects

#Test ESS levels for the random effects for each model 
# Apply the autocorr.diag function to each element of the list
ESS_list_CS_rand <- lapply(mcmc_list_CS_1st, function(x) {
  effectiveSize(x$VCV)
})

# Convert the list to a vector
ESS_vec_CS_rand <- unlist(ESS_list_CS_rand)
quantile(ESS_vec_CS_rand) #lowest ESS = 609.8032; generally values are around 1000
# View(as.data.frame(ESS_vec_CS_rand))

##################

####################################################
#Variation in worker size
####################################################

##################