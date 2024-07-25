######################## AntEnv: BPMM analysis diagnostic checks over 400 models ######################
###For both colony size and variation in worker size response variables
##Run diagnostic checks: ESS, Gelman-Rubin scale reduction factor test (SRF), autocorrelation
#Louis Bell-Roberts
#25/10/2023


library(tidyverse)
library(MCMCglmm)
library(ggplot2)
library(mulTree)

##
#####Custom functions
##


#Calculates SRF (Gelman-Rubin) estimate
##Input the convergence models from read.mulTree with convergence set to TRUE
# srf_function <- function(multree_conv) {
#   
#   # Initialize an empty vector to store the SRFs
#   list_psrf <- c()
#   
#   # Loop through the list of gelman.diag objects
#   for (i in 1:length(multree_conv)) {
#     # Get the first column of the ith object and calculate its range
#     ith_psrf <- multree_conv[[i]]$psrf[,2]
#     
#     list_psrf <- c(list_psrf, ith_psrf)
#   }
#   
#   #PSR statistics
#   return(quantile(list_psrf))
#   # View(as.data.frame(list_psrf_CS))
#   return(length(list_psrf))
# }

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





#Write a function that reads in all of the individual models
model_reader_function <- function() {
  
  ## get list of all rda files in the working directory
  model_files <- list.files(pattern = ".rds")
  ## read in all models
  # Create an empty list to store the loaded data
  mcmc_list <- list()
  
  #Read in each MCMCglmm model
  mcmc_list <- lapply(model_files, function(file) {
    readRDS(file)
  }) #mcmc_list is the list of all of the models read in
  return(mcmc_list)
}






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
#CS
####################################################

##################




#############
#Global
#############

#Set working directory
setwd("/...")

##
#Read in all of the individual models
##
##MAKE SURE WORKING DIRECTORY HAS BEEN SET CORRECTLY ALREADY
mcmc_list_CS_global_1st <- model_reader_function() #Are these read in in the correct order?

#Set working directory again
setwd("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_signal/Outputs/CS/1M/2nd_run/")

mcmc_list_CS_global_2nd <- model_reader_function()


#######
#Calculate whether scale-reduction factor (Gelman-Rubin) is <1.1 for all 400 models
#######
#For srf function to work, it needs:
##gelman.diag() function
###gelman.diag function takes lists of mcmc objects - therefore, need to pair the different chain runs (this will give me 400 pairs of MCMCglmm objects)
####the mcmc objects can be either: 1. the Sols (fixed effects and nodes/animals), 2. the VCVs
#####Write code to calculate srf for each

#Create a list of paired mcmc objects

# Create a list to store the paired lists
paired_lists_CS_global_Sol <- list()

# Iterate over the indices of the lists
for (i in 1:length(mcmc_list_CS_global_1st)) {
  # Select the "$Sol" element from each object in list_a and list_b
  CS_global_pair_Sol <- list(mcmc_list_CS_global_1st[[i]]$Sol, mcmc_list_CS_global_2nd[[i]]$Sol)
  
  # Add the pair to the list of paired lists
  paired_lists_CS_global_Sol[[i]] <- CS_global_pair_Sol
}

#Run gelman function on a single pair from the list
# gelman.diag(paired_lists_CS_Sol[[1]])
# View(as.data.frame(unlist(gelman.diag(paired_lists_CS_Sol[[1]]))))

#Run gelman.diag function over the list of paired chain runs
gelman_results_CS_global_Sol <- lapply(paired_lists_CS_global_Sol, gelman.diag) #Takes a while to run

quantile(unlist(gelman_results_CS_global_Sol))
# View(as.data.frame(unlist(gelman_results_CS_Sol[[4]]))) #View for a single model
# View(as.data.frame(unlist(gelman_results_CS_Sol))) #View all together

#Quantify whether each Sol passes srf test
##Create a single large data frame which holds all of the srf scores, and can then check whether they exceed 1.1

# Initialize an empty data frame
CS_global_Sol_df <- data.frame()
# Combine matrices into a single data frame
for (i in 1:length(gelman_results_CS_global_Sol)) {
  # Convert the matrix to a data frame
  CS_global_Sol_temporary_df <- as.data.frame(gelman_results_CS_global_Sol[[i]]$psrf)
  
  # Add an identifier column to keep track of the original matrix
  CS_global_Sol_temporary_df$matrix_id <- i
  
  # Append the temporary data frame to the combined data frame
  CS_global_Sol_df <- rbind(CS_global_Sol_df, CS_global_Sol_temporary_df)
} #Passes test
#All Sols <1.1 NEED TO LOOK AT UPPER CI


#######
##Repeat process for the VCV objects
#######


#Create a list of paired mcmc objects

paired_lists_CS_global_VCV <- list()

# Iterate over the indices of the lists
for (i in 1:length(mcmc_list_CS_global_1st)) {
  # Select the "$Sol" element from each object in list_a and list_b
  CS_global_pair_VCV <- list(mcmc_list_CS_global_1st[[i]]$VCV, mcmc_list_CS_global_2nd[[i]]$VCV)
  
  # Add the pair to the list of paired lists
  paired_lists_CS_global_VCV[[i]] <- CS_global_pair_VCV
}

#Calculate srf scores for each mcmc object pair in the list
CS_global_gelman_results_VCV <- lapply(paired_lists_CS_global_VCV, gelman.diag)
quantile(unlist(CS_global_gelman_results_VCV))

#Quantify whether each VCV passes srf test
##Create a single large data frame which holds all of the srf scores, and can then check whether they exceed 1.1

# Initialize an empty data frame
CS_global_VCV_df <- data.frame()
# Combine matrices into a single data frame
for (i in 1:length(CS_global_gelman_results_VCV)) {
  # Convert the matrix to a data frame
  CS_global_VCV_temporary_df <- as.data.frame(CS_global_gelman_results_VCV[[i]]$psrf)
  
  # Add an identifier column to keep track of the original matrix
  CS_global_VCV_temporary_df$matrix_id <- i
  
  # Append the temporary data frame to the combined data frame
  CS_global_VCV_df <- rbind(CS_global_VCV_df, CS_global_VCV_temporary_df)
} #All VCVs <1.1 NEED TO LOOK AT UPPER CI
# Good enough - only one value over 1.1 (1.118662)


########
#Test levels of autocorrelation for the fixed and random effects for each model 
########



##Fixed effects

# Apply the autocorr.diag function to each element of the list
autocorr_list_CS_global_fix <- lapply(mcmc_list_CS_global_1st, function(x) {
  autocorr.diag(x$Sol)
})

# Use lapply to extract the second row of each matrix
second_row_list_fix_CS_global <- lapply(autocorr_list_CS_global_fix, function(x) {
  x[2,]
})

# Convert the list to a vector
second_row_vec_fix_CS_global <- unlist(second_row_list_fix_CS_global)
quantile(second_row_vec_fix_CS_global) # -0.105209364 to 0.098130070
# View(as.data.frame(second_row_vec_fix_CS_global))

##Random effects

# Apply the autocorr.diag function to each element of the list
autocorr_list_CS_global_rand <- lapply(mcmc_list_CS_global_1st, function(x) {
  autocorr.diag(x$VCV)
})

# Use lapply to extract the second row of each matrix
second_row_list_rand_CS_global <- lapply(autocorr_list_CS_global_rand, function(x) {
  x[2,]
})

# Convert the list to a vector
second_row_vec_rand_CS_global <- unlist(second_row_list_rand_CS_global)
quantile(second_row_vec_rand_CS_global) # -0.1284119764 to 0.0847181148
# View(as.data.frame(second_row_vec_rand_CS_global))

#######
#Calculate effective sample size for both the fixed and random effects for each model
#######

##Fixed effects

#Test levels of autocorrelation for the fixed effects for each model 
# Apply the autocorr.diag function to each element of the list
ESS_list_CS_global_fix <- lapply(mcmc_list_CS_global_1st, function(x) {
  effectiveSize(x$Sol)
})

# Convert the list to a vector
ESS_vec_CS_global_fix <- unlist(ESS_list_CS_global_fix)
quantile(ESS_vec_CS_global_fix) #lowest ESS = 566.224
# View(as.data.frame(ESS_vec_CS_global_fix))


##Random effects

#Test levels of autocorrelation for the random effects for each model 
# Apply the autocorr.diag function to each element of the list
ESS_list_CS_global_rand <- lapply(mcmc_list_CS_global_1st, function(x) {
  effectiveSize(x$VCV)
})

# Convert the list to a vector
ESS_vec_CS_global_rand <- unlist(ESS_list_CS_global_rand)
quantile(ESS_vec_CS_global_rand) #lowest ESS = 638.8511
# View(as.data.frame(ESS_vec_CS_global_rand))


###########
#Tropical
###########

#Set working directory
setwd("/...")

##
#Read in all of the individual models
##
##MAKE SURE WORKING DIRECTORY HAS BEEN SET CORRECTLY ALREADY
mcmc_list_CS_tropical_1st <- model_reader_function() #Are these read in in the correct order?

#Set working directory again
setwd("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_signal/Outputs/CS/1M/2nd_run/")

mcmc_list_CS_tropical_2nd <- model_reader_function()


#######
#Calculate whether scale-reduction factor (Gelman-Rubin) is <1.1 for all 400 models
#######
#For srf function to work, it needs:
##gelman.diag() function
###gelman.diag function takes lists of mcmc objects - therefore, need to pair the different chain runs (this will give me 400 pairs of MCMCglmm objects)
####the mcmc objects can be either: 1. the Sols (fixed effects and nodes/animals), 2. the VCVs
#####Write code to calculate srf for each

#Create a list of paired mcmc objects

# Create a list to store the paired lists
paired_lists_CS_tropical_Sol <- list()

# Iterate over the indices of the lists
for (i in 1:length(mcmc_list_CS_tropical_1st)) {
  # Select the "$Sol" element from each object in list_a and list_b
  CS_tropical_pair_Sol <- list(mcmc_list_CS_tropical_1st[[i]]$Sol, mcmc_list_CS_tropical_2nd[[i]]$Sol)
  
  # Add the pair to the list of paired lists
  paired_lists_CS_tropical_Sol[[i]] <- CS_tropical_pair_Sol
}

#Run gelman function on a single pair from the list
# gelman.diag(paired_lists_CS_Sol[[1]])
# View(as.data.frame(unlist(gelman.diag(paired_lists_CS_Sol[[1]]))))

#Run gelman.diag function over the list of paired chain runs
gelman_results_CS_tropical_Sol <- lapply(paired_lists_CS_tropical_Sol, gelman.diag) #Takes a while to run

quantile(unlist(gelman_results_CS_tropical_Sol))
# View(as.data.frame(unlist(gelman_results_CS_Sol[[4]]))) #View for a single model
# View(as.data.frame(unlist(gelman_results_CS_Sol))) #View all together

#Quantify whether each Sol passes srf test
##Create a single large data frame which holds all of the srf scores, and can then check whether they exceed 1.1

# Initialize an empty data frame
CS_tropical_Sol_df <- data.frame()
# Combine matrices into a single data frame
for (i in 1:length(gelman_results_CS_tropical_Sol)) {
  # Convert the matrix to a data frame
  CS_global_Sol_temporary_df <- as.data.frame(gelman_results_CS_tropical_Sol[[i]]$psrf)
  
  # Add an identifier column to keep track of the original matrix
  CS_global_Sol_temporary_df$matrix_id <- i
  
  # Append the temporary data frame to the combined data frame
  CS_global_Sol_df <- rbind(CS_global_Sol_df, CS_global_Sol_temporary_df)
} #Passes test
#All Sols <1.1 NEED TO LOOK AT UPPER CI


#######
##Repeat process for the VCV objects
#######


#Create a list of paired mcmc objects

paired_lists_CS_global_VCV <- list()

# Iterate over the indices of the lists
for (i in 1:length(mcmc_list_CS_global_1st)) {
  # Select the "$Sol" element from each object in list_a and list_b
  CS_global_pair_VCV <- list(mcmc_list_CS_global_1st[[i]]$VCV, mcmc_list_CS_global_2nd[[i]]$VCV)
  
  # Add the pair to the list of paired lists
  paired_lists_CS_global_VCV[[i]] <- CS_global_pair_VCV
}

#Calculate srf scores for each mcmc object pair in the list
CS_global_gelman_results_VCV <- lapply(paired_lists_CS_global_VCV, gelman.diag)
quantile(unlist(CS_global_gelman_results_VCV))

#Quantify whether each VCV passes srf test
##Create a single large data frame which holds all of the srf scores, and can then check whether they exceed 1.1

# Initialize an empty data frame
CS_global_VCV_df <- data.frame()
# Combine matrices into a single data frame
for (i in 1:length(CS_global_gelman_results_VCV)) {
  # Convert the matrix to a data frame
  CS_global_VCV_temporary_df <- as.data.frame(CS_global_gelman_results_VCV[[i]]$psrf)
  
  # Add an identifier column to keep track of the original matrix
  CS_global_VCV_temporary_df$matrix_id <- i
  
  # Append the temporary data frame to the combined data frame
  CS_global_VCV_df <- rbind(CS_global_VCV_df, CS_global_VCV_temporary_df)
} #All VCVs <1.1 NEED TO LOOK AT UPPER CI
# Good enough - only one value over 1.1 (1.118662)


########
#Test levels of autocorrelation for the fixed and random effects for each model 
########



##Fixed effects

# Apply the autocorr.diag function to each element of the list
autocorr_list_CS_global_fix <- lapply(mcmc_list_CS_global_1st, function(x) {
  autocorr.diag(x$Sol)
})

# Use lapply to extract the second row of each matrix
second_row_list_fix_CS_global <- lapply(autocorr_list_CS_global_fix, function(x) {
  x[2,]
})

# Convert the list to a vector
second_row_vec_fix_CS_global <- unlist(second_row_list_fix_CS_global)
quantile(second_row_vec_fix_CS_global) # -0.105209364 to 0.098130070
# View(as.data.frame(second_row_vec_fix_CS_global))

##Random effects

# Apply the autocorr.diag function to each element of the list
autocorr_list_CS_global_rand <- lapply(mcmc_list_CS_global_1st, function(x) {
  autocorr.diag(x$VCV)
})

# Use lapply to extract the second row of each matrix
second_row_list_rand_CS_global <- lapply(autocorr_list_CS_global_rand, function(x) {
  x[2,]
})

# Convert the list to a vector
second_row_vec_rand_CS_global <- unlist(second_row_list_rand_CS_global)
quantile(second_row_vec_rand_CS_global) # -0.1284119764 to 0.0847181148
# View(as.data.frame(second_row_vec_rand_CS_global))

#######
#Calculate effective sample size for both the fixed and random effects for each model
#######

##Fixed effects

#Test levels of autocorrelation for the fixed effects for each model 
# Apply the autocorr.diag function to each element of the list
ESS_list_CS_global_fix <- lapply(mcmc_list_CS_global_1st, function(x) {
  effectiveSize(x$Sol)
})

# Convert the list to a vector
ESS_vec_CS_global_fix <- unlist(ESS_list_CS_global_fix)
quantile(ESS_vec_CS_global_fix) #lowest ESS = 566.224
# View(as.data.frame(ESS_vec_CS_global_fix))


##Random effects

#Test levels of autocorrelation for the random effects for each model 
# Apply the autocorr.diag function to each element of the list
ESS_list_CS_global_rand <- lapply(mcmc_list_CS_global_1st, function(x) {
  effectiveSize(x$VCV)
})

# Convert the list to a vector
ESS_vec_CS_global_rand <- unlist(ESS_list_CS_global_rand)
quantile(ESS_vec_CS_global_rand) #lowest ESS = 638.8511
# View(as.data.frame(ESS_vec_CS_global_rand))

