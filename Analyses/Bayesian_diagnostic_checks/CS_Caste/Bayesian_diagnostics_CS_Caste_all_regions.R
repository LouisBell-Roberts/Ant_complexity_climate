######################## AntEnv: BPMM analysis diagnostic checks over 400 models ######################
###For analysis of the relationship between colony size and worker polymorphism (monomorphic/polymorphic) across all 4 geographic scales
##Run diagnostic checks: ESS, Gelman-Rubin scale reduction factor test (SRF), autocorrelation
#Louis Bell-Roberts
#20/03/2024

.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))

#Models pass all diagnostic tests. ESS is >500 for all models


library(tidyverse)
library(MCMCglmm)
library(ggplot2)
library(gtools)

############
#Custom functions
############

###
#Write a function that reads in all of the individual models
model_reader_function <- function(region, run) {
  # Construct directory path
  directory_path <- file.path("/drives/4tb/Louis/AntEnv/Analyses/Regressions/MCMCglmm/CS_Caste/Outputs", region, run)
  
  # Get list of all rda files in the specified directory
  model_files <- list.files(path = directory_path, pattern = ".rds")
  
  # Sort the file names alphanumerically
  sorted_files <- mixedsort(model_files)
  
  # Read in all models
  # Create an empty list to store the loaded data
  mcmc_list <- list()
  
  # Read in each MCMCglmm model
  mcmc_list <- lapply(sorted_files, function(file) {
    readRDS(file.path(directory_path, file))
  })  # mcmc_list is the list of all the models read in
  
  return(mcmc_list)
}



###
#Convergence test function
gelman_rubin_function <- function(mcmc_list_1, mcmc_list_2, component = "Sol") {
  # Create a list to store the paired lists
  paired_lists <- list()
  
  # Iterate over the indices of the lists
  for (i in seq_along(mcmc_list_1)) {
    if (component == "Sol") {
      # Select the "$Sol" element from each object in list_a and list_b
      pair_component <- list(mcmc_list_1[[i]]$Sol, mcmc_list_2[[i]]$Sol)
    } else if (component == "VCV") {
      # Select the "$VCV" element from each object in list_a and list_b
      pair_component <- list(mcmc_list_1[[i]]$VCV[, 1], mcmc_list_2[[i]]$VCV[, 1]) #Needed to select the first column (the column for the additive genetic variance), rather than select both columns, as the second column is for the residual variance which is fixed at 1. Using values for the residual variance, which are all fixed at one, causes errors down the line, plus we do not need to check their convergence properties.
    } else {
      stop("Invalid component argument. Use 'Sol' or 'VCV'.")
    }
    
    # Add the pair to the list of paired lists
    paired_lists[[i]] <- pair_component
  }
  
  # Run gelman.diag function over the list of paired chain runs
  gelman_results <- lapply(paired_lists, gelman.diag)
  
  # Initialize an empty data frame
  component_df <- data.frame()
  
  # Combine matrices into a single data frame
  for (i in seq_along(gelman_results)) {
    # Convert the matrix to a data frame
    component_temporary_df <- as.data.frame(gelman_results[[i]]$psrf)
    
    # Add an identifier column to keep track of the original matrix
    component_temporary_df$matrix_id <- i
    
    # Append the temporary data frame to the combined data frame
    component_df <- rbind(component_df, component_temporary_df)
  }
  
  # Return a list containing both the component data frame and gelman_results
  return(list(component_df = component_df, gelman_results = gelman_results))
}



###
#Autocorrelation test function
autocorr_analysis_function <- function(mcmc_list, component = "Sol") {
  # Apply the autocorr.diag function to each element of the list
  autocorr_list <- lapply(mcmc_list, function(x) {
    if (component == "Sol") {
      autocorr.diag(x$Sol)
    } else if (component == "VCV") {
      autocorr.diag(x$VCV[, 1]) #Again, needed to select the first column (the column for the additive genetic variance), rather than select both columns, as the second column is for the residual variance which is fixed at 1.
    } else {
      stop("Invalid component argument. Use 'Sol' or 'VCV'.")
    }
  })
  
  # Use lapply to extract the second row of each matrix
  second_row_list <- lapply(autocorr_list, function(x) {
    x[2,]
  })
  
  # Convert the list to a vector
  second_row_vec <- unlist(second_row_list)
  
  # Return the quantiles of the second row vector
  return(quantile(second_row_vec))
}


###
#ESS calculation function
calculate_ess_function <- function(mcmc_list, component = "Sol") {
  # Apply the effectiveSize function to each element of the list
  ess_list <- lapply(mcmc_list, function(x) {
    if (component == "Sol") {
      effectiveSize(x$Sol)
    } else if (component == "VCV") {
      effectiveSize(x$VCV[, 1])
    } else {
      stop("Invalid component argument. Use 'Sol' or 'VCV'.")
    }
  })
  
  # Convert the list to a vector
  ess_vec <- unlist(ess_list)
  
  # Return the quantiles of the ESS vector
  return(quantile(ess_vec))
}

###
# Function that runs the entire suite of diagnostic checks on a model
perform_analysis <- function(region) {
  # Read in the models
  mcmc_list_1 <- model_reader_function(region, "1st_run")
  mcmc_list_2 <- model_reader_function(region, "2nd_run")
  
  # Convergence tests
  sol_conv <- gelman_rubin_function(mcmc_list_1, mcmc_list_2, component = "Sol")
  vcv_conv <- gelman_rubin_function(mcmc_list_1, mcmc_list_2, component = "VCV")
  
  cat("Quantiles of Gelman-Rubin convergence tests for Sol component:\n")
  quantiles_sol_conv <- quantile(unlist(sol_conv$gelman_results))
  print(quantiles_sol_conv)
  
  cat("Quantiles of Gelman-Rubin convergence tests for VCV component:\n")
  quantiles_vcv_conv <- quantile(unlist(vcv_conv$gelman_results))
  print(quantiles_vcv_conv)
  
  cat("Gelman-Rubin component data frames for Sol and VCV components:\n")
  sol_conv_df <- sol_conv$component_df
  vcv_conv_df <- vcv_conv$component_df
  
  # Autocorrelation test
  sol_autocorr <- autocorr_analysis_function(mcmc_list_1, component = "Sol")
  vcv_autocorr <- autocorr_analysis_function(mcmc_list_1, component = "VCV")
  
  cat("Quantiles of autocorrelation tests for Sol component:\n")
  quantiles_sol_autocorr <- quantile(unlist(sol_autocorr))
  print(quantiles_sol_autocorr)
  
  cat("Quantiles of autocorrelation tests for VCV component:\n")
  quantiles_vcv_autocorr <- quantile(unlist(vcv_autocorr))
  print(quantiles_vcv_autocorr)
  
  # ESS test
  sol_ESS <- calculate_ess_function(mcmc_list = mcmc_list_1, component = "Sol")
  vcv_ESS <- calculate_ess_function(mcmc_list = mcmc_list_1, component = "VCV")
  
  cat("Quantiles of ESS tests for Sol component:\n")
  quantiles_sol_ESS <- quantile(unlist(sol_ESS))
  print(quantiles_sol_ESS)
  
  cat("Quantiles of ESS tests for VCV component:\n")
  quantiles_vcv_ESS <- quantile(unlist(vcv_ESS))
  print(quantiles_vcv_ESS)
  
  # Return the results as a list
  return(list(
    quantiles_sol_conv = quantiles_sol_conv,
    quantiles_vcv_conv = quantiles_vcv_conv,
    sol_conv_df = sol_conv_df,
    vcv_conv_df = vcv_conv_df,
    quantiles_sol_autocorr = quantiles_sol_autocorr,
    quantiles_vcv_autocorr = quantiles_vcv_autocorr,
    quantiles_sol_ESS = quantiles_sol_ESS,
    quantiles_vcv_ESS = quantiles_vcv_ESS
  ))
}

############################################################################################################################################################


############
#Global
############
CS_Caste_Global <- perform_analysis(region = "Global")

############
#Tropical#
############
CS_Caste_Tropical <- perform_analysis("Tropical")

############
#Temperate#
############
CS_Caste_Temperate <- perform_analysis("Temperate")

############
#Both#
############
CS_Caste_Both <- perform_analysis("Both")
