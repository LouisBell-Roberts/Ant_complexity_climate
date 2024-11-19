library(tidyverse)

# Read the dataset from CSV
ant_occ_pre_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Kass_2022/Kass_data_with_island.csv", header = T)
length(ant_occ_pre_data$lon_opt) #1479293

# Replace all "." with "_" in the species column
ant_occ_pre_data$valid_species_name <- gsub("\\.", "_", ant_occ_pre_data$valid_species_name)

#Remove ant geographic occurrences on islands
ant_occ_pre_data <- ant_occ_pre_data %>% filter(Island == 0)

#Rename columns
ant_occ_data <- ant_occ_pre_data %>% rename(dec_lat = lat_opt)

# Create a new column named 'Region' and initialize it with 'NA'
ant_occ_data$Region <- NA

# Update the 'Region' column based on latitude values
ant_occ_data$Region[ant_occ_data$dec_lat >= -23.5 & ant_occ_data$dec_lat <= 23.5] <- "Tropical"
ant_occ_data$Region[is.na(ant_occ_data$Region)] <- "Temperate"

# Aggregate the occurrence data based on species name
aggregated_data <- aggregate(Region ~ valid_species_name, data = ant_occ_data, function(x) {
  if (all(x == "Tropical")) {
    "Tropical"
  } else if (all(x == "Temperate")) {
    "Temperate"
  } else {
    "Both"
  }
})


setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/")

# Create a new dataframe for binary classification
region_2cat <- data.frame(species = aggregated_data$valid_species_name, Region = ifelse(aggregated_data$Region == "Tropical", "Tropical", "Temperate_both"))

Tropical_sp <- filter(region_2cat, Region == 'Tropical')
write.csv(Tropical_sp, file = "AntEnv_tropical_sp_241103.csv")

# Create 2 dataframes: 1) only species that occur in 'Temperate' environments; 2) only species that occur in 'Both' environments
Temperate_sp <- filter(aggregated_data, Region == 'Temperate')
Both_sp <- filter(aggregated_data, Region == 'Both')
write.csv(Temperate_sp, file = "AntEnv_temperate_sp_241103.csv")
write.csv(Both_sp, file = "AntEnv_both_sp_241103.csv")




###########################
#Add in the species trait data for CS, caste number and climatic variables
###########################

#Read in ant trait data
ant_trait_data <- read.csv(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/Post_koos_AntEnv_merged_climate_trait_data.csv", header = T)

tropical_species_data <- read.csv(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/AntEnv_tropical_sp_241103.csv") #read in just tropical species
temperate_species_data <- read.csv(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/AntEnv_temperate_sp_241103.csv") #read in just temperate species
both_species_data <- read.csv(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/AntEnv_both_sp_241103.csv") #read in just 'both' species


# Filter the ant_trait_data dataset to include only species that are tropical
ant_traits_tropical <- filter(ant_trait_data, species %in% tropical_species_data$species)
length(ant_traits_tropical$species) #104

# Filter the ant_trait_data dataset to include only species that are temperate
ant_traits_temperate <- filter(ant_trait_data, species %in% temperate_species_data$valid_species_name)
length(ant_traits_temperate$species) #171

# Filter the ant_trait_data dataset to include only species that are 'both'
ant_traits_both <- filter(ant_trait_data, species %in% both_species_data$valid_species_name)
length(ant_traits_both$species) #226


#Write databases to csv, then move to the pruning script, should end up with 358 and 116 species for the tropical and temp_both
write.csv(ant_traits_tropical, file = "AntTrait_Tropical_241103.csv", row.names = F)
write.csv(ant_traits_temperate, file = "AntTrait_Temperate_241103.csv", row.names = F)
write.csv(ant_traits_both, file = "AntTrait_Both_241103.csv", row.names = F)










