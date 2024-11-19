########### Generate Ant Environment final dataset ###########
##Merge ant trait data and climatic data
##Prune phylogenies
##Prune database to match phylogenies
#Louis Bell-Roberts
#26/10/2024

library(tidyverse)
library(ape)
library(phytools)

##############################
#Merge ant trait data and climate data
##############################

#Read in the ant trait data
ant_trait_pre_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/Post_Koos_AntEnv_data_FINAL.csv")

#Filter to only include species with complete data
ant_trait_data <- ant_trait_pre_data %>% filter(complete.cases(colony.size), complete.cases(poly_id))

#Read in the climatic data
climate_pre_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/Kass_Ant_Climate_Data_1015.csv")

#Rename the species column
climate_data <- climate_pre_data %>% rename(
  species = valid_species_name
)

#Replace "." with "_" within species names
climate_data$species <- gsub(pattern = "\\.", replacement = "_", climate_data$species)

#Merge datasets
data_combined <- merge(ant_trait_data, climate_data, by = "species", all.x = TRUE)


##############################
#Prune phylogenies to match database
##############################

#Read in 400 ant trees
ant_trees <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/Economo_2018_400.tre")

#Read in the 4 MCCtrees
NCuniform_crown <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15K_NCuniform_crown_mcc.tre")

NCuniform_stem <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15k_NCuniform_stem_mcc.tre")

FBD_crown <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15K_FBD_crown_mcc.tre")

FBD_stem <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15K_FBD_stem_mcc.tre")

MCC_combined <- c(NCuniform_crown, NCuniform_stem, FBD_crown, FBD_stem)

#Prune 400 trees
ant_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, data_combined$species))

#Prune MCC trees
MCC_combined_pruned <- drop.tip.multiPhylo(MCC_combined, setdiff(MCC_combined[[1]]$tip.label, data_combined$species))


##############################
#Prune database to match trees
##############################

# Filter species in data_combined to match the species in ant_trees_pruned[[1]]$tip.label
data_combined_final <- data_combined %>%
  filter(species %in% ant_trees_pruned[[1]]$tip.label)



##############################
#Write database and trees to file
##############################

#Write final merged dataset to csv file
write.csv(data_combined_final, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/Post_koos_AntEnv_merged_climate_trait_data.csv", row.names = F)

#Write 400 trees to file
write.tree(phy = ant_trees_pruned, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/Post_Koos/400_trees/400trees_pruned_271024.tre")

#Write MCC trees to file
write.tree(phy = MCC_combined_pruned, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/Post_Koos/MCC_trees/MCCtrees_pruned_271024.tre")
