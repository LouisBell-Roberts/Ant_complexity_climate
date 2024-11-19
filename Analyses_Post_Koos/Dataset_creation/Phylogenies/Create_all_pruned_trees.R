#Prune databases and prune trees for AntEnv analyses
library(tidyverse)
library(ape)
library(phytools)

#Read in data file
ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/Final_datasets/Post_koos_AntEnv_merged_climate_trait_data.csv") #After species with unusual biology removed = 501 species

tropical_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/Final_datasets/AntTrait_Tropical_241103.csv", header = T)
temperate_sp_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/Final_datasets/AntTrait_Temperate_241103.csv", header = T)
both_sp_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/Final_datasets/AntTrait_Both_241103.csv", header = T)

# Replace all "_" with "." in the species column
ant_data$species <- gsub("_", ".", ant_data$species)
tropical_data$species <- gsub("_", ".", tropical_data$species)
temperate_sp_data$species <- gsub("_", ".", temperate_sp_data$species)
both_sp_data$species <- gsub("_", ".", both_sp_data$species)


#Read in the 4 MCC trees
anttree_NCuniform_stem <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Economo_2018/Dryad_archive/15k_all_ant_trees/15k_NCuniform_stem_mcc.tre")
anttree_NCuniform_crown <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Economo_2018/Dryad_archive/15k_all_ant_trees/15K_NCuniform_crown_mcc.tre")
anttree_FBD_stem <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Economo_2018/Dryad_archive/15k_all_ant_trees/15K_FBD_stem_mcc.tre")
anttree_FBD_crown <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Economo_2018/Dryad_archive/15k_all_ant_trees/15K_FBD_crown_mcc.tre")

#Read in the 400 trees
anttree_NCuniform_stem_100 <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Economo_2018/Dryad_archive/15k_all_ant_trees/15k_NCuniform_stem_posterior.tre")
anttree_NCuniform_crown_100 <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Economo_2018/Dryad_archive/15k_all_ant_trees/15k_NCuniform_crown_posterior.tre")
anttree_FBD_stem_100 <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Economo_2018/Dryad_archive/15k_all_ant_trees/15K_FBD_stem_posterior.tre")
anttree_FBD_crown_100 <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Economo_2018/Dryad_archive/15k_all_ant_trees/15K_FBD_crown_posterior.tre")

ant_400_trees <- c(anttree_NCuniform_stem_100, anttree_NCuniform_crown_100, anttree_FBD_stem_100, anttree_FBD_crown_100)



#####################
#Prune the MCC trees
#####################

#Prune trees - all 501
pruned_all_sp_NCuniform_stem <-drop.tip(anttree_NCuniform_stem, setdiff(anttree_NCuniform_stem$tip.label, ant_data$species))
pruned_all_sp_NCuniform_crown <-drop.tip(anttree_NCuniform_crown, setdiff(anttree_NCuniform_crown$tip.label, ant_data$species))
pruned_all_sp_FBD_stem <-drop.tip(anttree_FBD_stem, setdiff(anttree_FBD_stem$tip.label, ant_data$species))
pruned_all_sp_FBD_crown <-drop.tip(anttree_FBD_crown, setdiff(anttree_FBD_crown$tip.label, ant_data$species))

#Prune trees - tropical
pruned_tropical_NCuniform_stem <-drop.tip(anttree_NCuniform_stem, setdiff(anttree_NCuniform_stem$tip.label, tropical_data$species))
pruned_tropical_NCuniform_crown <-drop.tip(anttree_NCuniform_crown, setdiff(anttree_NCuniform_crown$tip.label, tropical_data$species))
pruned_tropical_FBD_stem <-drop.tip(anttree_FBD_stem, setdiff(anttree_FBD_stem$tip.label, tropical_data$species))
pruned_tropical_FBD_crown <-drop.tip(anttree_FBD_crown, setdiff(anttree_FBD_crown$tip.label, tropical_data$species))

#Prune trees - temperate
pruned_temperate_NCuniform_stem <-drop.tip(anttree_NCuniform_stem, setdiff(anttree_NCuniform_stem$tip.label, temperate_sp_data$species))
pruned_temperate_NCuniform_crown <-drop.tip(anttree_NCuniform_crown, setdiff(anttree_NCuniform_crown$tip.label, temperate_sp_data$species))
pruned_temperate_FBD_stem <-drop.tip(anttree_FBD_stem, setdiff(anttree_FBD_stem$tip.label, temperate_sp_data$species))
pruned_temperate_FBD_crown <-drop.tip(anttree_FBD_crown, setdiff(anttree_FBD_crown$tip.label, temperate_sp_data$species))

#Prune trees - both
pruned_both_NCuniform_stem <-drop.tip(anttree_NCuniform_stem, setdiff(anttree_NCuniform_stem$tip.label, both_sp_data$species))
pruned_both_NCuniform_crown <-drop.tip(anttree_NCuniform_crown, setdiff(anttree_NCuniform_crown$tip.label, both_sp_data$species))
pruned_both_FBD_stem <-drop.tip(anttree_FBD_stem, setdiff(anttree_FBD_stem$tip.label, both_sp_data$species))
pruned_both_FBD_crown <-drop.tip(anttree_FBD_crown, setdiff(anttree_FBD_crown$tip.label, both_sp_data$species))



#####################
#Prune the 400 trees
#####################

#Prune trees - all 501
pruned_all_sp_400trees <-drop.tip.multiPhylo(ant_400_trees, setdiff(ant_400_trees[[1]]$tip.label, ant_data$species))

#Prune trees - tropical 104
pruned_tropical_400trees <-drop.tip.multiPhylo(ant_400_trees, setdiff(ant_400_trees[[1]]$tip.label, tropical_data$species))

#Prune trees - temperate 171
pruned_temperate_400trees <-drop.tip.multiPhylo(ant_400_trees, setdiff(ant_400_trees[[1]]$tip.label, temperate_sp_data$species))

#Prune trees - both 226
pruned_both_400trees <-drop.tip.multiPhylo(ant_400_trees, setdiff(ant_400_trees[[1]]$tip.label, both_sp_data$species))




#####################
#Replace "." with "_" in all trees
#####################

###MCC trees

# List of pruned trees
MCCtrees_vector <- list(
  pruned_all_sp_NCuniform_stem,
  pruned_all_sp_NCuniform_crown,
  pruned_all_sp_FBD_stem,
  pruned_all_sp_FBD_crown,
  pruned_tropical_NCuniform_stem,
  pruned_tropical_NCuniform_crown,
  pruned_tropical_FBD_stem,
  pruned_tropical_FBD_crown,
  pruned_temperate_NCuniform_stem,
  pruned_temperate_NCuniform_crown,
  pruned_temperate_FBD_stem,
  pruned_temperate_FBD_crown,
  pruned_both_NCuniform_stem,
  pruned_both_NCuniform_crown,
  pruned_both_FBD_stem,
  pruned_both_FBD_crown
)

# Replace all instances of "." with "_" in the tip labels
MCCtrees_vector <- lapply(MCCtrees_vector, function(tree) {
  tree$tip.label <- gsub("\\.", "_", tree$tip.label)
  return(tree)
})

###400 trees

#Global
# Replace "." with "_" in all tip labels within a multiPhylo object
pruned_all_sp_400trees <- lapply(pruned_all_sp_400trees, function(tree) {
  tree$tip.label <- sub("\\.", "_", tree$tip.label)
  tree
})

#Tropical
pruned_tropical_400trees <- lapply(pruned_tropical_400trees, function(tree) {
  tree$tip.label <- sub("\\.", "_", tree$tip.label)
  tree
})

#Temperate
pruned_temperate_400trees <- lapply(pruned_temperate_400trees, function(tree) {
  tree$tip.label <- sub("\\.", "_", tree$tip.label)
  tree
})

#Both
pruned_both_400trees <- lapply(pruned_both_400trees, function(tree) {
  tree$tip.label <- sub("\\.", "_", tree$tip.label)
  tree
})


#####################
#Write trees to file
#####################

###
#MCC trees

# Set working directory
setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/Post_Koos/MCC_trees/")

# Create a vector of filenames corresponding to each tree
filenames <- c(
  "NCuniform_stem_Global_pruned_MCC.tre",
  "NCuniform_crown_Global_pruned_MCC.tre",
  "FBD_stem_Global_pruned_MCC.tre",
  "FBD_crown_Global_pruned_MCC.tre",
  "NCuniform_stem_tropical_pruned_MCC.tre",
  "NCuniform_crown_tropical_pruned_MCC.tre",
  "FBD_stem_tropical_pruned_MCC.tre",
  "FBD_crown_tropical_pruned_MCC.tre",
  "NCuniform_stem_temperate_pruned_MCC.tre",
  "NCuniform_crown_temperate_pruned_MCC.tre",
  "FBD_stem_temperate_pruned_MCC.tre",
  "FBD_crown_temperate_pruned_MCC.tre",
  "NCuniform_stem_both_pruned_MCC.tre",
  "NCuniform_crown_both_pruned_MCC.tre",
  "FBD_stem_both_pruned_MCC.tre",
  "FBD_crown_both_pruned_MCC.tre"
)

# Write phylogenetic trees to file in a loop
mapply(write.tree, MCCtrees_vector, filenames)


###
##400 trees

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/Post_Koos/400_trees/")

###All species
# write.tree(pruned_all_sp_400trees, file = "Global_ant_trees_pruned.tre")

###Tropical species
write.tree(pruned_tropical_400trees, file = "Tropical_ant_trees_pruned.tre")

###Temperate species
write.tree(pruned_temperate_400trees, file = "Temperate_ant_trees_pruned.tre")

###Both species
write.tree(pruned_both_400trees, file = "Both_ant_trees_pruned.tre")
