############################### Create Phylogenetic PCA analysis data set for each MCC tree ########################################
##Louis Bell-Roberts
#16/11/2023


# Packages
library(data.table)
library(ape)
library(ggtree) ## Basic ggplot package for trees
library(ggtreeExtra) ## For fancy tree plotting
library(ufs) ## For ggqq plot
library(phytools)
library(phylolm)
library(rr2)
library(ggplot2)
library(ggpubr)
library(dplyr)
theme_set(theme_minimal())

ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/AntClm_CSize_all_230710.csv")

#Transform temperature and rainfall with Z transformation
ant_data <- ant_data %>% 
  mutate(across(c(TMPavg_avg, PREavg_avg), .fns = ~scale(.x), .names = "Z{.col}"))

#Select only the variables that I want to be in the PCA analysis
ant_data_cs_pca<- dplyr::select(ant_data, ZTMPavg_avg, ZPREavg_avg)

#Must be in matrix format
ant_focal<- as.matrix(ant_data_cs_pca)

#Requires species names as row names
rownames(ant_focal) <- ant_data$valid_species_name

#Read in phylogenies, one for each MCC tree - required for phylogenetic PCA
NCuniform_stem_tree <- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/MCC_trees/pruned_474_sp_NCuniform_stem_MCC.tre")
NCuniform_crown_tree <- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/MCC_trees/pruned_474_sp_NCuniform_crown_MCC.tre")
FBD_stem_tree<- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/MCC_trees/pruned_474_sp_FBD_stem_MCC.tre")
FBD_crown_tree<- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/MCC_trees/pruned_474_sp_FBD_crown_MCC.tre")

#Perform phylogenetic PCA on each tree
pca_mean_NCuniform_stem <- phyl.pca(NCuniform_stem_tree, ant_focal, method="BM", mode="cov")
pca_mean_NCuniform_crown <- phyl.pca(NCuniform_crown_tree, ant_focal, method="BM", mode="cov")
pca_mean_FBD_stem <- phyl.pca(FBD_stem_tree, ant_focal, method="BM", mode="cov")
pca_mean_FBD_crown <- phyl.pca(FBD_crown_tree, ant_focal, method="BM", mode="cov")

## Plot and check
#Loadings: positive loading indicates a positive association between the original variable and the value of the corresponding principal component, and vice versa.

#NC_uniform_stem
biplot(pca_mean_NCuniform_stem, col=c("gray","black"), cex= c(0.5,1),
       xlim=c(-4,4), ylim=c(-4,4))
pca_mean_NCuniform_stem$L #Matrix of loadings: PC1: TMP: -VE; PRE: -VE; PC2: TMP: +VE; PRE: -VE
summary(pca_mean_NCuniform_stem) #PC1 explains 67% of variance

pca_mean_NCuniform_stem_df <- as.data.frame(pca_mean_NCuniform_stem$L)
predictors <- c("Average Temperature", "Average Rainfall")
pca_mean_NCuniform_stem_df <- data.frame(Predictors = predictors, pca_mean_NCuniform_stem_df)
pca_mean_NCuniform_stem_var_df <- summary(pca_mean_NCuniform_stem)
# write.csv(pca_mean_NCuniform_stem_df, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_loadings_NCuniform_stem.csv", row.names = F)
# write.csv(pca_mean_NCuniform_stem_var_df$importance, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_variation_explained_NCuniform_stem.csv", row.names = F)

#NC_uniform_crown
# biplot(pca_mean_NCuniform_crown, col=c("gray","black"), cex= c(0.5,1),
#        xlim=c(-4,4), ylim=c(-4,4))
pca_mean_NCuniform_crown$L #Matrix of loadings: PC1: TMP: +VE; PRE: +VE; PC2: TMP: -VE; PRE: +VE 
summary(pca_mean_NCuniform_crown) #PC1 explains 77% of variance

pca_mean_NCuniform_crown_df <- as.data.frame(pca_mean_NCuniform_crown$L)
pca_mean_NCuniform_crown_df <- data.frame(Predictors = predictors, pca_mean_NCuniform_crown_df)
pca_mean_NCuniform_crown_var_df <- summary(pca_mean_NCuniform_crown)
# write.csv(pca_mean_NCuniform_crown_df, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_loadings_NCuniform_crown.csv", row.names = F)
# write.csv(pca_mean_NCuniform_crown_var_df$importance, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_variation_explained_NCuniform_crown.csv", row.names = F)


#FBD_stem
# biplot(pca_mean_FBD_stem, col=c("gray","black"), cex= c(0.5,1),
#        xlim=c(-4,4), ylim=c(-4,4))
pca_mean_FBD_stem$L #Matrix of loadings: PC1: TMP: -VE; PRE: -VE; PC2: TMP: +VE; PRE: -VE
summary(pca_mean_FBD_stem) #PC1 explains 64% of variance

pca_mean_FBD_stem_df <- as.data.frame(pca_mean_FBD_stem$L)
pca_mean_FBD_stem_df <- data.frame(Predictors = predictors, pca_mean_FBD_stem_df)
pca_mean_FBD_stem_var_df <- summary(pca_mean_FBD_stem)
# write.csv(pca_mean_FBD_stem_df, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_loadings_FBD_stem.csv", row.names = F)
# write.csv(pca_mean_FBD_stem_var_df$importance, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_variation_explained_FBD_stem.csv", row.names = F)

#FBD_crown
# biplot(pca_mean_FBD_crown, col=c("gray","black"), cex= c(0.5,1),
#        xlim=c(-4,4), ylim=c(-4,4))
pca_mean_FBD_crown$L #Matrix of loadings: PC1: TMP: -VE; PRE: -VE; PC2: TMP: +VE; PRE: -VE
summary(pca_mean_FBD_crown) #PC1 explains 67% of variance

pca_mean_FBD_crown_df <- as.data.frame(pca_mean_FBD_crown$L)
pca_mean_FBD_crown_df <- data.frame(Predictors = predictors, pca_mean_FBD_crown_df)
pca_mean_FBD_crown_var_df <- summary(pca_mean_FBD_crown)
# write.csv(pca_mean_FBD_crown_df, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_loadings_FBD_crown.csv", row.names = F)
# write.csv(pca_mean_FBD_crown_var_df$importance, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_variation_explained_FBD_crown.csv", row.names = F)


#Create combined table
pca_load_NCuniform_stem <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_loadings_NCuniform_stem.csv")
pca_load_NCuniform_crown <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_loadings_NCuniform_crown.csv")
pca_load_FBD_stem <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_loadings_FBD_stem.csv")
pca_load_FBD_crown <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_loadings_FBD_crown.csv")

MCC_tree <- c(rep("NCuniform_stem", 4), rep("NCuniform_crown", 4), rep("FBD_stem", 4), rep("FBD_crown", 4))
pca_load_combined <- rbind(pca_load_NCuniform_stem, pca_load_NCuniform_crown, pca_load_FBD_stem, pca_load_FBD_crown)
pca_load_combined <- data.frame(pca_load_combined, MCC_tree = MCC_tree)
pca_load_combined <- rename(pca_load_combined, Variable = Predictors)
pca_load_combined$PC1 <- round(pca_load_combined$PC1, digits = 2)
pca_load_combined$PC2 <- round(pca_load_combined$PC2, digits = 2)

write.csv(pca_load_combined, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Supplementary_materials/Tables/pca_loadings_combined.csv", row.names = F)

## Make the dataset, there is only 2 PCs because of two variables

################
#Global datasets
################

#NC_uniform_stem
NCuniform_stem_PCA_data <- data.frame(ant_data,
                       PC1= pca_mean_NCuniform_stem$S[,1][match(ant_data$valid_species_name, rownames(pca_mean_NCuniform_stem$S))],
                       PC2= pca_mean_NCuniform_stem$S[,2][match(ant_data$valid_species_name, rownames(pca_mean_NCuniform_stem$S))])
# write.csv(x = NCuniform_stem_PCA_data, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/PCA_MCC_trees/NCuniform_stem_all_PCA_231116.csv", row.names = F)

#NC_uniform_crown
NCuniform_crown_PCA_data <- data.frame(ant_data,
                       PC1= pca_mean_NCuniform_crown$S[,1][match(ant_data$valid_species_name, rownames(pca_mean_NCuniform_crown$S))],
                       PC2= pca_mean_NCuniform_crown$S[,2][match(ant_data$valid_species_name, rownames(pca_mean_NCuniform_crown$S))])
# write.csv(x = NCuniform_crown_PCA_data, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/PCA_MCC_trees/NCuniform_crown_all_PCA_231116.csv", row.names = F)
#FBD_stem
FBD_stem_PCA_data <- data.frame(ant_data,
                       PC1= pca_mean_FBD_stem$S[,1][match(ant_data$valid_species_name, rownames(pca_mean_FBD_stem$S))],
                       PC2= pca_mean_FBD_stem$S[,2][match(ant_data$valid_species_name, rownames(pca_mean_FBD_stem$S))])
# write.csv(x = FBD_stem_PCA_data, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/PCA_MCC_trees/FBD_stem_all_PCA_231116.csv", row.names = F)
#FBD_crown
FBD_crown_PCA_data <- data.frame(ant_data,
                       PC1= pca_mean_FBD_crown$S[,1][match(ant_data$valid_species_name, rownames(pca_mean_FBD_crown$S))],
                       PC2= pca_mean_FBD_crown$S[,2][match(ant_data$valid_species_name, rownames(pca_mean_FBD_crown$S))])
# write.csv(x = FBD_crown_PCA_data, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/PCA_MCC_trees/FBD_crown_all_PCA_231116.csv", row.names = F)


##################
#Regional datasets
##################

#Read in regional datasets so that I can match species
ant_data_tropical<- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/Tropical_ant_data.csv")
ant_data_temperate<- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/Temperate_ant_data_230804.csv")
ant_data_both<- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/Both_ant_data_230804.csv")

###Tropical
#NC_uniform_stem
ant_data_tropical_PCA_NCuniform_stem <- NCuniform_stem_PCA_data %>% dplyr::filter(valid_species_name %in% ant_data_tropical$valid_species_name)
#NC_uniform_crown
ant_data_tropical_PCA_NCuniform_crown <- NCuniform_crown_PCA_data %>% dplyr::filter(valid_species_name %in% ant_data_tropical$valid_species_name)
#FBD_stem
ant_data_tropical_PCA_FBD_stem <- FBD_stem_PCA_data %>% dplyr::filter(valid_species_name %in% ant_data_tropical$valid_species_name)
#FBD_crown
ant_data_tropical_PCA_FBD_crown <- FBD_crown_PCA_data %>% dplyr::filter(valid_species_name %in% ant_data_tropical$valid_species_name)

###Temperate
#NC_uniform_stem
ant_data_temperate_PCA_NCuniform_stem <- NCuniform_stem_PCA_data %>% dplyr::filter(valid_species_name %in% ant_data_temperate$valid_species_name)
#NC_uniform_crown
ant_data_temperate_PCA_NCuniform_crown <- NCuniform_crown_PCA_data %>% dplyr::filter(valid_species_name %in% ant_data_temperate$valid_species_name)
#FBD_stem
ant_data_temperate_PCA_FBD_stem <- FBD_stem_PCA_data %>% dplyr::filter(valid_species_name %in% ant_data_temperate$valid_species_name)
#FBD_crown
ant_data_temperate_PCA_FBD_crown <- FBD_crown_PCA_data %>% dplyr::filter(valid_species_name %in% ant_data_temperate$valid_species_name)


###Both
#NC_uniform_stem
ant_data_both_PCA_NCuniform_stem <- NCuniform_stem_PCA_data %>% dplyr::filter(valid_species_name %in% ant_data_both$valid_species_name)
#NC_uniform_crown
ant_data_both_PCA_NCuniform_crown <- NCuniform_crown_PCA_data %>% dplyr::filter(valid_species_name %in% ant_data_both$valid_species_name)
#FBD_stem
ant_data_both_PCA_FBD_stem <- FBD_stem_PCA_data %>% dplyr::filter(valid_species_name %in% ant_data_both$valid_species_name)
#FBD_crown
ant_data_both_PCA_FBD_crown <- FBD_crown_PCA_data %>% dplyr::filter(valid_species_name %in% ant_data_both$valid_species_name)


###Save files
##Tropical
# write.csv(x = ant_data_tropical_PCA_NCuniform_stem, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/PCA_MCC_trees_data/Tropical/NCuniform_stem_tropical_PCA_231117.csv") #NC_uniform_stem
# write.csv(x = ant_data_tropical_PCA_NCuniform_crown, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/PCA_MCC_trees_data/Tropical/NCuniform_crown_tropical_PCA_231117.csv") #NC_uniform_crown
# write.csv(x = ant_data_tropical_PCA_FBD_stem, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/PCA_MCC_trees_data/Tropical/FBD_stem_tropical_PCA_231117.csv") #FBD_stem
# write.csv(x = ant_data_tropical_PCA_FBD_crown, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/PCA_MCC_trees_data/Tropical/FBD_crown_tropical_PCA_231117.csv") #FBD_crown

##Temperate
# write.csv(x = ant_data_temperate_PCA_NCuniform_stem, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/PCA_MCC_trees_data/Temperate/NCuniform_stem_temperate_PCA_231117.csv") #NC_uniform_stem
# write.csv(x = ant_data_temperate_PCA_NCuniform_crown, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/PCA_MCC_trees_data/Temperate/NCuniform_crown_temperate_PCA_231117.csv") #NC_uniform_crown
# write.csv(x = ant_data_temperate_PCA_FBD_stem, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/PCA_MCC_trees_data/Temperate/FBD_stem_temperate_PCA_231117.csv") #FBD_stem
# write.csv(x = ant_data_temperate_PCA_FBD_crown, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/PCA_MCC_trees_data/Temperate/FBD_crown_temperate_PCA_231117.csv") #FBD_crown

##Both
# write.csv(x = ant_data_both_PCA_NCuniform_stem, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/PCA_MCC_trees_data/Both/NCuniform_stem_both_PCA_231117.csv") #NC_uniform_stem
# write.csv(x = ant_data_both_PCA_NCuniform_crown, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/PCA_MCC_trees_data/Both/NCuniform_crown_both_PCA_231117.csv") #NC_uniform_crown
# write.csv(x = ant_data_both_PCA_FBD_stem, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/PCA_MCC_trees_data/Both/FBD_stem_both_PCA_231117.csv") #FBD_stem
# write.csv(x = ant_data_both_PCA_FBD_crown, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/PCA_MCC_trees_data/Both/FBD_crown_both_PCA_231117.csv") #FBD_crown



