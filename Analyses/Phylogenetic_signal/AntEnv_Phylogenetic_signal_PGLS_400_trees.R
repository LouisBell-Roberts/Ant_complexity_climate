########AntEnv phylogenetic signal in climatic variables - analysing species globally, tropical, temperate and both######
#####Using PGLS and R2
####Only analysing average climatic variables
###Written for Louis' laptop with parallel computation
##Data and trees have already been pruned
#Louis Bell-Roberts
#15/09/2023

# Loading essential packages ####
library(data.table)
library(ape)
library(phylolm)
library(doParallel)
library(phytools)
registerDoParallel(8)

setwd("/Users/louis.bell-roberts/Documents/Github/Division-of-labour/Phylogenetic_analysis/Ant_environmental_data/Analyses/Phylogenetic_signal/Outputs/")
##############
#Global
##############

ant_data_global<- fread("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/AntClm_CSize_all_230710.csv")
# For PhyloLM data structure
ant_data2_global<- ant_data_global
rownames(ant_data2_global)<- ant_data2_global$valid_species_name
# Read in the trees
tree_set<- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/400_trees/400trees_pruned_CSize_230712.tre") 
# Temp variables
iterations<- length(tree_set)

global_TMP_lm<- lm(TMPavg_avg~ 1, data=ant_data2_global)
global_PRE_lm<- lm(sqrt(PREavg_avg)~ 1, data=ant_data2_global)
global_DTR_lm<- lm(log(DTRavg_avg)~ 1, data=ant_data2_global)

global_out<- foreach(i=1:iterations, .combine=rbind) %dopar%{
  focal_tree<- tree_set[[i]]
  
  global_TMP <- phylolm(TMPavg_avg ~ 1, 
                          data=ant_data2_global, phy=focal_tree, model = "lambda") #
  global_PRE <- phylolm(sqrt(PREavg_avg) ~ 1, 
                          data=ant_data2_global, phy=focal_tree, model = "lambda") #
  global_DTR <- phylolm(log(DTRavg_avg) ~ 1, 
                          data=ant_data2_global, phy=focal_tree, model = "lambda") #
  #
  print(c(i,
          round(R2_lik(mod= global_TMP, mod.r= global_TMP_lm),3),
          round(R2_lik(mod= global_PRE, mod.r= global_PRE_lm),3),
          round(R2_lik(mod= global_DTR, mod.r= global_DTR_lm),3)
  ))
}
#
colnames(global_out) = c('TreeID','TMP_RandomR2','PRE_RandomR2','DTR_RandomR2')
global_out<- as.data.frame(global_out)
# identify the model with the lowest average AIC value
colMeans(global_out) #FullMod is the best model
#
fwrite(global_out, "PGLS_global_out_R2.txt")
#

###############################################################################

##############
#Tropical
##############

ant_data_tropical<- fread("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/Tropical_ant_data.csv")
# For PhyloLM data structure
ant_data2_tropical<- ant_data_tropical
rownames(ant_data2_tropical)<- ant_data2_tropical$valid_species_name
# Read in the trees
tree_set<- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/Tropical_vs_temperate/400_trees/Tropical_ant_trees_pruned.tre") 
# Temp variables
iterations<- length(tree_set)

tropical_TMP_lm<- lm(TMPavg_avg~ 1, data=ant_data2_tropical)
tropical_PRE_lm<- lm(sqrt(PREavg_avg)~ 1, data=ant_data2_tropical)
tropical_DTR_lm<- lm(log(DTRavg_avg)~ 1, data=ant_data2_tropical)

tropical_out<- foreach(i=1:iterations, .combine=rbind) %dopar%{
  focal_tree<- tree_set[[i]]
  
  tropical_TMP <- phylolm(TMPavg_avg ~ 1, 
                           data=ant_data2_tropical, phy=focal_tree, model = "lambda") #
  tropical_PRE <- phylolm(sqrt(PREavg_avg) ~ 1, 
                           data=ant_data2_tropical, phy=focal_tree, model = "lambda") #
  tropical_DTR <- phylolm(log(DTRavg_avg) ~ 1, 
                           data=ant_data2_tropical, phy=focal_tree, model = "lambda") #
  #
  print(c(i,
          round(R2_lik(mod= tropical_TMP, mod.r= tropical_TMP_lm),3),
          round(R2_lik(mod= tropical_PRE, mod.r= tropical_PRE_lm),3),
          round(R2_lik(mod= tropical_DTR, mod.r= tropical_DTR_lm),3)
  ))
}
#
colnames(tropical_out) = c('TreeID','TMP_RandomR2','PRE_RandomR2','DTR_RandomR2')
tropical_out<- as.data.frame(tropical_out)
# identify the model with the lowest average AIC value
colMeans(tropical_out) #FullMod is the best model
#
fwrite(tropical_out, "PGLS_tropical_out_R2.txt")
#

###############################################################################

##############
#Temperate
##############

ant_data_temperate<- fread("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/Temperate_ant_data_230804.csv")
# For PhyloLM data structure
ant_data2_temperate<- ant_data_temperate
rownames(ant_data2_temperate)<- ant_data2_temperate$valid_species_name
# Read in the trees
tree_set<- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/Tropical_vs_temperate/400_trees/Temperate_ant_trees_pruned.tre") 
# Temp variables
iterations<- length(tree_set)

temperate_TMP_lm<- lm(TMPavg_avg~ 1, data=ant_data2_temperate)
temperate_PRE_lm<- lm(sqrt(PREavg_avg)~ 1, data=ant_data2_temperate)
temperate_DTR_lm<- lm(log(DTRavg_avg)~ 1, data=ant_data2_temperate)

temperate_out<- foreach(i=1:iterations, .combine=rbind) %dopar%{
  focal_tree<- tree_set[[i]]
  
  temperate_TMP <- phylolm(TMPavg_avg ~ 1, 
                      data=ant_data2_temperate, phy=focal_tree, model = "lambda") #
  temperate_PRE <- phylolm(sqrt(PREavg_avg) ~ 1, 
                      data=ant_data2_temperate, phy=focal_tree, model = "lambda") #
  temperate_DTR <- phylolm(log(DTRavg_avg) ~ 1, 
                      data=ant_data2_temperate, phy=focal_tree, model = "lambda") #
  #
  print(c(i,
          round(R2_lik(mod= temperate_TMP, mod.r= temperate_TMP_lm),3),
          round(R2_lik(mod= temperate_PRE, mod.r= temperate_PRE_lm),3),
          round(R2_lik(mod= temperate_DTR, mod.r= temperate_DTR_lm),3)
  ))
}
#
colnames(temperate_out) = c('TreeID','TMP_RandomR2','PRE_RandomR2','DTR_RandomR2')
temperate_out<- as.data.frame(temperate_out)
# identify the model with the lowest average AIC value
colMeans(temperate_out) #FullMod is the best model
#
fwrite(temperate_out, "PGLS_temperate_R2.txt")
#

###############################################################################

##############
#Both
##############

ant_data_both<- fread("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/Both_ant_data_230804.csv")
# For PhyloLM data structure
ant_data2_both<- ant_data_both
rownames(ant_data2_both)<- ant_data2_both$valid_species_name
# Read in the trees
tree_set<- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/Tropical_vs_temperate/400_trees/Both_ant_trees_pruned.tre") 
# Temp variables
iterations<- length(tree_set)

both_TMP_lm<- lm(TMPavg_avg~ 1, data=ant_data2_both)
both_PRE_lm<- lm(sqrt(PREavg_avg)~ 1, data=ant_data2_both)
both_DTR_lm<- lm(log(DTRavg_avg)~ 1, data=ant_data2_both)

both_out<- foreach(i=1:iterations, .combine=rbind) %dopar%{
  focal_tree<- tree_set[[i]]
  
  both_TMP <- phylolm(TMPavg_avg ~ 1, 
                       data=ant_data2_both, phy=focal_tree, model = "lambda") #
  both_PRE <- phylolm(sqrt(PREavg_avg) ~ 1, 
                      data=ant_data2_both, phy=focal_tree, model = "lambda") #
  both_DTR <- phylolm(log(DTRavg_avg) ~ 1, 
                      data=ant_data2_both, phy=focal_tree, model = "lambda") #
  #
  print(c(i,
          round(R2_lik(mod= both_TMP, mod.r= both_TMP_lm),3),
          round(R2_lik(mod= both_PRE, mod.r= both_PRE_lm),3),
          round(R2_lik(mod= both_DTR, mod.r= both_DTR_lm),3)
  ))
}
#
colnames(both_out) = c('TreeID','TMP_RandomR2','PRE_RandomR2','DTR_RandomR2')
both_out<- as.data.frame(both_out)
# identify the model with the lowest average AIC value
colMeans(both_out) #FullMod is the best model
#
fwrite(both_out, "PGLS_both_R2.txt")
#

###############################################################################