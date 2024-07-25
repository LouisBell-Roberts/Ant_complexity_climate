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
theme_set(theme_minimal())

ant_data<- fread("data/AntEnv_allCS_230608.csv")
ant_data2<- dplyr::mutate(ant_data, remove= dplyr::case_when(
  Supercolony > 0 | Parasitic > 0 | hybridisation > 0 | Parthenogenesis > 0 | Gamergates > 0 ~ 1,
  Supercolony < 1 | Parasitic < 1 | hybridisation < 1 | Parthenogenesis < 1 | Gamergates < 1 ~ 0
))
ant_data2$remove<- as.factor(ant_data2$remove)

ant_data_cs_pca<- dplyr::select(ant_data2, colony.size.y, pre_avg_avg, DTR_avg_avg,temp_sd_avg) 
rownames(ant_data_cs_pca)<- ant_data2$valid_species_name # For some reason not working. : (
#
CS_tree<- read.tree("data/Sample_pruned_allCS_230608.tre")
#
test_pca<- phyl.pca(CS_tree, ant_data_cs_pca, method="BM", mode="cov")
#
