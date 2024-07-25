########## R2 calculations for caste and colony size using MCMCglmm - Shinichi method ###########
####Binary and continuous data
###Result of script - R2 values are quite high 70/80%. Much greater than when using phylolm/glm or phyr
##20/10/2023
#Louis Bell-Roberts


#Loading essential packages ####
library(data.table)
library(ape)
library(phylolm)
library(doParallel)
library(phytools)
library(MCMCglmm)
registerDoParallel(8)
#
setwd('/Users/louis.bell-roberts/Documents/Github/Division-of-labour/Phylogenetic_analysis/Ant_environmental_data/Analyses/Regressions/Caste_response/Model_selection_AICaveraging_method/230816_avg_models_only_model_selection/All_sp_474/Outputs/')
#
ant_data2<- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/AntClm_CSize_all_230710.csv")

# Read in the trees
tree_set<- read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/400_trees/400trees_pruned_CSize_230712.tre") 

ant_data2 <- ant_data2 %>% dplyr::rename(animal = valid_species_name)
ant_data2 <- as.data.frame(ant_data2)

#MCMCglmm model
prior1 <- list(R = list(V = 1, fix = 1), G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

mm0 <- MCMCglmm(poly_id ~ 1, random = ~animal, 
         data = ant_data2, family="categorical", pedigree = tree_set[[1]], nitt = 110000, burnin = 10000, thin = 100, prior = prior1, verbose = F)
summary(mm0)

mmF <- MCMCglmm(poly_id~ TMPavg_avg+ sqrt(PREavg_avg)+ log(DTRavg_avg),
         data = ant_data2, random = ~animal, family = "categorical", pedigree = tree_set[[1]], nitt = 110000, burnin = 10000, thin = 100, prior = prior1, verbose = F)

summary(testMCMC)

c2 <- (16 * sqrt(3)/(15 * pi))^2
mmF$Sol1 <- mmF$Sol/sqrt(1+c2 * mmF$VCV[, "units"]) 
mmF$VCV1 <- mmF$VCV/(1+c2 * mmF$VCV[, "units"]) 

mFixed <- mean(mmF$Sol1[,2]) * mmF$X[, 2] + mean(mmF$Sol1[, 3]) * mmF$X[, 3] 

mVarF<- var(mFixed)

# MCMCglmm - marginal
mVarF/(mVarF+sum(apply(mmF$VCV1,2,mean)))
# alternative with crebile intervals
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(mmF$Sol1[i,] %*% t(mmF$X)))
  vmVarF[i]<-Var}

# getting R2m
(mVarF)/(mVarF+sum(apply(mmF$VCV1,2,mean)[-3])+ pi^2/3)
R2m<-vmVarF/(vmVarF+mmF$VCV1[,1]+mmF$VCV1[,2] + pi^2/3)
mean(R2m)
posterior.mode(R2m)
HPDinterval(R2m)

# MCMCglmm - conditional
(mVarF+sum(apply(mmF$VCV1,2,mean)[-3]))/(mVarF+sum(apply(mmF$VCV1,2,mean)[-3])+ pi^2/3)
# alternative with crebile intervals
R2c<-(vmVarF+mmF$VCV1[,1]+mmF$VCV1[,2])/(vmVarF+mmF$VCV1[,1]+mmF$VCV1[,2]+ pi^2/3)
mean(R2c)
posterior.mode(R2c)
HPDinterval(R2c)





#Colony size
prior1 <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)

mm0 <- MCMCglmm(log10(colony.size.y) ~ 1, random = ~animal, 
                data = ant_data2, family="gaussian", pedigree = tree_set[[1]], nitt = 110000, burnin = 10000, thin = 100, prior = prior1, verbose = F)
summary(mm0)

mmF <- MCMCglmm(log10(colony.size.y)~ TMPavg_avg+ sqrt(PREavg_avg)+ log(DTRavg_avg),
                data = ant_data2, random = ~animal, family = "gaussian", pedigree = tree_set[[1]], nitt = 110000, burnin = 10000, thin = 100, prior = prior1, verbose = F)
summary(mmF)

mFixed <- mean(mmF$Sol[,2]) * mmF$X[, 2] + mean(mmF$Sol[, 3]) * mmF$X[, 3] + mean(mmF$Sol[ ,4]) * mmF$X[, 4]

mVarF<- var(mFixed)

# MCMCglmm - marginal
mVarF/(mVarF+sum(apply(mmF$VCV,2,mean)))
# alternative with crebile intervals
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(mmF$Sol[i,] %*% t(mmF$X)))
  vmVarF[i]<-Var}

# getting R2m
R2m<-vmVarF/(vmVarF+mmF$VCV[,1])
mean(R2m)
posterior.mode(R2m)
HPDinterval(R2m)

# MCMCglmm - conditional
(mVarF+sum(apply(mmF$VCV,2,mean)[-2]))/(mVarF+sum(apply(mmF$VCV,2,mean)))
# alternative with crebile intervals
R2c<-(vmVarF+mmF$VCV[,1])/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2])

mean(R2c)
posterior.mode(R2c)
HPDinterval(R2c)

