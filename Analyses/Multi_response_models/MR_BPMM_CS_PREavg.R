########################Calculate measure of phylogenetic correlation between ant climatic and colony size varaibles ######################
###Calculate over initially 1 tree and then 400 trees using MCMCglmm
##Run MCMCglmm models in parallel
#Louis Bell-Roberts
#13/08/2023

#**************************************************************************************
#Packages ####
#**************************************************************************************
pacman::p_load("MCMCglmm")

#Load and prepare data and phylogenies
pre_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/AntClm_CSize_all_230710.csv") #After species with unusual biology removed = 474 species

#Read in ant trees
ant_trees <-read.tree("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/400_trees/400trees_pruned_CSize_230712.tre")
ant_tree <- ant_trees[[301]]

pre_data$valid_species_name <- gsub("_", ".", pre_data$valid_species_name, fixed = T)
is.ultrametric(ant_trees) # true
is.binary(ant_trees)

#Rename valid_species_name column
pre_data <- pre_data %>% 
  rename(animal=valid_species_name,
         PREavg=PREavg_avg,
         TMPsd=TMPsd_avg,
         colony.size=colony.size.y)

#Select data
data_select <- dplyr::select(pre_data, animal, PREavg, colony.size)

#**************************************************************************************
#Model settings ----
#**************************************************************************************
#To test if models run
nitt<-10 # number of iterations
burnin<-1 # burnin period
thin<-1 #thinning interval

#*******************
#Longer runs
nitt<-110000
burnin<-10000
thin<-100


#**************************************************************************************
#Single response Poisson model ----
#**************************************************************************************
#Need to create a matrix from the tree that will then be fed into the model
phyloA<-inverseA(ant_tree)$Ainv

# #Specify the priors for the random effects: 
# prior1 =list(R = list(V = diag(1),nu=0.002), #R = residual variation.
#              G = list(G1=list(V = diag(1),nu=0.002))) #G = between unit variance. In this example there is only one (G1) random effect linked to the phylogeny
# 
# M1<-MCMCglmm(trait ~ 1, #No fixed effects, just an intercept "1"
#              random=~species, 
#              ginverse=list(species=phyloA), #Links data with tree: all species in df need to be in the tree, but tips in tree can have missing data
#              family ="poisson", #specify distribution of data 
#              data = df, 
#              prior=prior1, 
#              nitt=nitt, burnin=burnin,thin=thin) 
# #Look at output
# plot(M1)
# summary(M1)
# 
# #Estimates of fixed effects
# posterior.mode(M1$Sol)
# #95% Credible intervals for fixed effects
# HPDinterval(M1$Sol)
# 
# #Estimates of random effects
# posterior.mode(M1$VCV)
# #95% Credible intervals for random effects
# HPDinterval(M1$VCV)

#**************************************************************************************
#Bivariate response model with one Poisson and one Gaussian trait ----
#**************************************************************************************
prior2 =list(R = list(V = diag(2),nu=1), #now need to specify a 2x2 matrix for priors as there is a variance for each response + their covariance. Matrix produced tells the model that the residual variance is independent between each trait. Print the diag(2) function and one will see that 1 to 1 = 1, and 2 to 2 = 1. i.e. assumes there is no covariance between the two traits/not considering covariance between two traits.
             #If I were to swap the 1s and 0s in the matrix, that would be saying, 'there would only be an effect when both traits are there'
             G = list(G1=list(V = diag(2), nu = 1))) #Same for phylogenetic effects - variance for each response + their covariance

M2<-MCMCglmm(cbind(log10(colony.size), sqrt(PREavg))~trait-1,
             random=~us(trait):animal, 
             ginverse=list(animal=phyloA),
             rcov = ~us(trait):units,
             family =c("gaussian","gaussian"), 
             data = data_select,prior=prior2, 
             nitt=nitt, burnin=burnin,thin=thin,verbose=F) 

#Look at output
plot(M2)
summary(M2)

#Estimates of fixed effects
posterior.mode(M2$Sol)
#95% Credible intervals for fixed effects
HPDinterval(M2$Sol)

#Estimates of random effects - you should now see 4 estimates for residual ('units') and 4 for species.
posterior.mode(M2$VCV)
#95% Credible intervals for random effects
HPDinterval(M2$VCV)

#These 4 estimates = variance for trait1, cov(trait1,trait2), cov(trait2,trait1),variance for trait2
#you can then calculate phylogenetic correlation between trait1 and trait2: cor = cov1,2/sqrt(var1*var2)


#**************************************************************************************
#Compute heritability and phylogenetic correlation values
#**************************************************************************************

#Heritability
herit_CS <- M2$VCV[ , "traitcolony.size:traitcolony.size.animal"] / (M2$VCV[ , "traitcolony.size:traitcolony.size.animal"] + M2$VCV[ , "traitcolony.size:traitcolony.size.units"])

herit_PREavg <- M2$VCV[ , "traitPREavg:traitPREavg.animal"] / (M2$VCV[ , "traitPREavg:traitPREavg.animal"] + M2$VCV[ , "traitPREavg:traitPREavg.units"])

posterior.mode(herit_CS)
HPDinterval(herit_CS)
posterior.mode(herit_PREavg)
HPDinterval(herit_PREavg)

#Phylogenetic correlation
genetic_corr <- M2$VCV[ , "traitcolony.size:traitPREavg.animal"] / sqrt(M2$VCV[ , "traitcolony.size:traitcolony.size.animal"] * M2$VCV[ , "traitPREavg:traitPREavg.animal"])

resid_corr <- M2$VCV[ , "traitcolony.size:traitPREavg.units"] / sqrt(M2$VCV[ , "traitcolony.size:traitcolony.size.units"] * M2$VCV[ , "traitPREavg:traitPREavg.units"])

posterior.mode(genetic_corr)
HPDinterval(genetic_corr)
posterior.mode(resid_corr)
HPDinterval(resid_corr)


#**************************************************************************************
#Checking convergence ----
#**************************************************************************************
#run model 3 times e.g.
# 
# M2_1<-MCMCglmm(cbind(trait1,trait2)~trait-1,
#                random=~us(trait):species, 
#                ginverse=list(species=phyloA),
#                rcov = ~us(trait):units,
#                family =c("poisson","gaussian"), 
#                data = df,prior=prior2, 
#                nitt=nitt, burnin=burnin,thin=thin,verbose=F) 
# 
# M2_2<-MCMCglmm(cbind(trait1,trait2)~trait-1,
#                random=~us(trait):species, 
#                ginverse=list(species=phyloA),
#                rcov = ~us(trait):units,
#                family =c("poisson","gaussian"), 
#                data = df,prior=prior2, 
#                nitt=nitt, burnin=burnin,thin=thin,verbose=F) 
# 
# M2_3<-MCMCglmm(cbind(trait1,trait2)~trait-1,
#                random=~us(trait):species, 
#                ginverse=list(species=phyloA),
#                rcov = ~us(trait):units,
#                family =c("poisson","gaussian"), 
#                data = df,prior=prior2, 
#                nitt=nitt, burnin=burnin,thin=thin,verbose=F) 
# 
# 
# #Fixed effects
# M2.Sol<-mcmc.list(M2_1$Sol,M2_2$Sol,M2_3$Sol)
# 
# #Traces should overlap
# plot(M2.Sol)
# 
# #Random effects
# M2.VCV<-mcmc.list(M2_1$VCV,M2_2$VCV,M2_3$VCV)
# plot(Grp_1.VCV)
# 
# #Tests of convergence - PSR should be less than 1.1
# gelman.diag(M2.Sol,multivariate = FALSE)
# gelman.diag(M2.VCV,multivariate = FALSE)
