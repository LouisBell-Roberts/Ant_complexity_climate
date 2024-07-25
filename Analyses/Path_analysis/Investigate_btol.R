## Investigating  btol parameter ##
#01/11/2023

# Packages
library(ape)
library(phylolm)
library(ggplot2)
library(tidyverse)
library(phylopath)
library(phytools)

#Read in ant data
ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/AntClm_CSize_all_230710.csv") #After species with unusual biology removed = 474 species

#Assign rownames as species names
rownames(ant_data) <- ant_data$valid_species_name

#Log transform colony size variable
ant_data$colony.size.y <- log10(ant_data$colony.size.y)

#Read in trees
anttree_NCuniform_stem <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/MCC_trees/pruned_474_sp_NCuniform_stem_MCC.tre")
anttree_NCuniform_crown <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/MCC_trees/pruned_474_sp_NCuniform_crown_MCC.tre")
anttree_FBD_stem <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/MCC_trees/pruned_474_sp_FBD_stem_MCC.tre")
anttree_FBD_crown <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Trees/CS/474_species/MCC_trees/pruned_474_sp_FBD_crown_MCC.tre")

#Select columns of interest
ant_data_temp_avg_selected <- ant_data %>% dplyr::select(valid_species_name, poly_id, colony.size.y, TMPavg_avg)
ant_data_rain_avg_selected <- ant_data %>% dplyr::select(valid_species_name, poly_id, colony.size.y, PREavg_avg)

#Rename the variables used in the analysis
ant_data_temp_avg_selected_named <- ant_data_temp_avg_selected %>% 
  rename(Colony_size=colony.size.y,
         Poly_id=poly_id,
         TMPavg=TMPavg_avg)

ant_data_rain_avg_selected_named <- ant_data_rain_avg_selected %>%
  rename(Colony_size=colony.size.y,
         Poly_id=poly_id,
         PREavg=PREavg_avg)

# PREavg

#Run the models on each of the 4 different MCC trees
test_model <- phyloglm(Poly_id ~ Colony_size + PREavg, data = ant_data_rain_avg_selected_named, phy = anttree_NCuniform_stem, method = "logistic_MPLE", btol = 66)
summary(test_model)
sort(test_model$fitted.values)

#Plot - PRE only
plot(ant_data_rain_avg_selected_named$PREavg,jitter(ant_data_rain_avg_selected_named$Poly_id,factor=0,amount=0.02),
     xlab="PRE",ylab="Poly_id",xlim=c())

cc <- coef(test_model)
curve(plogis(cc[1]+cc[2]*x),col="red",add=TRUE)



#Plot single predictor model -does it match the fitted values?
#Can fitted values not lie near to 0 or 1? - yes
#Check if changing btol affects parameter estimates at all
#Check if any of the phylopath models have warnings worse than btol reached
#Check fitted values of the phylopath models with warnings

max(test$fitted.values)
# par(las=1,bty="l") ## cosmetic
plot(ant_data_rain_avg_selected_named$PREavg,jitter(as.numeric(ant_data_rain_avg_selected_named$Poly_id),factor=0,amount=0.02),
     xlab="PREavg",ylab="Poly_id",xlim=c())
cc <- coef(fit)
curve(plogis(cc[1]+cc[2]*x),col="red",add=TRUE)
#plot my data and logistic regression using ben bolker code - do i have complete separation? does the logistic curve get bounded at 0 or 1?













# Set a random seed for reproducibility
set.seed(42)

# Number of observations
n <- 100

# Simulate the binary variable (0 and 1)
binary_variable <- sample(0:1, n, replace = TRUE)

# Simulate the continuous variable in the range [10, 20]
continuous_variable <- runif(n, min = 10, max = 20)

# Create a data frame
simulated_data <- data.frame(BinaryVariable = binary_variable, ContinuousVariable = continuous_variable)

# Display the first few rows of the simulated data
head(simulated_data)

# Load the ggplot2 package
library(ggplot2)

# Create a scatterplot
ggplot(simulated_data, aes(x = ContinuousVariable, y = BinaryVariable)) +
  geom_point(aes(color = factor(BinaryVariable))) +
  labs(x = "Continuous Variable", y = "Binary Variable") +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()

##############################
PREavg_result_NCuniform_stem <- phylo_path(PREavg_models, data = ant_data_rain_avg_selected_named, tree = anttree_NCuniform_stem, model = "lambda", method = "logistic_MPLE")
test <- phyloglm(Poly_id ~ Colony_size + PREavg, data = ant_data_rain_avg_selected_named, phy = anttree_NCuniform_stem, method = "logistic_MPLE", btol = 4)
max(test$fitted.values)
par(las=1,bty="l") ## cosmetic
plot(ant_data_rain_avg_selected_named$PREavg,jitter(as.numeric(ant_data_rain_avg_selected_named$Poly_id)-1,factor=0,amount=0.02),
     xlab="PREavg",ylab="Poly_id",xlim=c())
cc <- coef(test)
curve(plogis(cc[1]+cc[3]*x),col="red",add=TRUE)
