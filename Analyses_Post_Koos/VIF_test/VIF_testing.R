########AntEnv vif test for CS response - new transformations on predictors######
###Identify how significant the problem of multicollinearity is
##Data and trees have already been pruned
#Louis Bell-Roberts
#13/11/2024


# Loading essential packages ####
library(car)
library(ape)
library(phylolm)
library(phytools)

ant_data <- read.csv("/Users/louis.bell-roberts/Documents/Github/Ant_complexity_climate/Supplementary_data/Data/AntTrait_Global_241027.csv") #After species with unusual biology removed = 501 species

ant_data_tropical<- read.csv("/Users/louis.bell-roberts/Documents/Github/Ant_complexity_climate/Supplementary_data/Data/AntTrait_Tropical_241103.csv")
ant_data_temperate<- read.csv("/Users/louis.bell-roberts/Documents/Github/Ant_complexity_climate/Supplementary_data/Data/AntTrait_Temperate_241103.csv")
ant_data_both<- read.csv("/Users/louis.bell-roberts/Documents/Github/Ant_complexity_climate/Supplementary_data/Data/AntTrait_Both_241103.csv")

#For global dataset
##Create model
avg_only<- lm(log10(colony.size) ~ TMP_avg + PRE_avg + DTR_avg, 
              data=ant_data) #Variables should not be transformed
vif(avg_only)


#For regional datasets
##Tropical
###Create model
avg_only_tropical<- lm(log10(colony.size) ~ TMP_avg + PRE_avg + DTR_avg, 
                       data=ant_data_tropical) #Variables should not be transformed
vif(avg_only_tropical)


##Temperate
###Create model
avg_only_temperate<- lm(log10(colony.size) ~ TMP_avg + PRE_avg + DTR_avg, 
                        data=ant_data_temperate) #Variables should not be transformed
vif(avg_only_temperate)


##Both
###Create model
avg_only_both<- lm(log10(colony.size) ~ TMP_avg + PRE_avg + DTR_avg, 
                   data=ant_data_both) #Variables should not be transformed
vif(avg_only_both)



#########
#Create csv tables for vif scores

vif_tab_global <- vif(avg_only)
vif_tab_trop <- vif(avg_only_tropical)
vif_tab_temp <- vif(avg_only_temperate)
vif_tab_both <- vif(avg_only_both)
region <- c("Global", "Tropical", "Temperate", "Both")

vif_tab <- data.frame(t(data.frame(Global = vif_tab_global, Tropical = vif_tab_trop, Temperate = vif_tab_temp, Both = vif_tab_both)))
vif_tab_round <- round(vif_tab, digits = 2)
# write.csv(vif_tab_round, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Figures/VIF/vif_table.csv")



