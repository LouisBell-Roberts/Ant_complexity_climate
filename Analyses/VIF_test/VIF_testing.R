########AntEnv vif test for CS response - new transformations on predictors######
###Identify how significant the problem of multicollinearity is
##Data and trees have already been pruned
#Louis Bell-Roberts
#11/08/2023


# Loading essential packages ####
library(car)
library(ape)
library(phylolm)
library(phytools)

ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/AntClm_CSize_all_230710.csv") #After species with unusual biology removed = 474 species

ant_data_tropical<- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/Tropical_ant_data.csv")
ant_data_temperate<- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/Temperate_ant_data_230804.csv")
ant_data_both<- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/Tropical_vs_temperate/Both_ant_data_230804.csv")

#For global dataset
##Create model
full_model<- lm(log10(colony.size.y) ~ TMPavg_avg + sqrt(PREavg_avg) + log(DTRavg_avg) + sqrt(TMPsd_avg) + log(PREsd_avg), 
                     data=ant_data)
vif(full_model)

noTMPsd<- lm(log10(colony.size.y) ~ TMPavg_avg + sqrt(PREavg_avg) + log(DTRavg_avg) + log(PREsd_avg), 
                data=ant_data)
vif(noTMPsd)

noPREsd<- lm(log10(colony.size.y) ~ TMPavg_avg + sqrt(PREavg_avg) + log(DTRavg_avg) + sqrt(TMPsd_avg), 
                data=ant_data)
vif(noPREsd)

# avg_only<- lm(log10(colony.size.y) ~ TMPavg_avg + sqrt(PREavg_avg) + log(DTRavg_avg), 
#                 data=ant_data)
avg_only<- lm(log10(colony.size.y) ~ TMPavg_avg + PREavg_avg + DTRavg_avg, 
              data=ant_data) #Variables should not be transformed
vif(avg_only)

sd_only<- lm(log10(colony.size.y) ~ log(DTRavg_avg) + sqrt(TMPsd_avg) + log(PREsd_avg), 
                data=ant_data)
vif(sd_only)

#Thoughts: VIF values are not crazily high even in the full model (max = 5.9). However, high VIF scores seem to be driven by variables of the same climatic type simultaneously being present e.g. TMPavg and TMPsd or PREavg and PREsd. Including DTR doesn't seem to influence VIF scores. If avg or sd is removed, VIF score decrease dramatically. Therefore, should consider doing an analysis with either only avg or sd variables.



#For regional datasets
##Tropical
###Create model
full_model_tropical<- lm(log10(colony.size.y) ~ TMPavg_avg + sqrt(PREavg_avg) + log(DTRavg_avg) + sqrt(TMPsd_avg) + log(PREsd_avg), 
                data=ant_data_tropical)
vif(full_model_tropical)

# avg_only_tropical<- lm(log10(colony.size.y) ~ TMPavg_avg + sqrt(PREavg_avg) + log(DTRavg_avg), 
#               data=ant_data_tropical)
avg_only_tropical<- lm(log10(colony.size.y) ~ TMPavg_avg + PREavg_avg + DTRavg_avg, 
                       data=ant_data_tropical) #Variables should not be transformed
vif(avg_only_tropical)

sd_only_tropical<- lm(log10(colony.size.y) ~ log(DTRavg_avg) + sqrt(TMPsd_avg) + log(PREsd_avg), 
             data=ant_data_tropical)
vif(sd_only_tropical)


##Temperate
###Create model
full_model_temperate<- lm(log10(colony.size.y) ~ TMPavg_avg + sqrt(PREavg_avg) + log(DTRavg_avg) + sqrt(TMPsd_avg) + log(PREsd_avg), 
                         data=ant_data_temperate)
vif(full_model_temperate)

# avg_only_temperate<- lm(log10(colony.size.y) ~ TMPavg_avg + sqrt(PREavg_avg) + log(DTRavg_avg), 
#                        data=ant_data_temperate)
avg_only_temperate<- lm(log10(colony.size.y) ~ TMPavg_avg + PREavg_avg + DTRavg_avg, 
                        data=ant_data_temperate) #Variables should not be transformed
vif(avg_only_temperate)

sd_only_temperate<- lm(log10(colony.size.y) ~ log(DTRavg_avg) + sqrt(TMPsd_avg) + log(PREsd_avg), 
                      data=ant_data_temperate)
vif(sd_only_temperate)


##Both
###Create model
full_model_both<- lm(log10(colony.size.y) ~ TMPavg_avg + sqrt(PREavg_avg) + log(DTRavg_avg) + sqrt(TMPsd_avg) + log(PREsd_avg), 
                          data=ant_data_both)
vif(full_model_both)

# avg_only_both<- lm(log10(colony.size.y) ~ TMPavg_avg + sqrt(PREavg_avg) + log(DTRavg_avg), 
#                         data=ant_data_both)
avg_only_both<- lm(log10(colony.size.y) ~ TMPavg_avg + PREavg_avg + DTRavg_avg, 
                   data=ant_data_both) #Variables should not be transformed
vif(avg_only_both)

sd_only_both<- lm(log10(colony.size.y) ~ log(DTRavg_avg) + sqrt(TMPsd_avg) + log(PREsd_avg), 
                       data=ant_data_both)
vif(sd_only_both)


#########
#Create csv tables for vif scores

vif_tab_global <- vif(avg_only)
vif_tab_trop <- vif(avg_only_tropical)
vif_tab_temp <- vif(avg_only_temperate)
vif_tab_both <- vif(avg_only_both)
region <- c("Global", "Tropical", "Temperate", "Both")

vif_tab <- data.frame(t(data.frame(Global = vif_tab_global, Tropical = vif_tab_trop, Temperate = vif_tab_temp, Both = vif_tab_both)))
vif_tab_round <- round(vif_tab, digits = 2)
write.csv(vif_tab_round, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/AntEnv/Figures/VIF/vif_table.csv")



