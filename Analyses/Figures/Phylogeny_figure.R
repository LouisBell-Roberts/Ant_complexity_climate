########################ASR with corHMM - multi-tree analysis; single run######################
###Designed to run on West server
##Runs ASR with data for just Caste number
#Louis Bell-Roberts
#21/04/2023

#Phylogenetic path analysis analysing variation in worker size
#Louis Bell-Roberts
#15/11/2023


library(tidyverse)
library(ape)
library(phylolm)
library(phytools)
library(ggplot2)
library(phylopath)
library(car)

#Read in data file
ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Following_review/ant_data.csv")

#Set variables so that they're in the correct structure and apply transformations
ant_data[ant_data == ""] <- NA #Replace blank by NA
class(ant_data$colony.size) #numeric
class(ant_data$queen.mating.frequency) #numeric
class(ant_data$queen.number.continuous) #numeric
class(ant_data$worker.size.variation) #numeric
class(ant_data$caste.number) #numeric
ant_data$queen.number.binary <- as.factor(ant_data$queen.number.binary)

ant_data$colony.size <- log10(ant_data$colony.size)
ant_data$queen.mating.frequency <- log10(ant_data$queen.mating.frequency)
ant_data$queen.number.continuous <- log10(ant_data$queen.number.continuous)
ant_data$worker.size.variation <- sqrt(ant_data$worker.size.variation)

#Set rownames as species names for phylolm package
rownames(ant_data) <- ant_data$species

#Read in phylogenetic trees
NCuniform_stem <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15k_NCuniform_stem_mcc.tre")
NCuniform_crown <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15K_NCuniform_crown_mcc.tre")
FBD_stem <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15K_FBD_stem_mcc.tre")
FBD_crown <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15K_FBD_crown_mcc.tre")

#Filter data
all_variables <- dplyr::filter(ant_data, complete.cases(worker.size.variation), complete.cases(queen.mating.frequency), complete.cases(colony.size), complete.cases(queen.number.continuous))

#Prune tree
NCuniform_stem_pruned <- drop.tip(NCuniform_stem, setdiff(NCuniform_stem$tip.label, all_variables$species))
NCuniform_crown_pruned <- drop.tip(NCuniform_crown, setdiff(NCuniform_crown$tip.label, all_variables$species))
FBD_stem_pruned <- drop.tip(FBD_stem, setdiff(FBD_stem$tip.label, all_variables$species))
FBD_crown_pruned <- drop.tip(FBD_crown, setdiff(FBD_crown$tip.label, all_variables$species))

#Prune database
all_variables <- filter(all_variables, species %in% NCuniform_stem_pruned$tip.label)


##########
#Plotting transitions
##########

#Transitions can be plotted using different colouring schemes
##Method 1 - transitions plotted for every type of increase in castes, with colour based on the caste number that results e.g. if a species ends up with 3 castes (regardless of the number of castes that it started with) it's colour yellow
##Method 2 - transition plotted for every type of increase in caste, with a different colour for each type of jump



##
#Method 1 - Subfamily labels without common names. Silhouettes with common names attaches dotted around the outside of the figure
##

#Change caste number column name
sxData_rename <- sxData %>% dplyr::rename(`Number of worker castes` = sxData.Caste3)
sxData_rename$`Number of worker castes` <- as.factor(sxData_rename$`Number of worker castes`)

#Add node labels to phylogeny - nodes start at 570 and there are 568 of them
ASR_models[[1]]$phy$node.label <- c(seq(570,(570+567)))

nicer_tree <- ggtree(ASR_models[[1]]$phy, layout="fan", branch.length = "none", open.angle = 25, size = 0.1) + geom_rootedge(TRUE, size = 0.1)# + geom_tiplab(size = 0.9, offset = 2) + geom_nodelab(size = 0.5)

#Useful code for plotting tip and node labels
# nicer_tree <- ggtree(ASR_models[[1]]$phy, layout="fan", branch.length = "none", open.angle = 25, size = 0.1) + geom_rootedge(TRUE, size = 0.1) + geom_tiplab(size = 0.9, offset = 2) + geom_nodelab(size = 0.5)

nicer_tree_rotated <- rotate_tree(nicer_tree, angle=90)

tree_styled <- nicer_tree_rotated +
  geom_fruit(data=sxData_rename, # Data
             geom=geom_tile, # Plots 'tiles' (squares) onto each tip
             width=1.6, #Alters the length of the tiles
             position=position_identityx(hexpand=30), # Adjusts how far away the tiles are from the tre
             mapping=aes(y=sxData.animal, fill=`Number of worker castes`), # Analogous to ggplot aes()
             color = "white", # Colour of outline round the tiles (white makes it look like a gap)
             lwd = 0.2, # Width of line between tiles
             linetype = 1, # Default - other numbers make the line dashes, dotted etc.
             axis.params=list( # Add label to the geom_tile - by adding an x-axis
               axis="x",
               text = " ", # Label to plot
               text.size = 3.5, # Size of text
               hjust = 0, # Adjust position of text relative to the geom_tile
               vjust = 0.8,
               text.angle=360,
               line.colour="white"), # Set to white so axis line is not visible
             offset = 4, # Fine-scale adjustment of space between tree & layers
             pwidth=10) + # Width of whole plot
  scale_fill_manual(values=c("#E5F5FF", "#A6CFFA", "#5278B7", "#000005"))
tree_styled


#Add trasitions from 1 to 2
tree_trans <- tree_styled + ggtree::geom_point2(aes(subset = node %in% one_to_two_trans),
                                                size = 5, colour = '#A6CFFA', shape = 20, alpha = 0.85) #When plotting for multiple nodes, %in% is crucial. we assigning node as == to something. Therefore, to make node == to multiple values, we need the %in% operator.

#Add transitions from 1 to 3
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% one_to_three_trans),
                                               size = 5, colour = '#5278B7', shape = 20, alpha = 0.75)

#Add transitions from 1 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% one_to_four_trans),
                                               size = 5, colour = '#000005', shape = 20, alpha = 0.6)

#Add transitions from 2 to 3
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% two_to_three_trans),
                                               size = 5, colour = '#5278B7', shape = 20, alpha = 0.75)

#Add transitions from 2 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% two_to_four_trans),
                                               size = 5, colour = '#000005', shape = 20, alpha = 0.6)
#Add transitions from 3 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% three_to_four_trans),
                                               size = 5, colour = '#000005', shape = 20, alpha = 0.6)
tree_trans

tree_clad_lab <- tree_trans + geom_cladelab(node=1122, label="Amblyoponinae\n & Apomyrminae", align=TRUE, fontsize = 3.5, angle="auto",
                                            offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=443, label="Anaeuretinae", align=TRUE, fontsize = 3.5, angle="auto",
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=557, label="Apomyrminae", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black') # this subfamily is contained within Amblyoponinae

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=983, label="Dolichoderinae", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1038, label="Dorylinae                   ", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=841, label="Ectatomminae", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=854, label="Formicinae", align=TRUE, fontsize = 3.5, angle=240,
                                               offset = 2, offset.text = 2, textcolor='black', barcolor='black', hjust = 9, vjust = 9)

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=569, label="Leptanillinae", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1031, label="Myrmeciinae", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=576, label="Myrmicinae", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=550, label="Paraponerinae", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1061, label="Ponerinae", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')


tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1132, label="Proceratiinae", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1017, label="Pseudomyrmecinae", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=671, label="Big-headed ants", align=TRUE, fontsize = 3.5, angle="auto",
                                               offset = 3, offset.text = 0.5, textcolor='#777777', barcolor='#777777')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=719, label="Turtle ants", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 3, offset.text = 0.5, textcolor='#777777', barcolor='#777777')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=728, label="Fungus-growing \nants (including \nAtta and\nAcromyrmex)", align=TRUE, fontsize = 3.5, angle="auto",
                                               offset = 3, offset.text = 0.5, textcolor='#777777', barcolor='#777777')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=907, label="Wood ants", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 3, offset.text = 0.5, textcolor='#777777', barcolor='#777777')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=941, label="Desert ants", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 3, offset.text = 0.5, textcolor='#777777', barcolor='#777777')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=861, label="Carpenter ants", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 3, offset.text = 0.5, textcolor='#777777', barcolor='#777777')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=767, label="Fire ants", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 3, offset.text = 0.5, textcolor='#777777', barcolor='#777777')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1038, label="(army ants)", align=TRUE, fontsize = 3.5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='#777777', barcolor='black', barsize = 0)

tree_clad_lab

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/Version_9/Figures/")
# ggsave("Phylogeny_common_names.pdf", plot = tree_clad_lab, width = 11, height = 9, units = "in", dpi = 300)


###
#Method 2 - subfamily labels with common names in brackets
###
# 
# #Change caste number column name
# sxData_rename <- sxData %>% dplyr::rename(`Number of worker castes` = sxData.Caste3)
# sxData_rename$`Number of worker castes` <- as.factor(sxData_rename$`Number of worker castes`)
# 
# #Add node labels to phylogeny - nodes start at 570 and there are 568 of them
# ASR_models[[1]]$phy$node.label <- c(seq(570,(570+567)))
# 
# nicer_tree <- ggtree(ASR_models[[1]]$phy, layout="fan", branch.length = "none", open.angle = 25, size = 0.1) + geom_rootedge(TRUE, size = 0.1)# + geom_tiplab(size = 0.9, offset = 2) + geom_nodelab(size = 0.5)
# 
# #Useful code for plotting tip and node labels
# # nicer_tree <- ggtree(ASR_models[[1]]$phy, layout="fan", branch.length = "none", open.angle = 25, size = 0.1) + geom_rootedge(TRUE, size = 0.1) + geom_tiplab(size = 0.9, offset = 2) + geom_nodelab(size = 0.5)
# 
# nicer_tree_rotated <- rotate_tree(nicer_tree, angle=90)
# 
# tree_styled <- nicer_tree_rotated +
#   geom_fruit(data=sxData_rename, # Data
#              geom=geom_tile, # Plots 'tiles' (squares) onto each tip
#              width=1.6, #Alters the length of the tiles
#              position=position_identityx(hexpand=30), # Adjusts how far away the tiles are from the tre
#              mapping=aes(y=sxData.animal, fill=`Number of worker castes`), # Analogous to ggplot aes()
#              color = "white", # Colour of outline round the tiles (white makes it look like a gap)
#              lwd = 0.2, # Width of line between tiles
#              linetype = 1, # Default - other numbers make the line dashes, dotted etc.
#              axis.params=list( # Add label to the geom_tile - by adding an x-axis
#                axis="x",
#                text = " ", # Label to plot
#                text.size = 3.5, # Size of text
#                hjust = 0, # Adjust position of text relative to the geom_tile
#                vjust = 0.8,
#                text.angle=360,
#                line.colour="white"), # Set to white so axis line is not visible
#              offset = 4, # Fine-scale adjustment of space between tree & layers
#              pwidth=10) + # Width of whole plot
#   scale_fill_manual(values=c("#E5F5FF", "#A6CFFA", "#5278B7", "#000005"))
# tree_styled
# 
# #Add trasitions from 1 to 2
# tree_trans <- tree_styled + ggtree::geom_point2(aes(subset = node %in% one_to_two_trans),
#                                                 size = 5, colour = '#A6CFFA', shape = 20, alpha = 0.85) #When plotting for multiple nodes, %in% is crucial. we assigning node as == to something. Therefore, to make node == to multiple values, we need the %in% operator.
# 
# #Add transitions from 1 to 3
# tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% one_to_three_trans),
#                                                size = 5, colour = '#5278B7', shape = 20, alpha = 0.75)
# 
# #Add transitions from 1 to 4
# tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% one_to_four_trans),
#                                                size = 5, colour = '#000005', shape = 20, alpha = 0.6)
# 
# #Add transitions from 2 to 3
# tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% two_to_three_trans),
#                                                size = 5, colour = '#5278B7', shape = 20, alpha = 0.75)
# 
# #Add transitions from 2 to 4
# tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% two_to_four_trans),
#                                                size = 5, colour = '#000005', shape = 20, alpha = 0.6)
# #Add transitions from 3 to 4
# tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% three_to_four_trans),
#                                                size = 5, colour = '#000005', shape = 20, alpha = 0.6)
# tree_trans
# 
# tree_clad_lab <- tree_trans + geom_cladelab(node=1122, label="Amblyoponinae\n & Apomyrminae", align=TRUE, fontsize = 3, angle="auto",
#                                             offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=443, label="Anaeuretinae", align=TRUE, fontsize = 3, angle="auto",
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# # tree_clad_lab <- tree_clad_lab + geom_cladelab(node=557, label="Apomyrminae", align=TRUE, fontsize = 3, angle='auto',
# #                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black') # this subfamily is contained within Amblyoponinae
# 
# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=983, label="Dolichoderinae", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1038, label="Dorylinae\n (army ants)", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=841, label="Ectatomminae", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=854, label="Formicinae\n (including the wood, \ncarpenter, weaver\n and desert ants)", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=569, label="Leptanillinae", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1031, label="Myrmeciinae", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=576, label="Myrmicinae \n(including the \nbig-headed, \nturtle, fire and \nleaf-cutting \nants)", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=550, label="Paraponerinae", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1061, label="Ponerinae", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# 
# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1132, label="Proceratiinae", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1017, label="Pseudomyrmecinae", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')
# 
# tree_clad_lab




###
#Method 3 - with Anna's plotting style. Circular
###

#Change caste number column name
sxData_rename <- sxData %>% dplyr::rename(`Number of worker castes` = sxData.Caste3)
sxData_rename$`Number of worker castes` <- as.factor(sxData_rename$`Number of worker castes`)

#Add node labels to phylogeny - nodes start at 570 and there are 568 of them
ASR_models[[1]]$phy$node.label <- c(seq(570,(570+567)))

nicer_tree <- ggtree(ASR_models[[1]]$phy, layout="fan", branch.length = "none", open.angle = 25, size = 0.1) + geom_rootedge(TRUE, size = 0.1)# + geom_tiplab(size = 0.9, offset = 2) + geom_nodelab(size = 0.5)

nicer_tree_rotated <- rotate_tree(nicer_tree, angle=90)

tree_styled <- nicer_tree_rotated +
  geom_fruit(data=sxData_rename, # Data
             geom=geom_tile, # Plots 'tiles' (squares) onto each tip
             width=1.6, #Alters the length of the tiles
             position=position_identityx(hexpand=30), # Adjusts how far away the tiles are from the tre
             mapping=aes(y=sxData.animal, fill=`Number of worker castes`), # Analogous to ggplot aes()
             color = "white", # Colour of outline round the tiles (white makes it look like a gap)
             lwd = 0.2, # Width of line between tiles
             linetype = 1, # Default - other numbers make the line dashes, dotted etc.
             axis.params=list( # Add label to the geom_tile - by adding an x-axis
               axis="x",
               text = " ", # Label to plot
               text.size = 3.5, # Size of text
               hjust = 0, # Adjust position of text relative to the geom_tile
               vjust = 0.8,
               text.angle=360,
               line.colour="white"), # Set to white so axis line is not visible
             offset = 4, # Fine-scale adjustment of space between tree & layers
             pwidth=100) + # Width of whole plot
  scale_fill_manual(values=c("#E5F5FF", "#A6CFFA", "#5278B7", "#000005")) +
  theme_tree(legend.position = "none") #This line removes the legend

tree_styled

#Add trasitions from 1 to 2
tree_trans <- tree_styled + ggtree::geom_point2(aes(subset = node %in% one_to_two_trans),
                                                size = 5, colour = '#A6CFFA', shape = 20, alpha = 0.85) #When plotting for multiple nodes, %in% is crucial. we assigning node as == to something. Therefore, to make node == to multiple values, we need the %in% operator.

#Add transitions from 1 to 3
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% one_to_three_trans),
                                               size = 5, colour = '#5278B7', shape = 20, alpha = 0.75)

#Add transitions from 1 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% one_to_four_trans),
                                               size = 5, colour = '#000005', shape = 20, alpha = 0.6)

#Add transitions from 2 to 3
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% two_to_three_trans),
                                               size = 5, colour = '#5278B7', shape = 20, alpha = 0.75)

#Add transitions from 2 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% two_to_four_trans),
                                               size = 5, colour = '#000005', shape = 20, alpha = 0.6)
#Add transitions from 3 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% three_to_four_trans),
                                               size = 5, colour = '#000005', shape = 20, alpha = 0.6)
tree_trans

tree_clad_lab <- tree_trans + geom_cladelab(node=671, label="Pheidole\n(big-headed ants)", align=TRUE, fontsize = 5, angle="auto",
                                            offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=728, label="Fungus-growing \nants (including \nAtta and\nAcromyrmex)", align=TRUE, fontsize = 5, angle="auto",
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

# tree_clad_lab <- tree_trans + geom_cladelab(node=666, label="Attini (fungus-growing\n ants)", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=606, label="Carebara", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=861, label="Camponotus &\nCalomyrmex \n(including \ncarpenter ants)", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=941, label="Cataglyphis\n(desert ants)", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=907, label="Formica (wood ants)", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=962, label="Lasius &\nMyrmecocystus", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=804, label="Pogonomyrmex", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=632, label="Temnothorax", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=586, label="Crematogaster", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1062, label="Ponerinae", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1038, label="Dorylinae (army ants)", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')


tree_clad_lab <- tree_clad_lab + geom_cladelab(node=983, label="Dolichoderinae", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=766, label="Solenopsidini \n(including fire ants)", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=773, label="Stenammini", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=825, label="Myrmicini", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=842, label="Ectatommini", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=889, label="Polyrhachis", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=978, label="Brachymyrmex", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1017, label="Pseudomyrmecini", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1124, label="Amblyoponini", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=603, label="Acanthomyrmex", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=617, label="Tetramorium", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=719, label="Cephalotes \n(turtle ants)", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/Version_9/Figures/")
# ggsave("Phylogeny_multi_tax_blue_no_leg.pdf", plot = tree_clad_lab, width = 13, height = 13, units = "in", dpi = 300)








###
#Method 4 - same plot as method 3, but without the subfamily labels so that they can be added in powerpoint
###

#Change caste number column name
sxData_rename <- sxData %>% dplyr::rename(`Number of worker castes` = sxData.Caste3)
sxData_rename$`Number of worker castes` <- as.factor(sxData_rename$`Number of worker castes`)

#Add node labels to phylogeny - nodes start at 570 and there are 568 of them
ASR_models[[1]]$phy$node.label <- c(seq(570,(570+567)))

nicer_tree <- ggtree(ASR_models[[1]]$phy, layout="fan", branch.length = "none", open.angle = 25, size = 0.1) + geom_rootedge(TRUE, size = 0.1)# + geom_tiplab(size = 0.9, offset = 2) + geom_nodelab(size = 0.5)

nicer_tree_rotated <- rotate_tree(nicer_tree, angle=90)

tree_styled <- nicer_tree_rotated +
  geom_fruit(data=sxData_rename, # Data
             geom=geom_tile, # Plots 'tiles' (squares) onto each tip
             width=1.6, #Alters the length of the tiles
             position=position_identityx(hexpand=30), # Adjusts how far away the tiles are from the tre
             mapping=aes(y=sxData.animal, fill=`Number of worker castes`), # Analogous to ggplot aes()
             color = "white", # Colour of outline round the tiles (white makes it look like a gap)
             lwd = 0.2, # Width of line between tiles
             linetype = 1, # Default - other numbers make the line dashes, dotted etc.
             axis.params=list( # Add label to the geom_tile - by adding an x-axis
               axis="x",
               text = " ", # Label to plot
               text.size = 3.5, # Size of text
               hjust = 0, # Adjust position of text relative to the geom_tile
               vjust = 0.8,
               text.angle=360,
               line.colour="white"), # Set to white so axis line is not visible
             offset = 4, # Fine-scale adjustment of space between tree & layers
             pwidth=100) + # Width of whole plot
  scale_fill_manual(values=c("#E5F5FF", "#A6CFFA", "#5278B7", "#000005")) +
  theme_tree(legend.position = "none") #This line removes the legend

tree_styled

#Add trasitions from 1 to 2
tree_trans <- tree_styled + ggtree::geom_point2(aes(subset = node %in% one_to_two_trans),
                                                size = 5, colour = '#A6CFFA', shape = 20, alpha = 0.85) #When plotting for multiple nodes, %in% is crucial. we assigning node as == to something. Therefore, to make node == to multiple values, we need the %in% operator.

#Add transitions from 1 to 3
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% one_to_three_trans),
                                               size = 5, colour = '#5278B7', shape = 20, alpha = 0.75)

#Add transitions from 1 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% one_to_four_trans),
                                               size = 5, colour = '#000005', shape = 20, alpha = 0.6)

#Add transitions from 2 to 3
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% two_to_three_trans),
                                               size = 5, colour = '#5278B7', shape = 20, alpha = 0.75)

#Add transitions from 2 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% two_to_four_trans),
                                               size = 5, colour = '#000005', shape = 20, alpha = 0.6)
#Add transitions from 3 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% three_to_four_trans),
                                               size = 5, colour = '#000005', shape = 20, alpha = 0.6)
tree_trans

tree_clad_lab <- tree_trans + geom_cladelab(node=671, label="", align=TRUE, fontsize = 5, angle="auto",
                                            offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=728, label="", align=TRUE, fontsize = 5, angle="auto",
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

# tree_clad_lab <- tree_trans + geom_cladelab(node=666, label="Attini (fungus-growing\n ants)", align=TRUE, fontsize = 3, angle='auto',
#                                                offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=606, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=861, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=941, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=907, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=962, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=804, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=632, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=586, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1062, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1038, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')


tree_clad_lab <- tree_clad_lab + geom_cladelab(node=983, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=766, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=773, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=825, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=842, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=889, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=978, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1017, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1124, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=603, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=617, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=719, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_no_lab <- tree_clad_lab

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Wordfiles/Paper drafting/Version_9/Figures/")
ggsave("Phylogeny_multi_tax_blue_no_leg_no_lab.pdf", plot = tree_clad_no_lab, width = 8, height = 8, units = "in", dpi = 300)



####################
#Estimating the ancestral number of worker castes at the root of the phylogeny across 400 trees - run this code after line 200
# int_list <- list()
# 
# for (i in seq_along(ASR_states_list)) {
#   ASR_states <- ASR_states_list[[i]]
#   
#   # Combine the probability scores for the two different rate classes for each caste number
#   ASR_states_summed <- data.frame(
#     caste_1 = ASR_states$`(1,R1)` + ASR_states$`(1,R2)`,
#     caste_2 = ASR_states$`(2,R1)` + ASR_states$`(2,R2)`,
#     caste_3 = ASR_states$`(3,R1)` + ASR_states$`(3,R2)`,
#     caste_4 = ASR_states$`(4,R1)` + ASR_states$`(4,R2)`
#   )
#   
#   int_list[[i]] <- ASR_states_summed
# }
# 
# #Create a vector that is all of the probability estimates across the 400 trees that the MRCA for all ants had a single worker caste
# values_vector <- sapply(int_list, function(df) df[1, 1])
# mean(values_vector) #0.991207%