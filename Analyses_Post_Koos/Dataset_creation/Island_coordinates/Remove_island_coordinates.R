########### Comparing the Geographic occurrence data between La Rich and Kass et al., 2022 ########### 
#Louis Bell-Roberts
#14/10/2024

library(rworldmap)
library(rgeos)
library(maptools)
library(cleangeo)  ## For clgeo_Clean()
library(raster)
library(sp)
library(sf)
library(rgdal)
library(tools)
library(maps)

#Aims
##Remove island coordinates from Kass
#Read in csv file - already has all duplicate rows removed
Kass_orgin <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Kass_2022/3A_forAnalysis_processed_database_SPECIES.csv")
length(Kass_orgin$lon_opt) #1479293

#############################################
### Remove islands ###
#############################################
#Rename columns
coords<-as.data.frame(cbind(Kass_orgin$lon_opt,Kass_orgin$lat_opt))
length(coords$V1) #1479293
colnames(coords)<-c("Long", "Lat")
coords<-na.omit(coords) #rows with NA have already been omitted
length(coords$Long) #1479293

#convert coordinates into spatial object
coords<-SpatialPoints(coords=coords,proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
st_crs(coords)
# http://rgdal.r-forge.r-project.org/articles/PROJ6_GDAL3.html
cont<-rgdal::readOGR("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/La_Richeliere_dryad/50m_cultural/ne_50m_admin_0_countries.shp") 

proj4string(coords)<- proj4string(cont) #What does this line of code do?

##read in the downloaded ocean file as a shape file
oceans <-rgdal::readOGR("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/La_Richeliere_dryad/ne_110m_ocean/ne_110m_ocean.shp", "ne_110m_ocean",stringsAsFactors=FALSE, verbose=T)

### Download ocean layer - done this already
# URL <- "http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_ocean.zip"
# fil <- basename(URL)
# if (!file.exists(fil)) download.file(URL, fil)
# fils <- unzip(fil)
# oceans <- rgdal::readOGR(grep("shp$", fils, value=TRUE), "ne_110m_ocean",
# stringsAsFactors=FALSE, verbose=FALSE)
#############################################

# remove points that lay in the ocean - before: 1479293 rows
OcenDAta <- sp::over(coords, oceans)
length(OcenDAta$scalerank) #Remains at 1479293 rows
coords <- coords[which(is.na(OcenDAta$scalerank)),] #Keeps only the rows where scalerank == NA
str(coords)
clean_coords<-data.frame(coordinates(coords)) #Convert to a dataframe object
length(clean_coords$Long) #1377144

# Remove Hawaii - is this step necessary? NO as the dataset already does not contain any data from Hawaii. Because no native ant species exist in Hawaii
clean_coords_H_removed <- clean_coords[!((round(clean_coords$Lat) %in% seq(18, 28, by = 1)) & (round(clean_coords$Long) %in% seq(-178, -154, by = 1))),] 

plot(Kass_orgin$lat_opt~Kass_orgin$lon_opt, pch = ".", cex = 1) #Original Kass occurrence data
plot(clean_coords$Lat~clean_coords$Long, pch = ".", cex = 1) #Without islands


#Plotting
# Load world map data
world_map <- map_data("world")

##################
# Kass original occurrences
original <- ggplot() +
  # Plot world map
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightblue", color = "white") +
  # Highlight region of Hawaii
  geom_rect(aes(xmin = -160, xmax = -155, ymin = 18.5, ymax = 20.5), fill = "black", alpha = 0.05) +
  # Plot coordinates from your dataframe
  geom_point(data = Kass_orgin, aes(x = lon_opt, y = lat_opt), color = "red", size = 0.05) +
  # Set map boundaries
  coord_fixed(1.3) +
  # Add labels and title
  labs(title = "Original data", x = "Longitude", y = "Latitude") +
  # Set longitude values to appear every 5 units
  scale_x_continuous(breaks = seq(-180, 180, by = 10)) +
  scale_y_continuous(breaks = seq(-90, 90, by = 10)) +
  theme_minimal()


##################
# Without islands
minus_islands <- ggplot() +
  # Plot world map
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightblue", color = "white") +
  # Plot coordinates from your dataframe
  geom_point(data = clean_coords, aes(x = Long, y = Lat), color = "red", size = 0.05) +
  # Set map boundaries
  coord_fixed(1.3) +
  # Add labels and title
  labs(title = "Without islands", x = "Longitude", y = "Latitude") +
  # Set longitude values to appear every 5 units
  scale_x_continuous(breaks = seq(-180, 180, by = 10)) +
  theme_minimal()

# ggsave("original.pdf", plot = original, path = "/Users/louis.bell-roberts/Desktop/")
# ggsave("minus_island.pdf", plot = minus_islands, path = "/Users/louis.bell-roberts/Desktop/")


length(Kass_orgin$valid_species_name) #1479293
length(clean_coords$Long) #1377144

#Identify the occurrences in Kass_orgin that are islands
Kass_orgin <- cbind(Kass_orgin, paste0(Kass_orgin$lon_opt, ",", Kass_orgin$lat_opt))
Kass_orgin <- Kass_orgin %>%
  rename(new_id = `paste0(Kass_orgin$lon_opt, ",", Kass_orgin$lat_opt)`)

clean_coords <- cbind(clean_coords, paste0(clean_coords$Long, ",", clean_coords$Lat))
clean_coords <- clean_coords %>%
  rename(new_id = `paste0(clean_coords$Long, ",", clean_coords$Lat)`)

# Create the 'Island' column in Kass_orgin
Kass_orgin$Island <- ifelse(Kass_orgin$new_id %in% clean_coords$new_id, 0, 1)

#############################################
### Compare to the La Rich datasets ###
#############################################

#Read in la rich data - all occurrence records
La_Rich_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/La_Richeliere_dryad/AntOccurrenceGABI_final.csv")

#Read in CS data
ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Subsetted/Primary/Data/CS/New_transformations/Post_Koos/Post_Koos_AntEnv_data_FINAL.csv")
ant_data$species <- gsub("_", ".", ant_data$species)

##Identify species present in CS data that are absent from Kass 2022
setdiff(ant_data$species, Kass_orgin$valid_species_name) #9 species where CS data is present but no climatic data - Kass
setdiff(ant_data$species, La_Rich_data$valid_species_name) #49 species where CS data is present but no climatic data - la rich

###Identify which database has more species - Kass has more
length(unique(Kass_orgin$valid_species_name)) #14347
length(unique(La_Rich_data$valid_species_name)) #10060

###Which has more occurrences per species - Kass has more: 103 vs 74
length(La_Rich_data$valid_species_name)/length(unique(La_Rich_data$valid_species_name)) #73.87783 occurrences per species
length(Kass_orgin$valid_species_name)/length(unique(Kass_orgin$valid_species_name)) #103.1082 occurrences per species


#Write Kass_2022 with island occurrences to file
# write.csv(Kass_orgin, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Ant_environment/Kass_2022/Kass_data_with_island.csv")

