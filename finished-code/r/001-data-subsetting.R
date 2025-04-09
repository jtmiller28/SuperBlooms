### 001-data-subsetting
### Project: Superblooms 
### Date: 03/28/2025
### Author: JT Miller

# Purpose: Take the phenovision data and subset it to the geographic bounds of the Sonoran & Mojave Deserts 

# load libraries
library(sf)
library(tidyverse)
library(data.table)
library(ggplot2)

# Load shp files of the Mojave & Sonoran Deserts, as well as sociopolitical boundaries
na_shapes <- sf::read_sf("/blue/guralnick/millerjared/SuperBlooms/data/raw/NA-ecoregions/NA_CEC_Eco_Level3.shp")

# filter to only include these two regions
sonoran_shp <- na_shapes %>% 
  filter(NA_L3NAME == "Sonoran Desert")
mojave_shp <- na_shapes %>% 
  filter(NA_L3NAME == "Mojave Basin and Range")

# create an equal area projection
PROJ_CRS <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
sf::sf_use_s2(FALSE) # turning off for flat spatial mapping
# Load a basemap of states/provinces for the study region
basemap_wgs <- sf::st_read("/blue/guralnick/millerjared/SuperBlooms/data/raw/10m_cultural/10m_cultural/ne_10m_admin_1_states_provinces.shp") %>%
  dplyr::filter(name %in% c("California", "Nevada", "Arizona",  "Utah",
                            "Sonora",
                            "Baja California", "Baja California Sur",
                            "Sinaloa"))
# transform the projections to follow equal area
basemap <- sf::st_transform(basemap_wgs, PROJ_CRS)
sonoran_shp <- sf::st_transform(sonoran_shp, PROJ_CRS)
mojave_shp <- sf::st_transform(mojave_shp, PROJ_CRS)
desert_shps <- st_union(sonoran_shp, mojave_shp) # required step, as we need to assure that the boundaries of the shapefiles do not clip hexcells later in the code. 

# Plot the study regions
ggplot() + 
  geom_sf(data = basemap, fill = "white", color = "black") +  # Basemap layer
  geom_sf(data = desert_shps, aes(), color = "black", fill = "darkorange") +  # Study regions with fill
  ggtitle("The Sonoran & Mojave Deserts") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"  # Adjust legend position if needed
  )

# Save as SVG
ggsave("/blue/guralnick/millerjared/SuperBlooms/figures/mojave-sonoran-study-region.svg", 
       width = 15, height = 10)
# create a hex-cell loop for different grid cell sizes 
hex_cell_sizes <- c(5,10,25,50,100)
hexed_maps_storage <- list()
for(i in 1:length(hex_cell_sizes)){
hex_map <- sf::st_make_grid(basemap,
                              cellsize = hex_cell_sizes[i]*1000, #  5km hex cells
                              square = FALSE,
                              flat_topped = FALSE) %>%
  sf::st_as_sf() %>%
  dplyr::mutate(!!paste0("hex", hex_cell_sizes[i], "_id") := row_number()) %>% 
  sf::st_intersection(desert_shps) %>%
  dplyr::mutate(area = sf::st_area(.))

# additionally, save hex-map as a rds object
saveRDS(hex_map, paste0("/blue/guralnick/millerjared/SuperBlooms/data/processed/", hex_cell_sizes[i], "km", "-hexed-map.rds"))
hexed_maps_storage[[i]] <- hex_map

# Make a gridded map of these regions 
ggplot() + 
  geom_sf(data = basemap, fill = "white", color = "black") +  # Basemap layer
  geom_sf(data = desert_shps, aes(), fill = "darkorange", color = "black") +  # Study regions with fill
  geom_sf(data = hex_map, mapping = aes(), alpha = 0.3) +
  ggtitle(paste0("The Sonoran & Mojave Deserts ", hex_cell_sizes[i], "km Grid")) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"  # Adjust legend position if needed
  )
ggsave(paste0("/blue/guralnick/millerjared/SuperBlooms/figures/", "desert-", hex_cell_sizes[i], "km-grid.svg"),width = 15, height = 10)
}



# Load the phenovision data
pheno_dt <- fread("/blue/guralnick/millerjared/SuperBlooms/data/raw/phenobase-annotations-01-2025.csv")
colnames(pheno_dt) = c("uuid", "datasource", "date", "doy", "year",
              "lat", "long", "coord_uncertainty", "family", "taxon_id",
              "genus", "taxon_name", "taxon_rank", "obs_type", "phenology_trait",
              "img_url", "photo_url", "obs_url", "confidence", "family_probability",
              "probability", "detected", "V23", "V24", "V25")

# make the data spatial 
pheno_sf <- sf::st_as_sf(pheno_dt, 
                             coords = c("long", "lat"),
                             crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
                             remove = FALSE) %>% 
  st_transform(PROJ_CRS) # tranform to equal area projection to match the shp files

# use shapes to filter out occurrences to only occur within the deserts range
desert_obs <- pheno_sf[desert_shps,]
desert_obs <- desert_obs %>% mutate(my_id = 1:n())

# clean data for only flowering records as these are what are relevant for the question at hand. 
desert_obs <- desert_obs %>% filter(phenology_trait == "flower")


# use a for-loop to assign obs to grid cells for each iteration of grid type w/label 

for(i in 1:length(hexed_maps_storage)){
  # join with hexed map to produce hexed observation data
  desert_obs_hexed <- sf::st_join(desert_obs, hexed_maps_storage[[i]], left = FALSE) # setting left to false means we discard non-hexed points (should be none)
  # select relevant info so we can bind with the dataset
  desert_obs_id_field <- desert_obs_hexed %>% select(paste0("hex", hex_cell_sizes[i], "_id"), my_id) %>% st_drop_geometry()
  # merge, use the id we set up prior as the key
  desert_obs <- merge(desert_obs, desert_obs_id_field, by = "my_id") 
  desert_obs <- desert_obs %>% distinct()
}

# drop spatial & added fields
desert_obs <- desert_obs %>% select(-my_id) %>% st_drop_geometry()

# write out the dataset for further wrangling 
fwrite(desert_obs, "/blue/guralnick/millerjared/SuperBlooms/data/processed/desert-observations-hexed.csv")
