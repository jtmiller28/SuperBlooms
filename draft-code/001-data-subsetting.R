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

# plot the regions as a visual of our study regions
ggplot() + 
  geom_sf(basemap, mapping = aes()) +
  geom_sf(sonoran_shp, mapping = aes(), fill = "darkorange") + 
  geom_sf(mojave_shp, mapping = aes(), fill = "goldenrod") + 
  ggtitle("The Sonoran & Mojave Deserts") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
  

# create a svg 
ggsave("/blue/guralnick/millerjared/SuperBlooms/figures/mojave-sonoran-study-region.svg", width = 15, height = 10)

# Load the phenovision data
pheno_dt <- fread("/blue/guralnick/millerjared/SuperBlooms/data/raw/phenobase-annotations-01-2025.csv", nrows = 100)
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
mojave_obs <- pheno_sf[mojave_shp, ]
sonoran_obs <- pheno_sf[sonoran_shp,]

# label these incase we care down the line what desert they are from, then drop spatial properties
mojave_obs <- mojave_obs %>% mutate(ecoregion = "mojave") %>% st_drop_geometry()
sonoran_obs <- sonoran_obs %>% mutate(ecoregion = "sonoran") %>% st_drop_geometry()

# bring datasets together
desert_obs <- rbind(mojave_obs, sonoran_obs)

# write out the dataset for further wrangling 
fwrite(desert_obs, "/blue/guralnick/millerjared/SuperBlooms/data/processed/desert-observations-annot.csv")
