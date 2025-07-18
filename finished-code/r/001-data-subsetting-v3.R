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
library(arrow)

# Load shp files of the Mojave & Sonoran Deserts, as well as sociopolitical boundaries
#na_shapes <- sf::read_sf("/blue/guralnick/millerjared/SuperBlooms/data/raw/NA-ecoregions/NA_CEC_Eco_Level3.shp")

# filter to only include these two regions
# sonoran_shp <- na_shapes %>%
#   filter(NA_L3NAME == "Sonoran Desert")
# mojave_shp <- na_shapes %>%
#   filter(NA_L3NAME == "Mojave Basin and Range")

# create an equal area projection
PROJ_CRS <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
sf::sf_use_s2(FALSE) # turning off for flat spatial mapping
# Load a basemap of states/provinces for the study region
basemap_wgs <- sf::st_read("./data/raw/10m_cultural/10m_cultural/ne_10m_admin_1_states_provinces.shp") %>%
  dplyr::filter(name %in% c("California"))
# Create a bounding box that covers all longitudes but only up to 35 degrees latitude
bbox <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 36), crs = sf::st_crs(basemap_wgs))
bbox_sf <- sf::st_as_sfc(bbox)
# Intersect your basemap_wgs with this bounding box
basemap_cut <- sf::st_intersection(basemap_wgs, bbox_sf)
# Now transform the clipped map to your target CRS
basemap <- sf::st_transform(basemap_cut, PROJ_CRS)
# desert_shps <- rbind(sonoran_shp, mojave_shp)
# desert_shps <- sf::st_transform(desert_shps, PROJ_CRS)
# desert_shps <- st_combine(desert_shps)
# Plot the study regions
ggplot() +
  geom_sf(data = basemap, fill = "white", color = "black") +  # Basemap layer
  #geom_sf(data = desert_shps, aes(), color = "black", fill = "darkorange") +  # Study regions with fill
  ggtitle("The Sonoran & Mojave Deserts") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"  # Adjust legend position if needed
  )

# Save as SVG
#ggsave("/blue/guralnick/millerjared/SuperBlooms/figures/mojave-sonoran-study-region.svg",
#width = 15, height = 10)

# Save as a png
#ggsave("/blue/guralnick/millerjared/SuperBlooms/figures/mojave-sonoran-study-region.png",
#width = 15, height = 10)
# create a hex-cell loop for different grid cell sizes
# hex_cell_sizes <- c(5,10,25,50,100)
# hexed_maps_storage <- list()
# for(i in 1:length(hex_cell_sizes)){
#   hex_map <- sf::st_make_grid(basemap,
#                               cellsize = hex_cell_sizes[i]*1000, #  5km hex cells
#                               square = FALSE,
#                               flat_topped = FALSE) %>%
#     sf::st_as_sf() %>%
#     dplyr::mutate(!!paste0("hex", hex_cell_sizes[i], "_id") := row_number()) %>%
#     #sf::st_intersection(desert_shps) %>%
#     dplyr::mutate(area = sf::st_area(.))
#   
#   # additionally, save hex-map as a rds object
#   saveRDS(hex_map, paste0("/blue/guralnick/millerjared/SuperBlooms/data/processed/", hex_cell_sizes[i], "km", "-hexed-map.rds"))
#   hexed_maps_storage[[i]] <- hex_map
#   
#   # Make a gridded map of these regions
#   ggplot() +
#     geom_sf(data = basemap, fill = "white", color = "black") +  # Basemap layer
#     geom_sf(data = desert_shps, aes(), fill = "darkorange", color = "black") +  # Study regions with fill
#     geom_sf(data = hex_map, mapping = aes(), alpha = 0.3) +
#     ggtitle(paste0("The Sonoran & Mojave Deserts ", hex_cell_sizes[i], "km Grid")) +
#     theme_bw() +
#     theme(
#       plot.title = element_text(hjust = 0.5),
#       legend.position = "right"  # Adjust legend position if needed
#     )
#   ggsave(paste0("/blue/guralnick/millerjared/SuperBlooms/figures/", "desert-", hex_cell_sizes[i], "km-grid.svg"),width = 15, height = 10)
#   ggsave(paste0("/blue/guralnick/millerjared/SuperBlooms/figures/", "desert-", hex_cell_sizes[i], "km-grid.png"),width = 15, height = 10)
# }

# load necessary hexmap 
hex_map <- readRDS("./data/processed/50km-Scal-hexed-map.rds")
hex_map_wgs84 <- sf::st_transform(hex_map, crs = 4326)

area_limits_longlat <- st_bbox(hex_map_wgs84)

# Load the phenovision data
# process the phenovision data into usable parquet
pheno_table <- open_dataset(
  sources = "/home/jt-miller/Gurlab/phenovision-data/phenobase-annotations-01-2025.csv",
  format = "csv"
)

pheno_table <- pheno_table |>
  select(
    uuid = "0cf717af-4366-48df-8909-ea5fe454e1a0", 
    datasource = iNaturalist, 
    date = "1806-12-31", 
    doy = "365", 
    year = "1806",
    lat = "-34.007", 
    long = "18.45", 
    coord_uncertainty = "110", 
    family = "Proteaceae", 
    taxon_id = "4339",
    genus = "Leucadendron", 
    taxon_name = "Leucadendron grandiflorum", 
    taxon_rank = "species", 
    obs_type = "MachineObservation", 
    phenology_trait = "flower", 
    img_url = "https://inaturalist-open-data.s3.amazonaws.com/photos/15571385/small.jpg", 
    photo_url = "https://www.inaturalist.org/photos/15571385",
    obs_url = "https://www.inaturalist.org/observations/c3400d48-98f9-41ca-9abf-e3b61a49dd62", 
    confidence = "High", 
    family_probability = "10.57967/hf/2763", 
    probability = "0.9658203125", 
    detected = "Detected", 
    V23 = "0.07167550126757317", 
    V24 = "0.9771598808341608",  
    V25 = "0.945379119612814"
  ) 

pheno_table |> write_dataset(
  path = "./data/processed/phenovision-annotations-01-2025.parquet",
  format = "parquet")

pheno_table <- open_dataset("./data/processed/phenovision-annotations-01-2025.parquet",
  format = "parquet"
)

# Restrict data to study area, and 2017-2023
pheno_table <- pheno_table |>
  filter(
    lat >= area_limits_longlat[2] & lat <= area_limits_longlat[4],
    long >= area_limits_longlat[1] & long <= area_limits_longlat[3],
    year >= 2017 & year <= 2023
  )
# load in precut table 
pheno_dt <- pheno_table |> collect()

# bring in wcvp trait data to limit our data to only annuals 
wcvp_names <- fread("./data/processed/wcvp_name_trait_info/wcvp_names.csv")
wcvp_names <- wcvp_names %>% 
  select(taxon_name, lifeform_description)
wcvp_names %>% 
  group_by(lifeform_description) %>%
  summarize(n = n())

pheno_dt_lifeforms <- pheno_dt %>%
  left_join(wcvp_names, by = c("taxon_name" = "taxon_name"))

pheno_dt_annuals <- pheno_dt_lifeforms %>%
  filter(str_detect(lifeform_description, regex("annual", ignore_case = TRUE)))

#pheno_dt <- fread("/home/jt-miller/Gurlab/phenovision-data/phenobase-annotations-01-2025.csv")
# colnames(pheno_dt) = c("uuid", "datasource", "date", "doy", "year",
#                        "lat", "long", "coord_uncertainty", "family", "taxon_id",
#                        "genus", "taxon_name", "taxon_rank", "obs_type", "phenology_trait",
#                        "img_url", "photo_url", "obs_url", "confidence", "family_probability",
#                        "probability", "detected", "V23", "V24", "V25")

# make the data spatial
pheno_sf <- sf::st_as_sf(pheno_dt_annuals,
                         coords = c("long", "lat"),
                         crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
                         remove = FALSE) %>%
  st_transform(PROJ_CRS) # tranform to equal area projection to match the shp files

# use shapes to filter out occurrences to only occur within the deserts range
# socal_obs <- pheno_sf[basemap,]
socal_obs <- pheno_sf[hex_map,]
socal_obs <- socal_obs %>% mutate(my_id = 1:n())

# clean data for only flowering records as these are what are relevant for the question at hand.
socal_obs <- socal_obs %>% filter(phenology_trait == "flower")


# use a for-loop to assign obs to grid cells for each iteration of grid type w/label

#for(i in 1:length(hexed_maps_storage)){
  # join with hexed map to produce hexed observation data
  socal_obs_hexed <- sf::st_join(socal_obs, hex_map, left = FALSE) # setting left to false means we discard non-hexed points (should be none)
  # select relevant info so we can bind with the dataset
  socal_obs_id_field <- socal_obs_hexed %>% select(paste0("hex", "50", "_id"), my_id) %>% st_drop_geometry()
  # merge, use the id we set up prior as the key
  socal_obs <- merge(socal_obs, socal_obs_id_field, by = "my_id")
  socal_obs <- socal_obs %>% distinct()
#}

# drop spatial & added fields
socal_obs <- socal_obs %>% select(-my_id) %>% st_drop_geometry()

# remove redundant annotated flowering data for unique observational records
socal_obs <- dplyr::distinct(socal_obs, obs_url, .keep_all = TRUE)

# standardize temporal fields for background obs data
socal_obs <- socal_obs %>%
  mutate(eventDateParsed = parse_date_time(date, orders = c("ymd_HMS", "ymd_HM", "ymd"))) %>%
  mutate(doy = lubridate::yday(eventDateParsed)) %>%
  mutate(year = lubridate::year(eventDateParsed)) %>%
  mutate(month = lubridate::month(eventDateParsed))

# remove city nature challenge
# remove city nature challenge
# Make a vector of all CNC dates
library(lubridate)
cnc_date_ranges <- list(
  seq(ymd("2016-04-14"), ymd("2016-04-21"), by = "days"),
  seq(ymd("2017-04-14"), ymd("2017-04-18"), by = "days"),
  seq(ymd("2018-04-27"), ymd("2018-04-30"), by = "days"),
  seq(ymd("2019-04-26"), ymd("2019-04-29"), by = "days"),
  seq(ymd("2020-04-24"), ymd("2020-04-27"), by = "days"),
  seq(ymd("2021-04-30"), ymd("2021-05-03"), by = "days"),
  seq(ymd("2022-04-29"), ymd("2022-05-02"), by = "days"),
  seq(ymd("2023-04-28"), ymd("2023-05-01"), by = "days"),
  seq(ymd("2024-04-26"), ymd("2024-04-29"), by = "days")
)
cnc_dates <- unlist(cnc_date_ranges)
cnc_dates <- as.Date(cnc_dates, origin = "1970-01-01")

# Remove CNC rows efficiently
socal_obs <- socal_obs %>%
  filter(!(as.Date(eventDateParsed) %in% cnc_dates))


# write out the dataset for further wrangling
fwrite(socal_obs, "/home/jt-miller/Gurlab/SuperBlooms/data/processed/soCal-annual-data.csv")
