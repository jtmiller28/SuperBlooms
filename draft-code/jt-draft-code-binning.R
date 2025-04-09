### 002 Observation Binning 
### Project: SuperBlooms
### Author: JT Miller
### Date: 03-31-2025


### Purpose: Using prior superblooms in 2016, 2019, and 2023-2024, visualize a how superblooms look according to the observation record

# load packages
library(ggplot2)
library(data.table)
library(tidyverse)
library(sf)
library(RColorBrewer)

# load the subsetted data
flower_data <- fread("./data/processed/desert-observations-hexed.csv")
flower_data <- flower_data %>% distinct(obs_url, .keep_all=TRUE)

# grab observations in 2016, 2019, and 2023, construct bins: weekly, monthly
flower_data_sbs <- flower_data %>% 
  filter(year %in% c(2016, 2019, 2023))

flower_data_wsbs <- flower_data_sbs %>% # weekly bins
  mutate(week_bin = cut(doy, breaks = seq(1, 366, by = 7), right = FALSE, labels = FALSE)) %>% 
  filter(!is.na(week_bin)) # 365 doesnt quite make it into these categories, so Im electing to throw them away for exploratory purposes

flower_data_msbs <- flower_data %>% 
  mutate(month_bin = cut(doy, breaks = seq(1, 366, by = 30), right = FALSE, labels = FALSE)) %>% 
  filter(!is.na(month_bin)) 

# summarize the bin and num of observations, as well as uniq num of species per bin
non_bin_summary <- flower_data_wsbs %>% 
  group_by(year, doy) %>% 
  summarize(
    num_observations = n(), 
    uniqueSpecies = n_distinct(taxon_name),
    .groups = "drop"
  )

week_bin_summary <- flower_data_wsbs %>% 
  group_by(year, week_bin) %>% 
  summarize(
    num_observations = n(), 
    uniqueSpecies = n_distinct(taxon_name),
    .groups = "drop"
  )

# create a histogram style plot...
ggplot(non_bin_summary, aes(x = doy, y = num_observations, fill = factor(year))) +
  geom_col(position = "dodge") +
  scale_x_continuous(breaks  = seq(1, 366, by = 5)) + 
  labs(
    x = "Day of Year",
    y = "Number of Observations", 
    title = "Observations Over daily bins") + 
  facet_wrap(~year) +
theme_minimal()

ggplot(week_bin_summary, aes(x = week_bin, y = num_observations, fill = factor(year))) +
  geom_col(position = "dodge") +
  scale_x_continuous(breaks  = seq(1, 52, by = 4)) + 
  labs(
    x = "Week Bin (7-day intervals)",
    y = "Number of Observations", 
    title = "Observations Over Week bins") + 
  facet_wrap(~year)
  theme_minimal()
  
## Create a static figure showing the number of obs & sp richness per grid cell throughout a super bloom year
# load in the data for the hexcell map
hex_50km <- read_rds("./data/processed/50km-hexed-map.rds")
# for static figure (all of a year)
hex50_bin_summary <- flower_data_wsbs %>% 
    group_by(year, hex50_id) %>% 
    summarize(
      num_observations = n(), 
      uniqueSpecies = n_distinct(taxon_name),
      .groups = "drop"
    )

# filter data to a known superbloom year
hex50_bin_summary_2019 <- filter(hex50_bin_summary, year == 2019)

# merge the data 
hex50_2019_tbl <- merge(hex_50km, hex50_bin_summary_2019, by = "hex50_id", all.x = TRUE)
# reclass the tbl back into spatial format 
hex50_2019_tbl <- st_as_sf(hex50_2019_tbl, crs = st_crs(hex_50km))
# create a static plot showing sampling intensity over our study region
ggplot() + 
  geom_sf(hex50_2019_tbl, mapping = aes(fill = log(num_observations+1))) + 
  theme_bw() + 
  scale_fill_gradientn(
    colours = colorRampPalette((RColorBrewer::brewer.pal(11, "YlOrRd")))(9), 
    na.value = "#808080" # The line that denotes NAs as grey
  )  + 
  ggtitle("Observation Intensity per 50km grid cell in year 2019 Across the Mojave & Sonoran Deserts")

ggsave("./figures/obs-intensity-50km2019.jpeg", height = 10, width = 11)
# create a static plot showing number of unique species over our study region
ggplot() + 
  geom_sf(hex50_2019_tbl, mapping = aes(fill = log(uniqueSpecies+1))) + 
  theme_bw() + 
  scale_fill_gradientn(
    colours = colorRampPalette((RColorBrewer::brewer.pal(11, "YlOrRd")))(9), 
    na.value = "#808080" # The line that denotes NAs as grey
  )  + 
  ggtitle("Number of Observed Flowering Species per 50km grid cell in year 2019 Across the Mojave & Sonoran Deserts")
ggsave("./figures/sp-richness-50km2019.jpeg", height = 10, width = 11)

# try filter out to joshua tree & Anza-Borrego
shp_nps<- sf::read_sf("/home/jt-miller/Soltis-lab/useful-code-depository/nps-shps/Administrative Boundaries of National Park System Units.shp")
shp_joshuaTree <- shp_nps[grep("Joshua", shp_nps$UNIT_NAME),]
shp_sd <- sf::read_sf("/home/jt-miller/Soltis-lab/useful-code-depository/parks_datasd/parks_datasd.shp")
shp_anza<- shp_sd[grep("Anza", shp_sd$name),]
hex_50km <- readRDS("/home/jt-miller/Gurlab/SuperBlooms/data/processed/50km-hexed-map.rds")

shp_joshuaTree <- st_transform(shp_joshuaTree, crs = st_crs(hex_50km))
shp_anza <- st_transform(shp_anza, crs = st_crs(hex_50km))
ggplot() + 
  geom_sf(hex_50km, mapping = aes()) +
  geom_sf(shp_anza, mapping = aes(), fill = "orange") + 
  geom_sf(shp_joshuaTree, mapping = aes(), fill = "red")

# intersect grid_ids 
hex_50km_joshuaTree <- hex_50km[shp_joshuaTree, ]
hex_50km_anza <- hex_50km[shp_anza, ]

# create a desert parks hexcell 
desert_parks_50km <- rbind(hex_50km_joshuaTree, hex_50km_anza)

# remove sonoran & mojave designations 
desert_parks_50km <- desert_parks_50km %>% 
  select(hex50_id, area, x) %>% 
  distinct() # remove the issue where hexids are duplicated for hexcells that overlap sonoran & mojave
# load the subsetted data
flower_data <- fread("./data/processed/desert-observations-hexed.csv")
flower_data <- flower_data %>% distinct(obs_url, .keep_all=TRUE) %>% distinct(uuid, img_url, photo_url, obs_url, hex50_id, .keep_all = TRUE) # deal with some duplication issues when creating grid cells 

flower_data_sb <- flower_data %>% 
  filter(year == 2019)

flower_data_wsb <- flower_data_sb %>% # weekly bins
  mutate(week_bin = cut(doy, breaks = seq(1, 366, by = 7), right = FALSE, labels = FALSE)) %>% 
  filter(!is.na(week_bin)) # 365 doesnt quite make it into these categories, so Im electing to throw them away for exploratory purposes

wk_hex50_bin_summary <- flower_data_wsb %>% 
  group_by(hex50_id, week_bin) %>% 
  summarize(
    num_observations = n(), 
    uniqueSpecies = n_distinct(taxon_name),
    .groups = "drop"
  )

# merge these data to create a spatial summary of obs & richness per cell per doy 
desert_parks_wk_tbl <- merge(desert_parks_50km, wk_hex50_bin_summary, by = "hex50_id", all.x = TRUE)
find_blooms <- merge(flower_data_wsb, desert_parks_50km, by = "hex50_id")
desert_parks_wk_hist_summary <- desert_parks_wk_tbl  %>% 
  sf::st_drop_geometry() %>% 
  dplyr::group_by(week_bin, hex50_id) %>% 
  dplyr::summarise(total_obs = sum(num_observations, na.rm = TRUE))

# plot the histogram
ggplot(desert_parks_wk_hist_summary, mapping = aes(x = week_bin, y = total_obs)) + 
  geom_col(fill = "steelblue", color = "black", width = 1) +
  scale_x_continuous(breaks = seq(1, 52, by = 1)) +
  labs(
    x = "Week bin", 
    y = "Number of Observations", 
    title = "2019 Flowering Observations Across the Year"
  ) + 
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~hex50_id)

# use prior knowledge of some taxa listed in superblooms:
superbloom_taxa <- c("Eschsocholzia californica", "Phacelia campanularia", "Encelia farinosa", "Lupinus",
                     "Abronia", "Geraea canescens", "Chylismia brevipes", "Plagiobothrys", "Hesperocallis")

# extract these records from the data 
test_blooms <- flower_data_wsb %>% 
  filter(str_detect(taxon_name, str_c("^", superbloom_taxa, collapse = "|")))

bloom_wk_tbl <- merge(hex_50km, wk_hex50_bin_summary, by = "hex50_id", all.x = TRUE)

test_blooms_wk_hist_summary <- bloom_wk_tbl  %>% 
  sf::st_drop_geometry() %>% 
  dplyr::group_by(week_bin, hex50_id) %>% 
  dplyr::summarise(total_obs = sum(num_observations, na.rm = TRUE))

# find top ten events
top_10_hexcells <- test_blooms_wk_hist_summary %>% 
  group_by(hex50_id) %>% 
  summarize(total_bin_obs = sum(total_obs)) %>% 
  arrange(desc(total_bin_obs)) %>% 
  head(10) %>% 
  select(hex50_id)

top_bloom_wk_hist_summary <- filter(test_blooms_wk_hist_summary, hex50_id %in% top_10_hexcells$hex50_id)

ggplot(top_bloom_wk_hist_summary, mapping = aes(x = week_bin, y = total_obs)) + 
  geom_col(fill = "goldenrod", color = "black", width = 1) +
  scale_x_continuous(breaks = seq(1, 52, by = 1)) +
  labs(
    x = "Week bin", 
    y = "Number of Observations", 
    title = "2019 Flowering Observations Across the Year"
  ) + 
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~hex50_id)

# find hexcells:
ggplot() + 
  geom_sf(hex_50km, mapping = aes()) + 
  geom_sf(filter(hex_50km, hex50_id %in% top_10_hexcells$hex50_id), mapping = aes(), fill = "goldenrod", alpha = 0.5) +
  geom_sf_text(filter(hex_50km, hex50_id %in% top_10_hexcells$hex50_id), mapping = aes(label = hex50_id), size = 3, color = "black")

# check out data in high obs density bins
test <- filter(test_blooms, hex50_id == 740)
test <- filter(test, week_bin == 12)
high_density_week <- desert_parks_wk_hist_summary[, max(desert_parks_wk_hist_summary$total_obs)]
flower_data_test <- flower_data_wsb %>% 
  filter(week_bin == max(desert_parks_wk_hist_summary))