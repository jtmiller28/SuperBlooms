### Time Series Animation
### Project: SuperBlooms
### Author: JT Miller
### Date: 04-07-2025

# A script to build gifs showing how flowering is observed across superbloom years within our study region. 

# load libraries
library(ggplot2)
library(data.table)
library(tidyverse)
library(sf)
library(RColorBrewer)
library(gganimate)

# load the subsetted data
flower_data <- fread("/home/jt-miller/Gurlab/SuperBlooms/data/processed/desert-observations-hexed.csv")

# filter data for a known superbloom year
flower_data_sb <- flower_data %>% 
  filter(year == 2019)

# bin pheno data into increments of 7 (weeks)
flower_data_wsb <- flower_data_sb %>% # weekly bins
  mutate(week_bin = cut(doy, breaks = seq(1, 366, by = 7), right = FALSE, labels = FALSE)) %>% 
  filter(!is.na(week_bin)) # 365 doesnt quite make it into these categories, so Im electing to throw them away for exploratory purposes


# load in gridded map of study region 
hex_50km <- read_rds("/home/jt-miller/Gurlab/SuperBlooms/data/processed/50km-hexed-map.rds")

# bin data by hexcell, and week. Then summarize the number of obs & unique sp occuring within those bins. 
wk_hex50_bin_summary <- flower_data_wsb %>% 
  group_by(hex50_id, week_bin) %>% 
  summarize(
    num_observations = n(), 
    uniqueSpecies = n_distinct(taxon_name),
    .groups = "drop"
  )

# merge these data to create a spatial summary of obs & richness per cell per doy 
hex_50_wk_tbl <- merge(hex_50km, wk_hex50_bin_summary, by = "hex50_id", all.x = TRUE)

# sum up the total obs across all spatial cells to build a histogram of phenology across the deserts (for histogram visual)
hex_50_wk_hist_summary <- hex_50_wk_tbl %>% 
  sf::st_drop_geometry() %>% 
  dplyr::group_by(week_bin) %>% 
  dplyr::summarise(total_obs = sum(num_observations, na.rm = TRUE))

# create a map
map_plot <- ggplot() +
  geom_sf(hex_50km, mapping = aes()) +
  geom_sf(hex_50_wk_tbl, mapping = aes(fill = log(num_observations + 1)), color = NA) +
  theme_bw() + # remove boring grey background 
  scale_fill_gradientn(
    colours = colorRampPalette(RColorBrewer::brewer.pal(11, "YlOrRd"))(9), 
    na.value = "#808080", 
    name = "log(N obs + 1)"
  ) +
  labs(
    title = "Weekly Aggregated Observation Intensity per 50km grid cell",
    subtitle = "Day of Year: {frame_time}", 
    caption = "iNaturalist Flowering Data annotated by PhenoVision", 
    x = NULL, y = NULL
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"), 
    plot.subtitle = element_text(size = 14),
    legend.position = "right"
  ) +
  transition_time(week_bin)

# create an animation a save. 
map_anim <- animate(map_plot, nframes = 366, fps = 5, units = "in", width = 8, height = 8, res = 100, renderer = magick_renderer())
anim_save("/home/jt-miller/Gurlab/SuperBlooms/figures/hex_obs_by_wk2019.gif", animation = map_anim)
