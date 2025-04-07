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
    title = "2019 Weekly Aggregated Observation Intensity \n per 50km grid cell",
    subtitle = "Week Bin: {frame_time}", 
    caption = "iNaturalist Flowering Data annotated by PhenoVision", 
    x = NULL, y = NULL
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"), 
    plot.subtitle = element_text(size = 14),
    legend.position = "right"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  transition_time(week_bin)

# create an animation a save. 
map_anim <- animate(map_plot, nframes = 52, fps = 1, units = "in", width = 8, height = 8, res = 100, renderer = magick_renderer())
anim_save("/home/jt-miller/Gurlab/SuperBlooms/figures/hex_obs_by_wk2019.gif", animation = map_anim)

# create an alternative field to transition by (as to keep the overall histogram static in image)
week_df <- data.frame(week_step_bin = hex_50_wk_hist_summary$week_bin)
# create a hist plot 
hist_plot <- ggplot(hex_50_wk_hist_summary, mapping = aes(x = week_bin, y = total_obs)) + 
  geom_col(fill = "steelblue", color = "black", width = 1) +
  geom_vline(week_df, mapping = aes(xintercept = week_step_bin), color = "red", linewidth = 1) +
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
  view_static() +
  transition_time(week_step_bin)

hist_anim <- animate(hist_plot, nframes = 52, fps = 1, units = "in", width = 8, height = 8, res = 100, renderer = magick_renderer())
anim_save("/home/jt-miller/Gurlab/SuperBlooms/figures/hist_obs_by_wk2019.gif", animation = hist_anim)

library(magick)
week_df <- filter(week_df, !is.na(week_step_bin))
# combine gifs into an object
comb_gif <- magick::image_append(c(map_anim[1], hist_anim[1]))

for(i in 2:max(week_df$week_step_bin)){
  combined <- image_append((c(map_anim[i], hist_anim[i])))
  comb_gif <- c(comb_gif, combined)
}

anim_save("/home/jt-miller/Gurlab/SuperBlooms/figures/hist-hex-obs-by-wk2019.gif", comb_gif )
