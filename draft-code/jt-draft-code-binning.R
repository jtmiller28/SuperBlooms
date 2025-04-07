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
library(gganimate)
library(patchwork)
library(gifski)
library(cowplot)

# load the subsetted data
flower_data <- fread("./data/processed/desert-observations-hexed.csv")

# filter data to only include flowering records (remove once subset finishes.)
flower_data <- flower_data %>% filter(phenology_trait == "flower")

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

## Create a animated graph showing how hexcells and obs change over a years span 
# load in the data for the hexcell map
hex_50km <- read_rds("./data/processed/50km-hexed-map.rds")
# for animated figure: bin the data by: year, hexcell, and doy. 
doy_hex50_bin_summary <- flower_data_wsbs %>% 
  group_by(year, hex50_id, doy) %>% 
  summarize(
    num_observations = n(), 
    uniqueSpecies = n_distinct(taxon_name),
    .groups = "drop"
  )

# filter to a known superbloom year
doy_hex50_bin_2019_summary <- filter(doy_hex50_bin_summary, year == 2019)

# join these summaries with the hexcell table 
hex_50_2019_doy_tbl <- merge(hex_50km, doy_hex50_bin_2019_summary, by = "hex50_id", all.x = TRUE)

# aggregate total obs per doy for a histogram visual 
hex50_2019_doy_summary <- hex_50_2019_doy_tbl %>% 
  sf::st_drop_geometry() %>% 
  group_by(doy) %>% 
  summarize(total_obs = sum(num_observations, na.rm = TRUE))

# using gganimate, illustrate how flower observation occurs over the span of a year.
hex_50_2019_doy_tbl <- hex_50_2019_doy_tbl %>% 
  mutate(doy_tracker = doy)
map_plot <- ggplot() + # feed ggplot the data
  geom_sf(hex_50km, mapping = aes()) +
  geom_sf(hex_50_2019_doy_tbl, mapping = aes(fill = log(num_observations + 1)), color = NA) + # assign aes fill to be by num of obs, allowing for 0s to not na, nas will be grey
  theme_bw() + # remove boring grey background
  scale_fill_gradientn(
    colours = colorRampPalette(RColorBrewer::brewer.pal(11, "YlOrRd"))(9), 
    na.value = "#808080", 
    name = "log(N obs + 1)"
  ) +
  labs(
    title = "Observation Intensity per 50km grid cell",
    subtitle = "Day of Year: {frame_time}", 
    caption = "iNaturalist Flowering Data annotated by PhenoVision", 
    x = NULL, y = NULL
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"), 
    plot.subtitle = element_text(size = 14),
    legend.position = "right"
  ) +
  transition_time(doy_tracker) #+ # assign doy as our animation unit
  #ease_aes("linear")


# create hist plot
doy_line_df <- data.frame(doy_tracker = 1:366)
hist_plot <- ggplot() + 
  geom_vline(doy_line_df, mapping = aes(xintercept = doy_tracker), color = "red", linewidth = 1) + # build doy tracking line
  geom_col(hex50_2019_doy_summary, mapping = aes(x = doy, y = total_obs), fill = "steelblue", width = 1) +
  scale_x_continuous(breaks = seq(0, 366, by = 30)) +
  labs(
    x = "Day of Year", 
    y = "Number of Observations", 
    title = "Flowering Observations Across the Year"
  ) + 
  theme_minimal() +
  transition_time(doy_tracker)
  
# combine plots 
combined_plot <- map_plot + hist_plot #+ 
  #plot_layout(heights = c(2,1)) + 
  #transition_time(doy_tracker) #+ 
 # ease_aes("linear")

test_anim <- animate(combined_plot, nframes = 366, fps =5, units = "in", width = 4, height = 4, res = 200, renderer = gifski_renderer())
anim_save("./figures/hex_obs_by_doy2019.gif", animation = test_anim)
# animate and save 
map_anim <- animate(map_plot, nframes = 366, fps = 5, units = "in", width = 4, height = 4, res = 200, renderer = magick_renderer())
anim_save("./figures/hex_obs_by_doy2019.gif", animation = map_anim)

hist_anim <- animate(hist_plot, nframes = 366, fps = 5, width = 1000, height = 800, renderer = gifski_renderer())
anim_save("./figures/hex_obs_by_doy2019-hist.gif")
  
# attempt 2! 
map_plot <- ggplot() +
  geom_sf(data = hex_50km, fill = NA, color = "grey60", size = 0.1) +
  geom_sf(data = hex_50_2019_doy_tbl, 
          mapping = aes(fill = log(num_observations + 1)), color = NA) +
  scale_fill_gradientn(
    colours = colorRampPalette(RColorBrewer::brewer.pal(11, "YlOrRd"))(9),
    na.value = "#808080",
    name = "log(N obs + 1)"
  ) +
  theme_bw() +
  labs(
    title = "Observation Intensity per 50km grid cell",
    subtitle = "Day of Year: {frame_time}",
    caption = "iNaturalist Flowering Data annotated by PhenoVision",
    x = NULL, y = NULL
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14),
    legend.position = "right"
  ) + 
  transition_time(doy) + # assign doy as our animation unit
  ease_aes("linear")

doy_line_df <- data.frame(doy = 1:366)
hist_plot <- ggplot(hex50_2019_doy_summary, aes(x = doy, y = total_obs)) + 
  geom_col(fill = "steelblue", width = 1) +
  geom_vline(data = doy_line_df, aes(xintercept = doy), color = "red", linewidth = 1) +
  scale_x_continuous(breaks = seq(0, 366, by = 30)) +
  labs(
    x = "Day of Year", 
    y = "Number of Observations", 
    title = "Flowering Observations Across the Year"
  ) + 
  theme_minimal() + 
  transition_time(doy) + # assign doy as our animation unit
  ease_aes("linear")

combined_plot <- map_plot / hist_plot +
  plot_layout(heights = c(2, 1)) +
  transition_time(doy) +
  ease_aes("linear")

anim <- animate(combined_plot, 
                nframes = 366, 
                fps = 5, 
                width = 1000, 
                height = 900, 
                renderer = gifski_renderer())

anim_save("./figures/hex_obs_by_doy2019.gif", anim)


# attemptign a purrr solution
library(ggplot2)
library(sf)
library(dplyr)
library(patchwork)
library(purrr)
library(magick)  # for stitching into gif
# OR gifski::save_gif() if you prefer

# Optional: make sure data has complete days per hex (avoid empty maps)
hex_50_2019_doy_tbl <- hex_50_2019_doy_tbl %>%
  tidyr::complete(doy = 1:366, nesting(hex50_id), fill = list(num_observations = 0)) %>% 
  st_as_sf(crs = st_crs(hex_50km))

# Directory to save frames
dir.create("./figures/frames", showWarnings = FALSE)

# Generate plots for each DOY
map(hist(1:366), function(current_doy) {
  
  # Filter map data for the current day
  map_data_day <- hex_50_2019_doy_tbl %>% 
    filter(doy == current_doy) %>% 
    st_as_sf
  
  # Map plot
  map_plot <- ggplot() +
    geom_sf(data = hex_50km, fill = NA, color = "grey60", size = 0.1) +
    geom_sf(data = map_data_day, mapping = aes(fill = log(num_observations + 1)), color = NA) +
    scale_fill_gradientn(
      colours = colorRampPalette(RColorBrewer::brewer.pal(11, "YlOrRd"))(9),
      na.value = "#808080",
      name = "log(N obs + 1)"
    ) +
    theme_bw() +
    labs(
      title = "Observation Intensity per 50km grid cell",
      subtitle = paste("Day of Year:", current_doy),
      caption = "iNaturalist Flowering Data annotated by PhenoVision",
      x = NULL, y = NULL
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14),
      legend.position = "right"
    )
  
  # Histogram plot (static bars + red line for current DOY)
  hist_plot <- ggplot(hex50_2019_doy_summary, mapping = aes(x = doy, y = total_obs)) +
    geom_col(fill = "steelblue", width = 1) +
    geom_vline(xintercept = current_doy, color = "red", linewidth = 1) +
    scale_x_continuous(breaks = seq(0, 366, by = 30)) +
    labs(
      x = "Day of Year",
      y = "Number of Observations",
      title = "Flowering Observations Across the Year"
    ) +
    theme_minimal()
  
  # Combine map + histogram
  combined <- map_plot / hist_plot + plot_layout(heights = c(2, 1))
  
  # Save to file
  ggsave(
    filename = sprintf("frames/frame_%03d.png", current_doy),
    plot = combined,
    width = 10,
    height = 9,
    dpi = 100
  )
})
