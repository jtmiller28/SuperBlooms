### Suplement figures


# load libarries 
library(ggplot2)
library(data.table)
library(tidyverse)

flower_data <- data.table::fread("/blue/guralnick/millerjared/SuperBlooms/data/processed/desert-observations-hexed.csv")

# remove redundant annotated flowering data for unique observational records
flower_data <- dplyr::distinct(flower_data, obs_url, .keep_all = TRUE)

flower_data_temporalbins <- flower_data %>% # weekly bins
  dplyr::mutate(week_bin = cut(doy, breaks = seq(1, 366, by = 7), right = FALSE, labels = FALSE)) %>% # monthly bins
  dplyr::mutate(month_bin = cut(doy, breaks = seq(1,366, by = 30), right = FALSE, labels = FALSE)) %>% 
  dplyr::filter(!is.na(week_bin)) %>% # remove records that failed to be binned appropriately 
  dplyr::filter(!is.na(month_bin))

# remove years with low obs count
year_obs_check <- flower_data_temporalbins %>% group_by(year) %>% summarize(n = n()) %>% arrange(desc(n))
year_obs_500 <- year_obs_check %>% filter(n >= 500) # set an arbitrary req
flower_data_temporalbins <- flower_data_temporalbins %>% filter(year %in% year_obs_500$year) %>% filter(year != 2024) # also remove year 2024, as this data is truncated by the 12th week

# remove non-spring weeks
flower_data_temporalbins <- flower_data_temporalbins %>% 
  filter(week_bin <= 22)

flower_data_summaries <- flower_data_temporalbins %>% 
  group_by(hex50_id) %>% 
  summarize(count = n())

hex50 <- readRDS("./data/processed/50km-hexed-map.rds")

flower_sf_merge <- merge(hex50, flower_data_summaries, by = "hex50_id", all.x = TRUE)




ggplot() +
  geom_sf(flower_sf_merge, mapping = aes(fill = log(count+1))) +
  theme_bw() + 
  scale_fill_gradientn(
    colours = colorRampPalette((RColorBrewer::brewer.pal(11, "YlOrRd")))(9), 
    na.value = "#808080" # The line that denotes NAs as grey
  )  + 
  ggtitle("Observation Intensity per 50km grid cell \n Across the Mojave & Sonoran Deserts") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("./figures/obs-intensity-50km.png", height = 10, width = 11)
