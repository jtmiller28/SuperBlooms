### 002 Observation Binning 
### Project: SuperBlooms
### Author: JT Miller
### Date: 03-31-2025


### Purpose: Using prior superblooms in 2016, 2019, and 2023-2024, visualize a how superblooms look according to the observation record

# load packages
library(ggplot2)
library(data.table)
library(tidyverse)

# load the subsetted data
flower_data <- fread("/blue/guralnick/millerjared/SuperBlooms/data/processed/desert-observations-hexed.csv")

# filter data to only include flowering records (remove once subset finishes.)
flower_data <- flower_data %>% filter(phenology_trait == "flower")

# grab observations in 2016, 2019, and 2023, construct bins: weekly, monthly
flower_data_sbs <- flower_data %>% 
  filter(year %in% c(2016, 2019, 2023))

flower_data_wsbs <- flower_data_sbs %>% # weekly bins
  mutate(week_bin = cut(doy, breaks = seq(1, 366, by = 7), right = FALSE, labels = FALSE)) %>% 
  filter(!is.na(week_bin)) %>% # 365 doesnt quite make it into these categories, so Im electing to throw them away for exploratory purposes

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

    

