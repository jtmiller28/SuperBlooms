### Model Residuals 
### Project: Superblooms
### Author: JT Miller
### Date: 04-09-2025

## Here we're going to try a modeling approach for detecting large positive residuals 
### I'm going to take two approaches, a) Will be using General Additive Models to create models of *expected* flowering, then a Pearsons residual test we will determine whether large positive residuals exist during key superbloom years.
### b) Bayesian Hierarchical Models - TBD

### Load libraries 
library(data.table)
library(tidyverse)

### Load the hexed flowering obs data
flower_data <- fread("./data/processed/desert-observations-hexed.csv")
# we only want unique observations (currently the data is structured so that multiple images of one obs event take up uniq rows)
flower_data <- dplyr::distinct(flower_data, obs_url, .keep_all = TRUE)

# bin the data into weekly bins summarizing counts
flower_data_wkbin<- flower_data %>% # weekly bins
  mutate(week_bin = cut(doy, breaks = seq(1, 366, by = 7), right = FALSE, labels = FALSE)) %>% 
  filter(!is.na(week_bin)) # 365 doesnt quite make it into these categories, so Im electing to throw them away for exploratory purposes

# sum up number of obs per year-hexcell-wkbin 
wkbin_hex_year_summary <- flower_data_wkbin %>% 
  group_by(year, hex50_id, week_bin) %>% 
  summarize(
    num_observations = n(), 
    uniqueSpecies = n_distinct(taxon_name),
    .groups = "drop"
  )
