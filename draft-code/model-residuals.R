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
library(lubridate)
library(mgcv)

### Load the hexed flowering obs data
flower_data <- fread("./data/processed/desert-observations-hexed.csv")
# we only want unique observations (currently the data is structured so that multiple images of one obs event take up uniq rows)
flower_data <- dplyr::distinct(flower_data, obs_url, .keep_all = TRUE)

### Load the hexed gbif RG inaturalist data 
obs_data <- fread("./data/processed/tracheophyte-desert-inats.csv")
obs_data <- obs_data %>% 
  mutate(eventDateParsed = parse_date_time(eventDate, orders = c("ymd_HMS", "ymd_HM", "ymd"))) %>% 
  mutate(doy = yday(eventDateParsed)) %>% 
  mutate(year = year(eventDateParsed))

# bin the data into weekly bins summarizing counts
flower_data_wkbin<- flower_data %>% # weekly bins
  mutate(week_bin = cut(doy, breaks = seq(1, 366, by = 7), right = FALSE, labels = FALSE)) %>% 
  filter(!is.na(week_bin)) %>% # 365 doesnt quite make it into these categories, so Im electing to throw them away for exploratory purposes
  filter(year %in% c(2017, 2018, 2019, 2020, 2021, 2022, 2023)) # years of high quantity data

obs_data_wkbin <- obs_data %>% 
  mutate(week_bin = cut(doy, breaks = seq(1, 366, by =7), right = FALSE, labels = FALSE)) %>% 
  filter(!is.na(week_bin))  %>% 
  filter(year %in% c(2017, 2018, 2019, 2020, 2021, 2022, 2023))

# sum up number of obs per year-hexcell-wkbin 
wkbin_flowering_hex_year_summary <- flower_data_wkbin %>% 
  group_by(year, hex50_id, week_bin) %>% 
  summarize(
    num_observations_flowering = n(), 
    uniqueFloweringSpecies = n_distinct(taxon_name),
    .groups = "drop"
  )

wkbin_obs_hex_year_summary <- obs_data_wkbin %>% 
  group_by(year, hex50_id, week_bin) %>% 
  summarize(
    num_observations_total = n(), 
    uniqueObsSpecies = n_distinct(verbatimScientificName), 
    .groups = "drop"
  )

# merge these two table
wkbin_hex_year_summary <- merge(wkbin_flowering_hex_year_summary, wkbin_obs_hex_year_summary, by = c("year", "hex50_id", "week_bin"), all.x = TRUE)

# as we are using inat RG records as a proxy, we're missing some data (about 200 year-hex-weekbins). Throw out wkbin-hex-years that are NA for modeling purposes
wkbin_hex_year_summary <- wkbin_hex_year_summary %>% filter(!is.na(num_observations_total))

# as a check, make sure theres ALWAYS more obs then flowering records 
wkbin_hex_year_summary <- wkbin_hex_year_summary %>% 
  mutate(obsFlowerRatioCheck = ifelse(num_observations_flowering <= num_observations_total, TRUE, FALSE))

# remove these, as they would mess with our model
wkbin_hex_year_summary <- wkbin_hex_year_summary %>% 
  filter(obsFlowerRatioCheck == TRUE)

# affix grid size to these data, so that we can use it as a fixed effect for the model 
hex_50km <- read_rds("/home/jt-miller/Gurlab/SuperBlooms/data/processed/50km-hexed-map.rds")
hex_info <- select(hex_50km, hex50_id, area) %>% sf::st_drop_geometry()
wkbin_hex_year_summary2 <- merge(hex_info, wkbin_hex_year_summary, by = "hex50_id")

## Given the week of the year, the hex cell, year, and how many total iNat obs were made, how many flowering records should I expect? 
mod <- mgcv::gam(num_observations_flowering ~
                 s(week_bin, bs = "cs") + # choosing cubic as temporal effects need more smoothing
                 s(hex50_id, bs = "re") + # choosing random effect as we want each bin to have random intercepts
                 s(year, bs = "cs", k = 7) +# smoothung trends across years in order to account for increased num of records per year, k = 8 as default of 10 wont fit
                 offset(log(num_observations_total)), # log of effort to control for more users and observations in later years
                 data = wkbin_hex_year_summary, 
                 family = "poisson")

mod <- gam(num_observations_flowering ~ 
             s(week_bin, bs = "cs") + 
             s(hex50_id, bs = "re") + 
             s(year, bs = "cs", k = 7) +
             s(as.numeric(area), bs = "tp") + # unsure about this, but I think area could be non-linear and it needs to be better at adapting to extremes so we're using thin plate splines
             log(num_observations_total), # removed offset in order to allow the model to learn from obs data 
           data = wkbin_hex_year_summary2, 
           family = "nb")
summary(mod)
plot(mod, pages = 1, residuals = TRUE, shade = TRUE, seWithMean = TRUE)
library(gratia)
draw(mod, residuals = TRUE) # prob a more informative plot

residuals_raw <- residuals(mod, type = "pearson")
hist(residuals_raw, breaks = 30, main = "Pearsons Risiduals Histogram")

plot(fitted(mod), residuals_raw,
     xlab = "Fitted values", ylab = "Pearson residuals",
     main = "Residuals vs. Fitted")
abline(h = 0, col = "red")

# simulation of residuals 
library(DHARMa)
sim_res <- simulateResiduals(mod, plot = TRUE)
sim_test_outliers <- testOutliers(sim_res, type = "bootstrap")
# check for overdispersion, especially if the test above is sig
testDispersion(sim_res)

# Extract deviations as superbloom candidates
wkbin_hex_year_summary$resid <- residuals(mod, type = "pearson")
wkbin_hex_year_summary$fitted <- fitted(mod)

# ID large pos residuals (z-score > 2)
wkbin_hex_year_summary$superBloomCandidate <- wkbin_hex_year_summary$resid > 2

test <- filter(wkbin_hex_year_summary, superBloomCandidate == TRUE)

# visualize the estimated effects of vars on the model 
summary(mod)
draw(mod, select = 1)
draw(mod, select = 1, residuals = TRUE)
draw(mod, select = 2)
draw(mod, select = 2, residuals = TRUE)
draw(mod, select = 3)
draw(mod, select = 3, residuals = TRUE)
draw(mod, select = 4)
draw(mod, select = 4, residuals = TRUE)
# build a plot IDing superblooms and plots where superblooms occurred in 2019
library(ggplot2)
# load in gridded map of study region 
hex_50km <- read_rds("/home/jt-miller/Gurlab/SuperBlooms/data/processed/50km-hexed-map.rds")
hex_50_wk_tbl <- merge(hex_50km, test, by = "hex50_id", all.y = TRUE)
# filter to just 2019 for simplicity 
hex_50_wk_tbl2 <- hex_50_wk_tbl %>% 
  filter(year == 2023)

ggplot() +
  geom_sf(hex_50km, mapping = aes()) + 
  geom_sf(hex_50_wk_tbl2, mapping = aes(), fill = "goldenrod") +
  ggtitle("Cells with candidate superblooms")
# sum up the total obs across all spatial cells to build a histogram of phenology across the deserts (for histogram visual)
hex_50_wk_hist_summary <- hex_50_wk_tbl %>% 
  sf::st_drop_geometry() %>% 
  dplyr::group_by(week_bin) %>% 
  dplyr::summarise(desert_week_bin_flowering_obs = sum(num_observations_flowering, na.rm = TRUE), 
                   superBloomCandidate )

