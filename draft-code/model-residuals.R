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
library(gratia)
library(DHARMa)
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

wkbin_hex_year_summary2 <- mutate(wkbin_hex_year_summary2, hex50_id = as.factor(hex50_id)) %>% mutate(log_area = log(as.numeric(area) + 1))
# remove highly skewed low area cells from the dataset, as these are causing fititng issues.
wkbin_hex_year_summary3 <- wkbin_hex_year_summary2 %>% filter(as.numeric(area) >= 84230390)
mod <- gam(num_observations_flowering ~ 
             s(week_bin, bs = "cs") + 
             #s(hex50_id, bs = "re") + 
             s(year, bs = "cs", k = 7) +
             s(as.numeric(area), bs = "tp", k = 10) + # unsure about this, but I think area could be non-linear and it needs to be better at adapting to extremes so we're using thin plate splines
             log(num_observations_total), # removed offset in order to allow the model to learn from obs data 
           data = wkbin_hex_year_summary2, 
           family = "poisson")

mod_nb <- gam(num_observations_flowering ~ 
                s(week_bin, bs = "cs", k = 15) + 
                s(hex50_id, bs = "re") +
                s(year, bs = "cs", k = 7) +
                #s(as.numeric(area), bs = "tp", k = 109) +
                s(log_area, bs = "tp", k = 100) +
                log(num_observations_total),
              method = "REML",
              data = wkbin_hex_year_summary3, 
              family = nb())

gam.check(mod_nb)

mod <- mod_nb
summary(mod)
plot(mod, pages = 1, residuals = TRUE, shade = TRUE, seWithMean = TRUE)
library(gratia)
draw(mod, residuals = TRUE) # prob a more informative plot
png(filename = "./figures/gam-50km-wkbin-plots-noresiduals.png", width = 800, height = 600, res = 100)
draw(mod)

dev.off()
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

# draw a rootogram 
rootogram(mod, max = 600) %>% draw()
# Extract deviations as superbloom candidates
wkbin_hex_year_summary3$resid <- residuals(mod, type = "pearson")
wkbin_hex_year_summary3$fitted <- fitted(mod)

# ID large pos residuals (z-score > 2)
wkbin_hex_year_summary3$superBloomCandidate <- wkbin_hex_year_summary3$resid > 2

test <- filter(wkbin_hex_year_summary3, superBloomCandidate == TRUE)
test %>% group_by(year) %>% summarize(n = n()) %>% arrange(desc(n))
# we know that superblooms are a spring occurrence, therefore lets remove candidates that happen after the last week of May : doy 151/7 = 21.5 ~ 22
test_spring <- filter(test , week_bin <= 22)
test_spring %>% group_by(year) %>% summarize(n = n()) %>% arrange(desc(n))
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



### More scratch work #########################################################################
# annotated & hexbinned flowering data
flower_data <- data.table::fread("./data/processed/desert-observations-hexed.csv")
# load background observation data for vascular plants of these regions (hexbinned prior)
obs_data <- fread("./data/processed/tracheophyte-desert-inats.csv")
flower_data <- dplyr::distinct(flower_data, obs_url, .keep_all = TRUE)
obs_data <- obs_data %>% 
  mutate(eventDateParsed = parse_date_time(eventDate, orders = c("ymd_HMS", "ymd_HM", "ymd"))) %>% 
  mutate(doy = lubridate::yday(eventDateParsed)) %>% 
  mutate(year = lubridate::year(eventDateParsed)) %>% 
  mutate(month = lubridate::month(eventDateParsed))
# create temporal bins 
flower_data_temporalbins <- flower_data %>% 
  # weekly bins
  dplyr::mutate(week_bin = cut(doy, breaks = seq(1, 366, by = 7), right = FALSE, labels = FALSE)) %>%
  # monthly bins
  dplyr::mutate(month_bin = cut(doy, breaks = seq(1,366, by = 30), right = FALSE, labels = FALSE)) %>% 
  # remove tail end data 
  dplyr::filter(!is.na(week_bin)) %>% 
  dplyr::filter(!is.na(month_bin))

# for obs data 
obs_data_temporalbins <- obs_data %>% 
  # weekly bins
  dplyr::mutate(week_bin = cut(doy, breaks = seq(1, 366, by = 7), right = FALSE, labels = FALSE)) %>%
  # monthly bins
  dplyr::mutate(month_bin = cut(doy, breaks = seq(1,366, by = 30), right = FALSE, labels = FALSE)) %>% 
  # remove tail end data 
  dplyr::filter(!is.na(week_bin)) %>% 
  dplyr::filter(!is.na(month_bin))

# set up spatial & temporal bin retrieval vectors
hex_map_files <- list.files("./data/processed", pattern = "\\.rds$", full.names = TRUE)
hex_map_res <- c(100, 10, 25, 50, 5)
temporal_bins <- c("doy", "week_bin", "month_bin")

# edit depending on what im curious to test. 
hex_grab <- 4
hex_cell_size <- hex_map_res[hex_grab]
temporal_name <- temporal_bins[2]
# create group by vars
hex_name <- paste0("hex", hex_cell_size, "_id")
group_col <- list(sym("year"), sym(hex_name), sym(temporal_name))
# assign bins
hex_flower_data_bin <- flower_data_temporalbins %>% 
  group_by(!!!group_col) %>% # !!! splice operator
  summarize( 
    num_observations_flowering = n(), 
    uniqueFloweringSpecies = n_distinct(taxon_name),
    .groups = "drop")

hex_obs_data_bin <- obs_data_temporalbins %>% 
  group_by(!!!group_col) %>% # use !!! splice operator     
  summarize(
    num_observations_total = n(), 
    uniqueObsSpecies = n_distinct(verbatimScientificName), 
    .groups = "drop"
  )
# merge flowreing and obs data 
hex_flower_obs_bin <- merge(hex_flower_data_bin, hex_obs_data_bin, by = c("year", hex_name, temporal_name), all.x = TRUE)
# remove data that has incorrect background
hex_flower_obs_bin <- hex_flower_obs_bin %>%       
  filter(!is.na(num_observations_total)) %>% # remove na obs
  mutate(obsFlowerRatioCheck = ifelse(num_observations_flowering <= num_observations_total, TRUE, FALSE)) %>% 
  filter(obsFlowerRatioCheck == TRUE)
# add area as a var 
hex_grid <- read_rds(paste0(hex_map_files[hex_grab]))
hex_info <- select(hex_grid, hex_name, area) %>% sf::st_drop_geometry()
hex_flower_obs_bin <- merge(hex_info, hex_flower_obs_bin , by = hex_name)
hex_flower_obs_bin <- hex_flower_obs_bin %>% 
  mutate(hex50_id = as.factor(hex50_id), 
         log_area = log(as.numeric(area) + 1))
# remove data from years of low data sampling
year_obs_check <- hex_flower_obs_bin %>% group_by(year) %>% summarize(n = sum(num_observations_flowering)) %>% arrange(desc(n))
year_obs_500 <- year_obs_check %>% filter(n >= 100) # set an arbitrary req
hex_flower_obs_bin <- hex_flower_obs_bin %>% filter(year %in% year_obs_500$year) %>% filter(year != 2024) # also remove year 2024, as this data is truncated by the 12th week
#hex_flower_obs_bin <- hex_flower_obs_bin %>% filter(num_observations_total != 1) # attempt removal of data (1s are wierdly skewing the data)
# Attempt to fit a GAM 
knots <- list(week_bins = c(0.5, 52.5))
mod_nb <- bam(num_observations_flowering ~ # Y (response var) is number of flowering records obs
                te(year, week_bin, bs = c("tp", "cc"), k = 10) + # for seasonal events, if we want them to vary overtime.
                s(hex50_id, bs = "re") + # code hexcell id as a random effect
                s(log_area, bs = "tp", k = 10) + # log of area with a thin-plate spline 
                log(num_observations_total),
              knots = knots, # where the cyclical spline should bridge data
              method = "fREML",
              data = hex_flower_obs_bin, 
              family = nb())
set.seed(100)
mod_nb <- gam(num_observations_flowering ~ # Y (response var) is number of flowering records obs
                s(year) + # overall trend
                s(week_bin, bs = "cc") + # seasonal effect
                # use tensor product to gather up interaction term?
                ti(year, week_bin, bs = c("tp", "cc"), k = 17) + # for seasonal events, if we want to make them interact and vary overtime.
                s(hex50_id, bs = "re") + # code hexcell id as a random effect
                s(log_area, bs = "tp", k = 10) + # log of area with a thin-plate spline 
                log(num_observations_total),
              knots = knots, # where the cyclical spline should bridge data
              method = "REML",
              data = hex_flower_obs_bin, 
              family = nb())

# check model results 
summary(mod_nb)

# check whether k-value basis needs adj
gam.check(mod_nb)

# draw gam
draw(mod_nb)
rootogram(mod_nb, max = 300) %>% draw()
# use DARHMa to check for regular glm assumptions being met
sim_res <- DHARMa::simulateResiduals(mod_nb, plot = TRUE)
sim_test_outliers <- testOutliers(sim_res, type = "bootstrap")
# check for overdispersion, especially if the test above is sig
testDispersion(sim_res)

# use a residuals test to pull out and label candidate superblooms
hex_flower_obs_bin$resid <- residuals(mod_nb, type = "pearson")
hex_flower_obs_bin$fitted <- fitted(mod_nb)

# ID large pos residuals (z-score > 2)
hex_flower_obs_bin$superBloomCandidate <- hex_flower_obs_bin$resid > 2

possibleSuperblooms <- filter(hex_flower_obs_bin, superBloomCandidate == TRUE)
possibleSuperblooms %>% group_by(year) %>% summarize(n = n()) %>% arrange(desc(n))
# we know that superblooms are a spring occurrence, therefore lets remove candidates that happen after the last week of May : doy 151/7 = 21.5 ~ 22
possibleSuperblooms_springOnly<- filter(possibleSuperblooms , week_bin <= 22)
possibleSuperblooms_springOnly %>% group_by(year) %>% summarize(n = n()) %>% arrange(desc(n))
fwrite(possibleSuperblooms_springOnly, "/blue/guralnick/millerjared/SuperBlooms/data/processed/possibleCandidateSuperblooms.csv")
##############################################################################################




















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

