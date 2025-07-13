library(data.table)
library(ggplot2)
library(MCMCpack)
library(tidyverse)
library(lme4)


# load simu
source("/home/jt-miller/Gurlab/SuperBlooms/draft-code/simu-script-v3.R")
#year_types <- c(rep("High", 10))
# Create three simulations with same seed, vary the detection probs 
set.seed(22011)       
year_types <- c(rep("Low",3), rep("Medium",4), rep("High",3))
out_100 <- sim.superbloom(
  n_years = 10, # num of years to simulate
  year_type_sequence = year_types, # vector of year types, e.g. c("Low", "Medium", "High") for each year
  days = 1:150, # interval of days to simulate
  mu_peak = 100, # average peak phenology 
  dur = 12, # duration of flowering period
  y_var = 10, # variance between years
  sigma_N = 0.5, # yearly variation in total abundance of flowering plants
  alpha_0 = 9, # base abundance for flowering plants
  alpha_1 = 2, # approx 7x abd of flowering plants during superbloom
  yearly_abd_varation = TRUE, # whether to allow yearly variation in abundance
  pbloom = 0.4, # probability of superbloom occurring during a Bernoulli trial
  pdetwend = 1, # probability of detection on weekends
  pdetwday = 1, # probability of detection on weekdays
  min_det = 1, # starting detection probability for beginning years in simu
  max_det = 1, # ending detection probability for end years in simu
  grid_n = 5, # number of grid cells by dim, so 5x5 = 25 cells
  dirichlet_conc = 1e6, # spatial cell heterozygosity (concentration parameter for Dirichlet distribution)
  site_heterogeniety = TRUE, # spatial heterogeneity of detection,
  plot_detection_trend = TRUE, # whether to plot a detection trend 
  plot_spatial_hotspots = TRUE, # whether to plot spatial hotspots 
  equal_site_detection = TRUE,  # whether to use equal detection across sites
  sdet = 1, # spatial cell detection prob, only applies if equal_site_detection is TRUE
  high_sdet = 1, # high site detection probability (skews when equal_site_detection is FALSE)        
  low_sdet = 1,  # low site detection probability (skews when equal_site_detection is FALSE)          
  prop_high_sdet = 0.5 # proportion of high detection sites (skews when equal_site_detection is FALSE)         
)
set.seed(22011)       
out_50 <- sim.superbloom(
  n_years = 10, # num of years to simulate
  year_type_sequence = year_types, # vector of year types, e.g. c("Low", "Medium", "High") for each year
  days = 1:150, # interval of days to simulate
  mu_peak = 100, # average peak phenology 
  dur = 12, # duration of flowering period
  y_var = 10, # variance between years
  sigma_N = 0.5, # yearly variation in total abundance of flowering plants
  alpha_0 = 9, # base abundance for flowering plants
  alpha_1 = 2, # approx 7x abd of flowering plants during superbloom
  yearly_abd_varation = TRUE, # whether to allow yearly variation in abundance
  pbloom = 0.4, # probability of superbloom occurring during a Bernoulli trial
  pdetwend = 0.5, # probability of detection on weekends
  pdetwday = 0.5, # probability of detection on weekdays
  min_det = 1, # starting detection probability for beginning years in simu
  max_det = 1, # ending detection probability for end years in simu
  grid_n = 5, # number of grid cells by dim, so 5x5 = 25 cells
  dirichlet_conc = 1e6, # spatial cell heterozygosity (concentration parameter for Dirichlet distribution)
  site_heterogeniety = TRUE, # spatial heterogeneity of detection,
  plot_detection_trend = TRUE, # whether to plot a detection trend 
  plot_spatial_hotspots = TRUE, # whether to plot spatial hotspots 
  equal_site_detection = TRUE,  # whether to use equal detection across sites
  sdet = 1, # spatial cell detection prob, only applies if equal_site_detection is TRUE
  high_sdet = 1, # high site detection probability (skews when equal_site_detection is FALSE)        
  low_sdet = 1,  # low site detection probability (skews when equal_site_detection is FALSE)          
  prop_high_sdet = 0.5 # proportion of high detection sites (skews when equal_site_detection is FALSE)         
)
set.seed(22011)       
out_0 <- sim.superbloom(
  n_years = 10, # num of years to simulate
  year_type_sequence = year_types, # vector of year types, e.g. c("Low", "Medium", "High") for each year
  days = 1:150, # interval of days to simulate
  mu_peak = 100, # average peak phenology 
  dur = 12, # duration of flowering period
  y_var = 10, # variance between years
  sigma_N = 0.5, # yearly variation in total abundance of flowering plants
  alpha_0 = 9, # base abundance for flowering plants
  alpha_1 = 2, # approx 7x abd of flowering plants during superbloom
  yearly_abd_varation = TRUE, # whether to allow yearly variation in abundance
  pbloom = 0.4, # probability of superbloom occurring during a Bernoulli trial
  pdetwend = 0, # probability of detection on weekends
  pdetwday = 0, # probability of detection on weekdays
  min_det = 1, # starting detection probability for beginning years in simu
  max_det = 1, # ending detection probability for end years in simu
  grid_n = 5, # number of grid cells by dim, so 5x5 = 25 cells
  dirichlet_conc = 1e6, # spatial cell heterozygosity (concentration parameter for Dirichlet distribution)
  site_heterogeniety = TRUE, # spatial heterogeneity of detection,
  plot_detection_trend = TRUE, # whether to plot a detection trend 
  plot_spatial_hotspots = TRUE, # whether to plot spatial hotspots 
  equal_site_detection = TRUE,  # whether to use equal detection across sites
  sdet = 1, # spatial cell detection prob, only applies if equal_site_detection is TRUE
  high_sdet = 1, # high site detection probability (skews when equal_site_detection is FALSE)        
  low_sdet = 1,  # low site detection probability (skews when equal_site_detection is FALSE)          
  prop_high_sdet = 0.5 # proportion of high detection sites (skews when equal_site_detection is FALSE)         
)

analyze_simus <- function(out, predictor_var){ # predictor_var needs to be either flowering_obs (Simulation Truth, should be correct) or predictor_var
  # Yearly-Site Analysis 
  ## create summary of flowering proportions by year and site.
  siteyear_features <- out %>%
    group_by(cell_id, year) %>%
    summarize(
      total = sum({{predictor_var}}),
      mean_obs = mean({{predictor_var}}),
      S_y_site = unique(S_y_site),
      .groups = 'drop'
    )
  ## ensure cell_id and year are factors
  siteyear_features$cell_id <- as.factor(siteyear_features$cell_id)
  siteyear_features$year <- as.factor(siteyear_features$year)
  ## construct a mixed-effects model predicting site S_y per year given the mean obs as fixed effect & sites as a random effect
  mod <- glmer(
    S_y_site ~ mean_obs + (1|cell_id), 
    data = siteyear_features,
    family = binomial(link = "logit")
  )
  ## attach probs & predict with a generic 0.5 threshold cutoff
  siteyear_features$superbloom_prob_glmm <- predict(mod, type = "response")
  siteyear_features$S_y_site_predict <- ifelse(siteyear_features$superbloom_prob_glmm > 0.5, 1, 0)
  confusion_matrix <- table(
    True = siteyear_features$S_y_site,
    Predicted = siteyear_features$S_y_site_predict
  )
  print(confusion_matrix)
  ggplot(siteyear_features, aes(x = as.factor(S_y_site), y = superbloom_prob_glmm)) +
    geom_boxplot() +
    labs(x = "True Superbloom Site-Year (S_y)", y = "Predicted Probability") + 
    ggtitle(label = "Superbloom Prediction by Site-Year") +
    theme_bw()
  
  ## Weekly-Year-Site Analysis
  ## create summary of flowering proportions by site-week-year
  siteweekyear_features <- out %>%
    mutate(
      week = ceiling(doy/7)
    ) %>%
    group_by(cell_id, year, week) %>%
    summarize(
      total = sum({{predictor_var}}),
      mean_obs = mean({{predictor_var}}),
      max_obs = max({{predictor_var}}),
      S_y_site = unique(S_y_site),
      .groups = 'drop'
    )
  ## ensure cell_id and year are factors
  siteweekyear_features$cell_id <- as.factor(siteweekyear_features$cell_id)
  siteweekyear_features$year <- as.factor(siteweekyear_features$year)
  ## construct a mixed-effects model predicting site S_y per week-year given the mean obs as fixed effect & sites as a random effect
  mod <- glmer(
    S_y_site ~ mean_obs + max_obs + (1|cell_id) + (1:year), 
    data = siteweekyear_features,
    family = binomial(link = "logit")
  )
  ## attach probs & predict with a generic 0.5 threshold cutoff
  siteweekyear_features$superbloom_prob_glmm <- predict(mod, type = "response")
  siteweekyear_features$S_y_site_predict <- ifelse(siteweekyear_features$superbloom_prob_glmm > 0.5, 1, 0)
  confusion_matrix <- table(
    True = siteweekyear_features$S_y_site,
    Predicted = siteweekyear_features$S_y_site_predict
  )
  print(confusion_matrix)
  ggplot(siteweekyear_features, aes(x = as.factor(S_y_site), y = superbloom_prob_glmm)) +
    geom_boxplot() +
    labs(x = "True Superbloom Site-Week-Year (S_y)", y = "Predicted Probability") + 
    ggtitle(label = "Superbloom Prediction by Site-Week-Year") +
    theme_bw()
  
  
}

### Non function v. for testing particular fits (Flowering Obs Simu TRuth)
# Yearly-Site Analysis 
## create summary of flowering proportions by year and site.
siteyear_features <- out %>%
  group_by(cell_id, year) %>%
  summarize(
    total = sum(flowering_obs),
    mean_obs = mean(flowering_obs),
    S_y_site = unique(S_y_site),
    .groups = 'drop'
  )
## ensure cell_id and year are factors
siteyear_features$cell_id <- as.factor(siteyear_features$cell_id)
siteyear_features$year <- as.factor(siteyear_features$year)
## construct a mixed-effects model predicting site S_y per year given the mean obs as fixed effect & sites as a random effect
mod <- glmer(
  S_y_site ~ mean_obs + (1|cell_id), 
  data = siteyear_features,
  family = binomial(link = "logit")#, 
  #control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))
)
## attach probs & predict with a generic 0.5 threshold cutoff
siteyear_features$superbloom_prob_glmm <- predict(mod, type = "response")
siteyear_features$S_y_site_predict <- ifelse(siteyear_features$superbloom_prob_glmm > 0.5, 1, 0)
confusion_matrix <- table(
  True = siteyear_features$S_y_site,
  Predicted = siteyear_features$S_y_site_predict
)
print(confusion_matrix)
ggplot(siteyear_features, aes(x = as.factor(S_y_site), y = superbloom_prob_glmm)) +
  geom_boxplot() +
  labs(x = "True Superbloom Site-Year (S_y)", y = "Predicted Probability") + 
  ggtitle(label = "Superbloom Prediction by Site-Year") +
  theme_bw()

## Weekly-Year-Site Analysis
## create summary of flowering proportions by site-week-year
siteweekyear_features <- out %>%
  mutate(
    week = ceiling(doy/7)
  ) %>%
  group_by(cell_id, year, week) %>%
  summarize(
    total = sum(flowering_obs),
    mean_obs = mean(flowering_obs),
    max_obs = max(flowering_obs),
    S_y_site = unique(S_y_site),
    .groups = 'drop'
  )
## ensure cell_id and year are factors
siteweekyear_features$cell_id <- as.factor(siteweekyear_features$cell_id)
siteweekyear_features$year <- as.factor(siteweekyear_features$year)
## construct a mixed-effects model predicting site S_y per week-year given the mean obs as fixed effect & sites as a random effect
mod <- glmer(
  S_y_site ~ mean_obs + max_obs + (1|cell_id) + (1:year), 
  data = siteweekyear_features,
  family = binomial(link = "logit")
)
## attach probs & predict with a generic 0.5 threshold cutoff
siteweekyear_features$superbloom_prob_glmm <- predict(mod, type = "response")
siteweekyear_features$S_y_site_predict <- ifelse(siteweekyear_features$superbloom_prob_glmm > 0.5, 1, 0)
confusion_matrix <- table(
  True = siteweekyear_features$S_y_site,
  Predicted = siteweekyear_features$S_y_site_predict
)
print(confusion_matrix)
ggplot(siteweekyear_features, aes(x = as.factor(S_y_site), y = superbloom_prob_glmm)) +
  geom_boxplot() +
  labs(x = "True Superbloom Site-Week-Year (S_y)", y = "Predicted Probability") + 
  ggtitle(label = "Superbloom Prediction by Site-Week-Year") +
  theme_bw()