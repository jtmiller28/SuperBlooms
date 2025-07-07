### Bootstrapping Simulation Predictive Accuracy
### Author: JT Miller

# load packages
library(lme4)
library(tidyverse)
library(MCMCpack)
library(ggplot2)
source("/home/jt-miller/Gurlab/SuperBlooms/draft-code/simulation-script.R") # simulation written

# Conceptual framework 
## for each observation (cell x day x year), the detection probability = pdet_day x pdet_spatial x pdet_year
### How does accuracy of superbloom prediction depend on the combination of these three components?

# Structure of Bootstrap
## for each iteration of 10 replicates, Simulate 10 years of flowering, Aggregate weekly sums and compute a mean_prop_flowering (weighted by observed flowering / observed background + observed flowering) 
## fit a mixed effects model predicting S_y (superbloom status) based on mean_prop_flowering as the predictor and the random effects of spatial cell and year
## predict per year superbloom probabilities, cut off 0.5 (arbitrary)
## compute a confusion matrix
## extract per row detection components pdet_day, pdet_spatial, pdet_year, and the product total_detection_probability
## Lable whether the year was a superbloom or not, and whether the prediction was correct or not


# Number of bootstrap iterations
n_iter <- 10

# define parameter combos 
param_grid <- expand.grid(
  pdetwend = c(0.2, 0.5), #pdetwend = c(0.2, 0.5, 0.8),
  pdetwday = c(0.001, 0.2), #pdetwday = c(0.1, 0.3, 0.6),
  min_det = c(0.05, 0.4),
  max_det = c(0.4, 0.6),
  high_sdet = c(0.3, 0.5),
  low_sdet = c(0.05, 0.3),
  prop_high_sdet = c(0.3, 0.7)
)
# store results in a list
all_iterations <- list() 
# For each parameter combination
for(row in 1:nrow(param_grid)) {
  
  params <- param_grid[row, ]
  message("Param set ", row, " / ", nrow(param_grid))
  
  # For replicates per parameter set
  for(replicate in 1){
    message(" Replicate ", replicate)
    
    # REPEAT until at least one superbloom
    repeat{
      sim_data <- sim.latent.with.background.spatial.weekend.dropout(
        n_years = 10,
        days = 1:150,
        mu_peak = 100,
        dur = 12,
        y_var = 10,
        alpha_0 = 9,
        alpha_1 = 1,
        pbloom = 0.2,
        pdetwend = params$pdetwend,
        pdetwday = params$pdetwday,
        min_det = params$min_det,
        max_det = params$max_det,
        grid_n = 5,
        dirichlet_conc = 10,
        small_constant = 2,
        bAbd_scaler = 1.5,
        equal_site_detection = FALSE,
        high_sdet = params$high_sdet,
        low_sdet = params$low_sdet,
        prop_high_sdet = params$prop_high_sdet,
        densityBasedMask = TRUE,
        plot_detection_trend = FALSE,
        plot_spatial_hotspots = FALSE,
        plot_dropout_surface = FALSE
      )
      if(any(sim_data$S_y == 1)) break
    }
    
    # Aggregate weekly
    sim_weekly <- sim_data %>% 
      mutate(week = ceiling(doy / 7)) %>%
      group_by(year, week, S_y, cell_id) %>% 
      summarize(
        flowering_obs_weekly_sum = sum(flowering_obs, na.rm = TRUE),
        mean_prop_flowering = mean(prop_observed_count_flowering, na.rm = TRUE),
        pdet_day = mean(pdet_day, na.rm = TRUE),
        pdet_spatial = mean(pdet_spatial, na.rm = TRUE),
        pdet_year = mean(pdet_year, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Fit GLMM
    model <- glmer(
      S_y ~ mean_prop_flowering + (1 | cell_id) + (1 | year),
      data = sim_weekly,
      family = binomial(),
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))
    )
    
    # Predict
    sim_weekly$pred_prob <- predict(model, type = "response")
    
    # Summarize per year
    year_preds <- sim_weekly %>%
      group_by(year) %>%
      summarize(
        mean_prob = mean(pred_prob),
        true_S_y = unique(S_y)
      ) %>%
      mutate(predicted_S_y = ifelse(mean_prob > 0.5, 1, 0))
    
    # Join predictions
    sim_weekly <- left_join(sim_weekly, year_preds, by = "year")
    
    # Add parameter metadata
    sim_weekly <- sim_weekly %>%
      mutate(
        param_set = row,
        replicate = replicate,
        pdetwend = params$pdetwend,
        pdetwday = params$pdetwday,
        min_det = params$min_det,
        max_det = params$max_det,
        high_sdet = params$high_sdet,
        low_sdet = params$low_sdet,
        prop_high_sdet = params$prop_high_sdet
      )
    
    # Store
    all_iterations[[length(all_iterations) + 1]] <- sim_weekly
  }
}

# Combine everything
results_all <- bind_rows(all_iterations)

# Summarize accuracy per parameter combination
accuracy_summary <- results_all %>%
  group_by(param_set, pdetwend, pdetwday, min_det, max_det, high_sdet, low_sdet, prop_high_sdet) %>%
  summarize(
    accuracy = mean(predicted_S_y == true_S_y),
    .groups = "drop"
  )

ggplot(accuracy_summary, aes(x=pdetwend, y=high_sdet, fill=accuracy)) +
  geom_tile() +
  facet_grid(min_det ~ max_det) +
  scale_fill_viridis_c() +
  labs(
    title="Accuracy of Superbloom Prediction",
    x="Weekend Detection Probability",
    y="High Site Detection Probability",
    fill="Accuracy"
  ) +
  theme_minimal()