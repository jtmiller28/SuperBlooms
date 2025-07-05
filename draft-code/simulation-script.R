# simu.sb <- function(
#     n_years = 30,
#     days = 1:365,
#     mu_peak = 100,
#     dur = 12,
#     y_var = 0,
#     alpha_0 = 9,
#     alpha_1 = 2.5,
#     pbloom = 0.2,
#     pdetwend = 0.8,
#     pdetwday = 0.5,
#     beta_0 = 3,
#     beta_1 = 0.5,
#     beta_2 = 0.3,
#     sigma_bg = 20,
#     densityObservable = 10,
#     dropout_prob_densityD = 0.3,
#     grid_n = 5,
#     sdet = 0.6,
#     dirichlet_conc = 1e6
# ) {
#   # Temporal hyperparameters
#   mu_mu <- rnorm(1, mean = mu_peak, sd = 10)
#   sigma_mu <- abs(rnorm(1, mean = y_var, sd = 5))
#   mu_sigma <- rnorm(1, mean = dur, sd = 3)
#   sigma_sigma <- abs(rnorm(1, mean = 0, sd = 3))
#   
#   # Superbloom probability
#   calibrate_bloom_beta <- function(mean, fixed_shape1 = 2) {
#     shape1 <- fixed_shape1
#     shape2 <- shape1 * (1 - mean) / mean
#     list(pshape1 = shape1, pshape2 = shape2)
#   }
#   bloom_params <- calibrate_bloom_beta(pbloom)
#   p_bloom <- rbeta(1, bloom_params$pshape1, bloom_params$pshape2)
#   
#   # Year-specific flowering parameters
#   mu_y <- rnorm(n_years, mu_mu, sigma_mu)
#   sigma_y <- abs(rnorm(n_years, mu_sigma, sigma_sigma))
#   S_y <- rbinom(n_years, 1, p_bloom)
#   log_N_y <- alpha_0 + alpha_1 * S_y
#   N_y <- exp(log_N_y)
#   
#   # Spatial detection heterogeneity
#   calibrate_beta <- function(mean, fixed_shape1 = 10) {
#     shape1 <- fixed_shape1
#     shape2 <- shape1 * (1 - mean) / mean
#     list(shape1 = shape1, shape2 = shape2)
#   }
#   beta_params <- calibrate_beta(sdet, fixed_shape1 = 2)
#   spatial_grid <- expand.grid(x = 1:grid_n, y = 1:grid_n) %>%
#     mutate(
#       cell_id = paste0("cell_", x, "_", y),
#       pdet_spatial = rbeta(n(), beta_params$shape1, beta_params$shape2)
#     )
#   n_cells <- nrow(spatial_grid)
#   
#   # Weekday/weekend detection
#   is_weekend <- function(doy) {
#     weekday <- (doy - 1) %% 7 + 1
#     weekday %in% c(6,7)
#   }
#   weekend_flag <- is_weekend(days)
#   pdet_day_vec <- ifelse(weekend_flag, pdetwend, pdetwday)
#   
#   # Generate data for each year
#   sim_obs <- purrr::map_dfr(1:n_years, function(y) {
#     doy_density <- dnorm(days, mu_y[y], sigma_y[y])
#     
#     # Dirichlet proportions across cells
#     prop_cells <- MCMCpack::rdirichlet(1, rep(dirichlet_conc, n_cells))[1,]
#     
#     # for each cell
#     daily_results <- purrr::map_dfr(1:n_cells, function(i) {
#       cell <- spatial_grid[i,]
#       prop_c <- prop_cells[i]
#       
#       # Expected flowering abundance
#       A_d_c <- N_y[y] * prop_c * doy_density
#       true_flowering <- rpois(length(A_d_c), A_d_c)
#       
#       # Observation probability (weekday/weekend * spatial)
#       pdet_day_cell <- pdet_day_vec * cell$pdet_spatial
#       obs_flowering <- rbinom(length(A_d_c), size = true_flowering, prob = pdet_day_cell)
#       
#       # Randomly zero out some counts (density-dependent dropout)
#       dropout_prob <- ifelse(obs_flowering <= densityObservable, 0.5, dropout_prob_densityD)
#       mask <- rbinom(length(obs_flowering), size=1, prob=1 - dropout_prob)
#       obs_flowering <- obs_flowering * mask
#       
#       # Background abundance (same across all cells for simplicity)
#       lambda_bg <- exp(
#         beta_0 +
#           beta_1 * 100 * dnorm(days, mu_y[y], sigma_bg) +
#           beta_2 * S_y[y]
#       )
#       bg_counts <- rpois(length(lambda_bg), lambda_bg)
#       
#       # Combine counts
#       total_obs <- obs_flowering + bg_counts
#       
#       tibble(
#         year = y,
#         doy = days,
#         S_y = S_y[y],
#         mu_y = mu_y[y],
#         sigma_y = sigma_y[y],
#         N_y = N_y[y],
#         abundance = A_d_c,
#         true_flowering = true_flowering,
#         flowering_obs = obs_flowering,
#         background_obs = bg_counts,
#         total_obs = total_obs,
#         prop_flowering = ifelse(total_obs > 0, obs_flowering/total_obs, 0),
#         x = cell$x,
#         y = cell$y,
#         cell_id = cell$cell_id,
#         pdet_day = pdet_day_cell,
#         prop_abundance = prop_c
#       )
#     })
#     # Compute total true count per year & doy
#     totals <- daily_results %>%
#       group_by(year, doy) %>%
#       summarize(
#         true_flowering_total = sum(true_flowering),
#         abundance_total = sum(abundance),
#         .groups = "drop"
#       )
#     
#     # Join back
#     daily_results %>%
#       left_join(totals, by = c("year", "doy"))
#     
#   })
#   return(sim_obs)
#   
# }

sim.latent.with.background.spatial.weekend.dropout <- function(
    n_years = 30, # num of years to simulate
    days = 1:365, # interval of days to simulate
    mu_peak = 100, # average peak phenology 
    dur = 12, # duration of flowering period
    y_var = 0, # variance between years
    alpha_0 = 9, # base abundance for flowering plants
    alpha_1 = 2.5, # additional superbloom effect on abundance, 2.5 = 12.5x (too high for most simus)
    pbloom = 0.2, # probability of superbloom occurring during a Bernoulli trial
    pdetwend = 0.8, # probability of detection on weekends
    pdetwday = 0.5, # probability of detection on weekdays
    min_det = 0.2, # starting detection probability for beginning years in simu
    max_det = 0.8, # ending detection probability for end years in simu
    densityObservable = 10, # threshold for density-dependent dropout (how many obs do you need to have to start recording)
    dropout_prob_densityD = 0.3, # probability that we drop it if this density dependence isnt met
    grid_n = 5, # number of grid cells by dim, so 5x5 = 25 cells
    sdet = 0.6, # spatial cell detection prob, only applies if equal_site_detection is TRUE
    dirichlet_conc = 1e6, # spatial cell heterozygosity (concentration parameter for Dirichlet distribution)
    bAbd_scaler = 0.2, # scaler value for background abundance 
    small_constant = 3, # intial value for background abundance 
    plot_detection_trend = FALSE, # whether to plot a detection trend 
    plot_spatial_hotspots = FALSE, # whether to plot spatial hotspots 
    equal_site_detection = TRUE,  # whether to use equal detection across sites
    high_sdet = 0.9, # high site detection probability (skews when equal_site_detection is FALSE)        
    low_sdet = 0.3,  # low site detection probability (skews when equal_site_detection is FALSE)          
    prop_high_sdet = 0.5 # proportion of high detection sites (skews when equal_site_detection is FALSE)         
) {
  # Temporal hyperparameters
  mu_mu <- rnorm(1, mean = mu_peak, sd = 10)
  sigma_mu <- abs(rnorm(1, mean = y_var, sd = 5))
  mu_sigma <- rnorm(1, mean = dur, sd = 3)
  sigma_sigma <- abs(rnorm(1, mean = 0, sd = 3))
  
  # Superbloom probability
  calibrate_bloom_beta <- function(mean, fixed_shape1 = 2) {
    shape1 <- fixed_shape1
    shape2 <- shape1 * (1 - mean) / mean
    list(pshape1 = shape1, pshape2 = shape2)
  }
  bloom_params <- calibrate_bloom_beta(pbloom)
  p_bloom <- rbeta(1, bloom_params$pshape1, bloom_params$pshape2)
  
  # Year-specific flowering parameters
  mu_y <- rnorm(n_years, mu_mu, sigma_mu)
  sigma_y <- abs(rnorm(n_years, mu_sigma, sigma_sigma))
  S_y <- rbinom(n_years, 1, p_bloom)
  log_N_y <- alpha_0 + alpha_1 * S_y
  N_y <- exp(log_N_y)
  
  # Spatial detection heterogeneity
  calibrate_beta <- function(mean, fixed_shape1 = 10) {
    shape1 <- fixed_shape1
    shape2 <- shape1 * (1 - mean) / mean
    list(shape1 = shape1, shape2 = shape2)
  }
  
  spatial_grid <- expand.grid(x = 1:grid_n, y = 1:grid_n) %>%
    dplyr::mutate(cell_id = paste0("cell_", x, "_", y))
  n_cells <- nrow(spatial_grid)
  
  if (equal_site_detection) {
    # All cells have the same detection centered at sdet
    beta_params <- calibrate_beta(sdet, fixed_shape1 = 2)
    spatial_grid <- spatial_grid %>%
      dplyr::mutate(
        pdet_spatial = rbeta(n(), beta_params$shape1, beta_params$shape2)
      )
  } else {
    # Split cells into high and low detection groups
    n_high <- round(prop_high_sdet * n_cells)
    n_low <- n_cells - n_high
    high_cells <- sample(1:n_cells, n_high)
    spatial_grid <- spatial_grid %>%
      dplyr::mutate(
        detection_group = ifelse(dplyr::row_number() %in% high_cells, "high", "low"),
        pdet_spatial = dplyr::case_when(
          detection_group == "high" ~ rbeta(1, calibrate_beta(high_sdet)$shape1, calibrate_beta(high_sdet)$shape2),
          detection_group == "low" ~ rbeta(1, calibrate_beta(low_sdet)$shape1, calibrate_beta(low_sdet)$shape2)
        )
      )
  }
  
  # Weekday/weekend detection
  is_weekend <- function(doy) {
    weekday <- (doy - 1) %% 7 + 1
    weekday %in% c(6, 7)
  }
  weekend_flag <- is_weekend(days)
  pdet_day_vec <- ifelse(weekend_flag, pdetwend, pdetwday)
  
  # Detection trend over years
  pdet_trend <- seq(min_det, max_det, length.out = n_years)
  
  # Generate data for each year
  sim_obs <- purrr::map_dfr(1:n_years, function(y) {
    doy_density <- dnorm(days, mu_y[y], sigma_y[y])
    prop_cells <- MCMCpack::rdirichlet(1, rep(dirichlet_conc, n_cells))[1,]
    
    daily_results <- purrr::map_dfr(1:n_cells, function(i) {
      cell <- spatial_grid[i, ]
      prop_c <- prop_cells[i]
      A_d_c <- N_y[y] * prop_c * doy_density
      true_flowering <- rpois(length(A_d_c), A_d_c)
      #pdet_day_cell <- pdet_day_vec * cell$pdet_spatial * pdet_trend[y]
      pdet_day_component <- pdet_day_vec
      pdet_spatial_component <- cell$pdet_spatial
      pdet_year_component <- pdet_trend[y]
      total_pdet <- pdet_day_component * pdet_spatial_component * pdet_year_component
      
      obs_flowering <- rbinom(length(A_d_c), size = true_flowering, prob = total_pdet)
      dropout_prob <- ifelse(obs_flowering <= densityObservable, 0.5, dropout_prob_densityD)
      mask <- rbinom(length(obs_flowering), size = 1, prob = 1 - dropout_prob)
      obs_flowering <- obs_flowering * mask
      
      lambda_bg <- small_constant + bAbd_scaler * A_d_c
      bg_true_counts <- rpois(length(lambda_bg), lambda_bg)
      bg_obs_counts <- rbinom(length(bg_true_counts), size = bg_true_counts, prob = total_pdet)
      bg_dropout_prob <- ifelse(bg_obs_counts <= densityObservable, 0.5, dropout_prob_densityD)
      bg_mask <- rbinom(length(bg_obs_counts), size = 1, prob = 1 - bg_dropout_prob)
      bg_obs_counts <- bg_obs_counts * bg_mask
      
      total_obs <- obs_flowering + bg_obs_counts
      
      
      tibble::tibble(
        year = y,
        doy = days,
        S_y = S_y[y],
        mu_y = mu_y[y],
        sigma_y = sigma_y[y],
        N_y = N_y[y],
        abundance = A_d_c,
        true_flowering = true_flowering,
        flowering_obs = obs_flowering,
        obs_background_count = bg_obs_counts,
        true_background_count = bg_true_counts,
        total_observed_count = total_obs,
        prop_observed_count_flowering = ifelse(total_obs > 0, obs_flowering / total_obs, 0),
        x = cell$x,
        y = cell$y,
        cell_id = cell$cell_id,
        pdet_day = pdet_day_component,
        pdet_year = pdet_year_component,
        pdet_spatial = pdet_spatial_component,
        total_detection_probability = total_pdet,
        detection_group = ifelse(equal_site_detection, "equal", cell$detection_group),
        prop_abundance = prop_c
        
      )
    })
    
    totals <- daily_results %>%
      dplyr::group_by(year, doy) %>%
      dplyr::summarize(
        true_flowering_total = sum(true_flowering),
        abundance_total = sum(abundance),
        .groups = "drop"
      )
    
    daily_results %>%
      dplyr::left_join(totals, by = c("year", "doy"))
  })
  
  # Optionally plot detection trend
  if (plot_detection_trend) {
    library(ggplot2)
    det_df <- tibble::tibble(year = 1:n_years, detection_rate = pdet_trend)
    print(
      ggplot(det_df, aes(x = year, y = detection_rate)) +
        geom_line(color = "steelblue", size = 1.2) +
        geom_point(color = "black") +
        theme_minimal() +
        labs(
          title = "Detection Probability Trend Over Years",
          y = "Detection Probability",
          x = "Year"
        )
    )
  }
  
  # Optionally plot spatial hotspots
  if (plot_spatial_hotspots) {
    library(ggplot2)
    cell_summary <- sim_obs %>%
      dplyr::group_by(cell_id, x, y) %>%
      dplyr::summarise(
        mean_total_observed = mean(total_observed_count, na.rm = TRUE),
        .groups = "drop"
      )
    print(
      ggplot(cell_summary, aes(x = x, y = y, fill = mean_total_observed)) +
        geom_tile() +
        scale_fill_viridis_c(option = "C") +
        theme_minimal() +
        coord_equal() +
        labs(
          title = "Spatial Hotspots: Mean Observed Counts per Cell",
          fill = "Mean Count",
          x = "Grid X",
          y = "Grid Y"
        )
    )
  }
  
  return(sim_obs)
}




#### including density obs as a probablistic dropout for masking
sim.latent.with.background.spatial.weekend.dropout <- function(
    n_years = 30, 
    days = 1:365, 
    mu_peak = 100, 
    dur = 12, 
    y_var = 0, 
    alpha_0 = 9, 
    alpha_1 = 2.5, 
    pbloom = 0.2, 
    pdetwend = 0.8, 
    pdetwday = 0.5, 
    min_det = 0.2, 
    max_det = 0.8, 
    grid_n = 5, 
    sdet = 0.6, 
    dirichlet_conc = 1e6, 
    bAbd_scaler = 0.2, 
    small_constant = 3, 
    plot_detection_trend = FALSE, 
    plot_spatial_hotspots = FALSE, 
    equal_site_detection = TRUE,  
    high_sdet = 0.9, 
    low_sdet = 0.3,  
    prop_high_sdet = 0.5,
    densityBasedMask = FALSE,
    beta0 = -1.5, 
    beta1 = 1.2, 
    beta2 = 1.0,
    plot_dropout_surface = FALSE
) {
  library(dplyr)
  library(tibble)
  library(purrr)
  
  # Temporal hyperparameters
  mu_mu <- rnorm(1, mean = mu_peak, sd = 10)
  sigma_mu <- abs(rnorm(1, mean = y_var, sd = 5))
  mu_sigma <- rnorm(1, mean = dur, sd = 3)
  sigma_sigma <- abs(rnorm(1, mean = 0, sd = 3))
  
  # Superbloom probability
  calibrate_bloom_beta <- function(mean, fixed_shape1 = 2) {
    shape1 <- fixed_shape1
    shape2 <- shape1 * (1 - mean) / mean
    list(pshape1 = shape1, pshape2 = shape2)
  }
  bloom_params <- calibrate_bloom_beta(pbloom)
  p_bloom <- rbeta(1, bloom_params$pshape1, bloom_params$pshape2)
  
  mu_y <- rnorm(n_years, mu_mu, sigma_mu)
  sigma_y <- abs(rnorm(n_years, mu_sigma, sigma_sigma))
  S_y <- rbinom(n_years, 1, p_bloom)
  log_N_y <- alpha_0 + alpha_1 * S_y
  N_y <- exp(log_N_y)
  
  calibrate_beta <- function(mean, fixed_shape1 = 10) {
    shape1 <- fixed_shape1
    shape2 <- shape1 * (1 - mean) / mean
    list(shape1 = shape1, shape2 = shape2)
  }
  
  spatial_grid <- expand.grid(x = 1:grid_n, y = 1:grid_n) %>%
    mutate(cell_id = paste0("cell_", x, "_", y))
  n_cells <- nrow(spatial_grid)
  
  if (equal_site_detection) {
    beta_params <- calibrate_beta(sdet, fixed_shape1 = 2)
    spatial_grid <- spatial_grid %>%
      mutate(
        pdet_spatial = rbeta(n(), beta_params$shape1, beta_params$shape2)
      )
  } else {
    n_high <- round(prop_high_sdet * n_cells)
    n_low <- n_cells - n_high
    high_cells <- sample(1:n_cells, n_high)
    spatial_grid <- spatial_grid %>%
      mutate(
        detection_group = ifelse(row_number() %in% high_cells, "high", "low"),
        pdet_spatial = case_when(
          detection_group == "high" ~ rbeta(1, calibrate_beta(high_sdet)$shape1, calibrate_beta(high_sdet)$shape2),
          detection_group == "low" ~ rbeta(1, calibrate_beta(low_sdet)$shape1, calibrate_beta(low_sdet)$shape2)
        )
      )
  }
  
  is_weekend <- function(doy) {
    weekday <- (doy - 1) %% 7 + 1
    weekday %in% c(6, 7)
  }
  weekend_flag <- is_weekend(days)
  pdet_day_vec <- ifelse(weekend_flag, pdetwend, pdetwday)
  pdet_trend <- seq(min_det, max_det, length.out = n_years)
  
  sim_obs <- purrr::map_dfr(1:n_years, function(y) {
    doy_density <- dnorm(days, mu_y[y], sigma_y[y])
    prop_cells <- MCMCpack::rdirichlet(1, rep(dirichlet_conc, n_cells))[1,]
    
    daily_results <- purrr::map_dfr(1:n_cells, function(i) {
      cell <- spatial_grid[i, ]
      prop_c <- prop_cells[i]
      A_d_c <- N_y[y] * prop_c * doy_density
      true_flowering <- rpois(length(A_d_c), A_d_c)
      
      pdet_day_component <- pdet_day_vec
      pdet_spatial_component <- cell$pdet_spatial
      pdet_year_component <- pdet_trend[y]
      total_pdet <- pdet_day_component * pdet_spatial_component * pdet_year_component
      
      obs_flowering <- rbinom(length(A_d_c), size = true_flowering, prob = total_pdet)
      
      if(densityBasedMask == TRUE){
      eps <- 1e-6
      p_mask <- plogis(beta0 + beta1 * log(obs_flowering + 1) + beta2 * log(total_pdet + eps))
      mask <- rbinom(length(p_mask), size = 1, prob = p_mask)
      obs_flowering <- obs_flowering * mask
      }
      
      lambda_bg <- small_constant + bAbd_scaler * A_d_c
      bg_true_counts <- rpois(length(lambda_bg), lambda_bg)
      bg_obs_counts <- rbinom(length(bg_true_counts), size = bg_true_counts, prob = total_pdet)
      if(densityBasedMask == TRUE){
      p_mask_bg <- plogis(beta0 + beta1 * log(bg_obs_counts + 1) + beta2 * log(total_pdet + eps))
      mask_bg <- rbinom(length(p_mask_bg), size = 1, prob = p_mask_bg)
      bg_obs_counts <- bg_obs_counts * mask_bg
      }
      total_obs <- obs_flowering + bg_obs_counts
      
      tibble(
        year = y,
        doy = days,
        S_y = S_y[y],
        mu_y = mu_y[y],
        sigma_y = sigma_y[y],
        N_y = N_y[y],
        abundance = A_d_c,
        true_flowering = true_flowering,
        flowering_obs = obs_flowering,
        obs_background_count = bg_obs_counts,
        true_background_count = bg_true_counts,
        total_observed_count = total_obs,
        prop_observed_count_flowering = ifelse(total_obs > 0, obs_flowering / total_obs, 0),
        x = cell$x,
        y = cell$y,
        cell_id = cell$cell_id,
        pdet_day = pdet_day_component,
        pdet_year = pdet_year_component,
        pdet_spatial = pdet_spatial_component,
        total_detection_probability = total_pdet,
        detection_group = ifelse(equal_site_detection, "equal", cell$detection_group),
        prop_abundance = prop_c
      )
    })
    
    totals <- daily_results %>%
      group_by(year, doy) %>%
      summarize(
        true_flowering_total = sum(true_flowering),
        abundance_total = sum(abundance),
        .groups = "drop"
      )
    
    daily_results %>%
      left_join(totals, by = c("year", "doy"))
  })
  
  if (plot_detection_trend) {
    library(ggplot2)
    det_df <- tibble(year = 1:n_years, detection_rate = pdet_trend)
    print(
      ggplot(det_df, aes(x = year, y = detection_rate)) +
        geom_line(color = "steelblue", size = 1.2) +
        geom_point(color = "black") +
        theme_minimal() +
        labs(
          title = "Detection Probability Trend Over Years",
          y = "Detection Probability",
          x = "Year"
        )
    )
  }
  
  if (plot_spatial_hotspots) {
    library(ggplot2)
    cell_summary <- sim_obs %>%
      group_by(cell_id, x, y) %>%
      summarize(
        sum_total_observed = sum(total_observed_count, na.rm = TRUE),
        .groups = "drop"
      )
    print(
      ggplot(cell_summary, aes(x = x, y = y, fill = sum_total_observed)) +
        geom_tile() +
        scale_fill_viridis_c(option = "C") +
        theme_minimal() +
        coord_equal() +
        labs(
          title = "Spatial Hotspots: Sum Observed Counts per Cell",
          fill = "Sum Count",
          x = "Grid X",
          y = "Grid Y"
        )
    )
  }
  
  if (plot_dropout_surface) {
    library(ggplot2)
    obs_range <- 1:50
    pdet_range <- seq(0.01, 1, by = 0.05)
    grid <- expand.grid(obs = obs_range, pdet = pdet_range)
    grid$p_mask <- plogis(beta0 + beta1 * log(grid$obs + 1) + beta2 * log(grid$pdet + 1e-6))
    print(
      ggplot(grid, aes(x = obs, y = p_mask, color = as.factor(round(pdet, 2)))) +
        geom_line() +
        labs(
          title = "Dropout Probability vs Observed Counts and Detection Probability",
          x = "Observed Counts",
          y = "Probability of Retaining Observation",
          color = "Detection Prob."
        ) +
        theme_minimal()
    )
  }
  
  return(sim_obs)
}
