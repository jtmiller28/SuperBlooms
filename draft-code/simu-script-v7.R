### V7: Removing weekned vs week after finding little effect down 


library(data.table)

# create a function that bootstraps known background data to generate simulated backgrounds 
gen.background <- function(sample.year.level = c("High", "Medium", "Low")){
  # load real sample background data
  background_data <- fread("/blue/guralnick/millerjared/SuperBlooms/data/processed/test50km-background-site-data.csv")
  # choose which sampling level of data to use
  if(sample.year.level == "High") {
    background_data <- background_data %>% filter(year == 2021)
  } else if(sample.year.level == "Medium") {
    background_data <- background_data %>% filter(year == 2019)
  } else if(sample.year.level == "Low") {
    background_data <- background_data %>% filter(year == 2018)
  }
  n_days <- 150 # set static, consider making flexible later
  # summarize total counts per doy 
  background_data <- background_data %>% 
    group_by(doy) %>% 
    summarize(n = n())
  # sample & bootstrap
  sim_obs <- sample(background_data$n, size = n_days, replace = TRUE)
  df_boot <- data.frame(doy = 1:n_days, obs = sim_obs)
  # return this bootstrapped simulated distribution
  return(df_boot)
}

sim.superbloom <- function(
    n_years = 10,
    year_type_sequence = NULL,
    days = 1:150,
    mu_peak = 100,
    dur = 12,
    y_var = 0,
    sigma_N = 0.5,
    alpha_0 = 9,
    alpha_1 = 2.5,
    yearly_abd_varation = FALSE,
    pbloom = 0.1,
    pdetwend = 0.8,
    pdetwday = 0.5,
    min_det = 0.2,
    max_det = 0.8,
    grid_n = 5,
    dirichlet_conc = 1e6,
    site_heterogeniety = FALSE,
    plot_detection_trend = FALSE,
    plot_spatial_hotspots = FALSE,
    equal_site_detection = TRUE,
    sdet = 0.6,
    high_sdet = 0.9,
    low_sdet = 0.3,
    prop_high_sdet = 0.5,
    randomizeSeed = TRUE, 
    seedSet = 202
) {
  requireNamespace("dplyr")
  requireNamespace("data.table")
  requireNamespace("purrr")
  requireNamespace("tibble")
  requireNamespace("MCMCpack")
  if(is.null(year_type_sequence)) {
    stop("Needs year_type_sequence (a vector of 'Low', 'Medium', 'High')")
  }
  if(length(year_type_sequence) != n_years) {
    stop("year_type_sequence must be the same length as n_years")
  }
  if(randomizeSeed == TRUE){
    set.seed(sample(1:1e6, 1)) # pick a randomized seed
  }
  # Temporal hyperparameters
  mu_mu <- rnorm(1, mean = mu_peak, sd = 10)
  sigma_mu <- abs(rnorm(1, mean = y_var, sd = 1))
  mu_sigma <- rnorm(1, mean = dur, sd = 3)
  sigma_sigma <- abs(rnorm(1, mean = 0, sd = 1))
  # Setup spatial grid and superbloom probabilities
  spatial_grid <- expand.grid(x = 1:grid_n, y = 1:grid_n) %>%
    dplyr::mutate(cell_id = paste0("cell_", x, "_", y))
  n_cells <- nrow(spatial_grid)
  # Assign "high superbloom probability" status to a portion of cells
  n_high_prob_cells <- round(0.6 * n_cells)
  high_prob_cells <- sample(1:n_cells, n_high_prob_cells)
  cell_superbloom_probs <- rep(0.15, n_cells) # low-prob for most cells
  cell_superbloom_probs[high_prob_cells] <- 0.5 # high-prob for special cells
  spatial_grid$superbloom_prob <- cell_superbloom_probs
  # Assign detection probabilities for sites
  calibrate_beta <- function(mean, fixed_shape1 = 10) {
    shape1 <- fixed_shape1
    shape2 <- shape1 * (1 - mean) / mean
    list(shape1 = shape1, shape2 = shape2)
  }
  if (equal_site_detection) {  # All same detection
    beta_params <- calibrate_beta(sdet, fixed_shape1 = 2)
    spatial_grid <- spatial_grid %>%
      dplyr::mutate(pdet_spatial = rbeta(n(), beta_params$shape1, beta_params$shape2))
    spatial_grid$detection_group <- "equal"
  } else {  # High/low detection cells
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
  
  # Weekday/weekend detection (REMOVE)
  # is_weekend <- function(doy) {
  #   weekday <- (doy - 1) %% 7 + 1
  #   weekday %in% c(6, 7)
  # }
  # weekend_flag <- is_weekend(days)
  # pdet_day_vec <- ifelse(weekend_flag, pdetwend, pdetwday)
  # Backgrounds
  gen.background <- function(sample.year.level = c("High", "Medium", "Low")){
    background_data <- data.table::fread("/blue/guralnick/millerjared/SuperBlooms/data/processed/test50km-background-site-data.csv")
    if(sample.year.level == "High") {
      background_data <- background_data %>% dplyr::filter(year == 2021)
    } else if(sample.year.level == "Medium") {
      background_data <- background_data %>% dplyr::filter(year == 2019)
    } else if(sample.year.level == "Low") {
      background_data <- background_data %>% dplyr::filter(year == 2018)
    }
    n_days <- 150
    background_data <- background_data %>%
      dplyr::group_by(doy) %>%
      dplyr::summarize(n = dplyr::n())
    sim_obs <- sample(background_data$n, size = n_days, replace = TRUE)
    df_boot <- data.frame(doy = 1:n_days, obs = sim_obs)
    return(df_boot)
  }
  year_backgrounds <- purrr::map(year_type_sequence, gen.background)
  bg_medians <- sapply(year_backgrounds, function(df) median(df$obs))
  bg_scaled <- (bg_medians - min(bg_medians)) / (max(bg_medians) - min(bg_medians) + 1e-8)
  pdet_trend <- min_det + bg_scaled * (max_det - min_det)
  # Year-specific flowering phenology shape
  mu_y <- rnorm(n_years, mu_mu, sigma_mu)
  sigma_y <- abs(rnorm(n_years, mu_sigma, sigma_sigma))
  sigma_y[sigma_y < 1] <- 1  # Lower bound for sd
  
  # Superbloom probability (allows just imputing a direct probability into the function for superbloom prob)
  calibrate_bloom_beta <- function(mean, fixed_shape1 = 2) {
    shape1 <- fixed_shape1
    shape2 <- shape1 * (1 - mean) / mean
    list(pshape1 = shape1, pshape2 = shape2)
  }
  bloom_params <- calibrate_bloom_beta(pbloom)
  p_bloom <- rbeta(1, bloom_params$pshape1, bloom_params$pshape2)
  
  # Draw superbloom years
  S_y <- rbinom(n_years, 1, pbloom)
  
  sim_obs <- purrr::map_dfr(1:n_years, function(y) {
    bg_this_year <- year_backgrounds[[y]]
    doy_density <- dnorm(days, mu_y[y], sigma_y[y])
    if (site_heterogeniety) {
      prop_cells <- MCMCpack::rdirichlet(1, rep(dirichlet_conc, n_cells))[1,]
    } else {
      prop_cells <- rep(1 / n_cells, n_cells)
    }
    site_results <- purrr::map_dfr(1:n_cells, function(i) {
      cell <- spatial_grid[i, ]
      prop_c <- prop_cells[i]
      # Superbloom logic: only assign S_y_site if S_y[y]==1
      if (S_y[y] == 1) {
        S_y_site <- rbinom(1, 1, cell$superbloom_prob)
      } else {
        S_y_site <- 0
      }
      if(yearly_abd_varation){
        log_N_y_cell <- alpha_0 + alpha_1 * S_y_site + rnorm(1, 0, sigma_N)
      } else {
        log_N_y_cell <- alpha_0 + alpha_1 * S_y_site
      }
      N_y_cell <- exp(log_N_y_cell)
      A_d_c <- N_y_cell * prop_c * doy_density
      true_flowering <- rpois(length(A_d_c), A_d_c)
      # pdet_day_component <- pdet_day_vec
      # pdet_spatial_component <- cell$pdet_spatial
      # pdet_year_component <- pdet_trend[y]
      # total_pdet <- pdet_day_component * pdet_spatial_component * pdet_year_component
      # obs_flowering <- rbinom(length(A_d_c), size = true_flowering, prob = total_pdet)
      epsilon <- 1e-6 # min val 
      clamp <- function(x) pmin(pmax(x, epsilon), 1 - epsilon)
      
      # use clamp to ensure probabilities are within bounds
      #pdet_day_component <- clamp(pdet_day_vec)
      pdet_spatial_component <- clamp(cell$pdet_spatial)
      pdet_year_component <- clamp(pdet_trend[y])
      
      # convert to logits 
      #logit_day <- qlogis(pdet_day_component)
      logit_spat <- qlogis(pdet_spatial_component)
      logit_year <- qlogis(pdet_year_component)
      
      # combine additively, vector by day
      total_logit <- logit_spat + logit_year
      total_pdet <- plogis(total_logit)
      
      # obs process 
      obs_flowering <- rbinom(length(A_d_c), size = true_flowering, prob = total_pdet)
      
      bg_obs_counts <- bg_this_year$obs
      total_obs <- obs_flowering + bg_obs_counts
      tibble::tibble(
        year = y,
        doy = days,
        S_y = S_y[y],
        S_y_site = S_y_site,
        mu_y = mu_y[y],
        sigma_y = sigma_y[y],
        N_y_site = N_y_cell,
        abundance = A_d_c,
        true_background_count = bg_obs_counts,
        true_flowering = true_flowering,
        flowering_obs = obs_flowering,
        total_obs_count = total_obs,
        prop_observed_count_flowering = ifelse(total_obs > 0, obs_flowering / total_obs, 0),
        x = cell$x,
        y = cell$y,
        cell_id = cell$cell_id,
        superbloom_prob = cell$superbloom_prob,
        pdet_year = pdet_year_component,
        pdet_spatial = pdet_spatial_component,
        total_detection_probability = total_pdet,
        detection_group = cell$detection_group,
        prop_abundance = prop_c
      )
    })
    totals <- site_results %>%
      dplyr::group_by(year, doy) %>%
      dplyr::summarize(
        true_flowering_total = sum(true_flowering),
        abundance_total = sum(abundance),
        .groups = "drop"
      )
    site_results %>%
      dplyr::left_join(totals, by = c("year", "doy"))
  })
  debug_years <- sim_obs %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      mu_y = mean(mu_y),
      sigma_y = mean(sigma_y),
      mean_N_y = mean(N_y_site),
      prop_superbloom_year = first(S_y),
      prop_sites_with_superbloom = mean(S_y_site),
      .groups="drop"
    )
  print(debug_years)
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
  if (plot_spatial_hotspots) {
    library(ggplot2)
    cell_summary <- sim_obs %>%
      dplyr::group_by(cell_id, x, y) %>%
      dplyr::summarise(
        cell_obs_count = sum(total_obs_count, na.rm = TRUE),
        .groups = "drop"
      )
    print(
      ggplot(cell_summary, aes(x = x, y = y, fill = cell_obs_count)) +
        geom_tile() +
        scale_fill_viridis_c(option = "C") +
        theme_minimal() +
        coord_equal() +
        labs(
          title = "Spatial Hotspots: Sum Observed Counts per Cell",
          fill = "Cell Sum Count",
          x = "Grid X",
          y = "Grid Y"
        )
    )
  }
  return(sim_obs)
}