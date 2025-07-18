### Proxy detection Probability Based on Empirical Data 

library(data.table)
library(tidyverse)
library(ggplot2)

# load empirical data
flowering_data <- fread("./data/processed/soCal-annual-data.csv")
# no background data, as we didnt have that for simulation?
#background_data <- fread("./data/processed/sCal-inat-background.csv")# build in pseudo-absences (days where detections didnt occur)

# build a generic background number of doys
all_doys <- 1:150
all_years <- unique(flowering_data$year)
all_hex50 <- unique(flowering_data$hex50_id)
# summarize as count type data
flowering_summary <- flowering_data %>% 
  group_by(year, doy, hex50_id) %>% 
  summarize(count = n_distinct(uuid), 
            .groups = 'drop') 

flowering_summary2 <- flowering_summary %>% 
  tidyr::complete(year = all_years, 
                  doy = all_doys, 
                  hex50_id = all_hex50,
                  fill = list(count = 0))

# Use data's P/A structure to figure out pdetwday/pdetwend
flowering_summary2$detected <- ifelse(flowering_summary2$count > 0, 1, 0)
# Weekday/weekend detection
is_weekend <- function(doy) {
  weekday <- (doy - 1) %% 7 + 1
  weekday %in% c(6, 7)
}
flowering_summary2$weekend_flag <- is_weekend(flowering_summary2$doy)
flowering_summary2 <- flowering_summary2 %>% 
  mutate(hex50_id = as.factor(hex50_id))

# use a mixed effects model 
library(lme4)
mod_glmm <- glmer(detected ~ is_weekend  + (1 | hex50_id),
                  data = flowering_summary2,
                  family = binomial)

# To extract site-specific detection probabilities for weekdays:
intercept <- fixef(mod_glmm)["(Intercept)"]
site_re  <- ranef(mod_glmm)$hex50_id[,1]
# All site log-odds for weekday
site_logodds <- intercept + site_re
site_probs <- plogis(site_logodds)
high_sdet <- max(site_probs)
low_sdet  <- min(site_probs)

# For weekends:
is_weekend_coef <- fixef(mod_glmm)["is_weekendTRUE"] # or ["is_weekend1"] depending on coding
site_logodds_weekend <- intercept + is_weekend_coef + site_re
site_probs_weekend <- plogis(site_logodds_weekend)
high_sdet_weekend <- max(site_probs_weekend)
low_sdet_weekend  <- min(site_probs_weekend)

# for years 
intercept <- fixef(mod_glmm)["(Intercept)"]
year_re   <- ranef(mod_glmm)$year[,1]
logodds_year <- intercept + year_re
det_probs_year <- plogis(logodds_year)
min_det <- min(det_probs_year)
max_det <- max(det_probs_year)