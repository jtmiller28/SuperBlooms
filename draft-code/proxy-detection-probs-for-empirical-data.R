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

# Use data's P/A structure 
flowering_summary2$detected <- ifelse(flowering_summary2$count > 0, 1, 0)


# Weekday/weekend detection
is_weekend <- function(doy) {
  weekday <- (doy - 1) %% 7 + 1
  weekday %in% c(6, 7)
}
flowering_summary2$weekend_flag <- is_weekend(flowering_summary2$doy)
flowering_summary2 <- flowering_summary2 %>% 
  mutate(hex50_id = as.factor(hex50_id)) %>% 
  mutate(year = as.factor(year))

# Use a mixed effects model to model detection based on whether the doy is a weekend vs weekday, including
library(lme4)
mod_glmm <- glmer(detected ~ is_weekend  + (1 | hex50_id) + (1 | year),
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
mean_sdet <- mean(site_probs)
print(mean_sdet)
# For weekends:
is_weekend_coef <- fixef(mod_glmm)["is_weekend"] # or ["is_weekend1"] depending on coding
site_logodds_weekend <- intercept + is_weekend_coef + site_re
site_probs_weekend <- plogis(site_logodds_weekend)
site_probs_weekday <- plogis(intercept + site_re)
high_sdet_weekend <- max(site_probs_weekend)
low_sdet_weekend  <- min(site_probs_weekend)
high_sdet_weekday <- max(site_probs_weekday)
low_sdet_weekday <- min(site_probs_weekday)
mean_sdet_weekend <- mean(site_probs_weekend)
mean_sdet_weekday <- mean(site_probs_weekday)

# for years 
intercept <- fixef(mod_glmm)["(Intercept)"]
year_re   <- ranef(mod_glmm)$year[,1]
probs_by_year_wkdays <- plogis(intercept + year_re) # for weekdays
probs_by_year_wkends <- plogis(intercept + is_weekend_coef + year_re) # for weekends
min_det_wkday <- min(probs_by_year_wkdays)
max_det_wkday <- max(probs_by_year_wkdays)
min_det_wkend <- min(probs_by_year_wkends)
max_det_wkend <- max(probs_by_year_wkends)
mean_det_wkday <- mean(probs_by_year_wkdays)
mean_det_wkend <- mean(probs_by_year_wkends)
print(mean_det_wkday)

# Model Coefficients 
# For each Site and Year
# Weekday logit: logit_p = intercept + site_effect + year_effect
# Weekend logit: logit_p = intercept + is_weekend_coef + site_effect + year_effect
intercept <- fixef(mod_glmm)["(Intercept)"] # extract intercept
is_weekend_coef <- fixef(mod_glmm)["is_weekend"] # extract weekend coef
site_re <- ranef(mod_glmm)$hex50_id # extract site re as a dataframe
year_re <- ranef(mod_glmm)$year # extract year re  as a dataframe

# set aside names
site_names <- rownames(site_re)
year_names <- rownames(year_re)

# grab effect
site_re_vec <- site_re[,1]
year_re_vec <- year_re[,1]

# expand combos 
probs_all <- expand.grid(site = site_names, year = year_names, is_weekend = c(FALSE, TRUE)) %>% 
  mutate(
    logit = fixef(mod_glmm)["(Intercept)"] + 
      ifelse(is_weekend, fixef(mod_glmm)["is_weekend"], 0) + 
      site_re_vec[site] + 
      year_re_vec[year], 
    prob = plogis(logit)
  )

low_sdet  <- min(probs_all$prob)
high_sdet <- max(probs_all$prob)
mean_sdet <- mean(probs_all$prob)

ggplot(data = probs_all, mapping = aes(x = site, y = prob)) + 
  geom_col() + 
  facet_wrap(~year)
