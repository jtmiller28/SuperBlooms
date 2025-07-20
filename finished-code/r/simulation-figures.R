### Simulation Figures
### Author: JT Miller

## Creating Figures For Simulation outputs for Presentation 

# load libraries
library(data.table)
library(ggplot2)
library(MCMCpack)
library(tidyverse)
library(lme4)
library(latex2exp)
source("/blue/guralnick/millerjared/SuperBlooms/draft-code/simu-script-v7.R") # simulation functions


### Create Figure Image Export 
year_types <- c(rep("Low", 3),rep("Medium", 3), rep("High", 4))
sim_data <- sim.superbloom(
  n_years = 10,
  year_type_sequence = year_types,
  yearly_abd_varation = TRUE, 
  days = 1:150,
  mu_peak = 100,
  dur = 12,
  y_var = 10,
  alpha_0 = 9,
  alpha_1 = 2,
  pbloom = 0.3,
  min_det = 1,
  max_det = 1,
  grid_n = 5,
  dirichlet_conc = 1,
  equal_site_detection = FALSE,
  high_sdet = 1,
  low_sdet = 1,
  prop_high_sdet = 0.5,
  plot_detection_trend = FALSE,
  plot_spatial_hotspots = FALSE,
  seedSet = 100
)

# find a superbloom Y and non-superbloom Y
superbloom_y <- filter(sim_data, S_y == 1) %>% pull(year) %>% unique()
superbloom_y <- max(superbloom_y)

nonsuperbloom_y <- filter(sim_data, S_y == 0) %>% pull(year) %>% unique()
nonsuperbloom_y <- max(nonsuperbloom_y)

out %>% filter(year %in% c(superbloom_y, nonsuperbloom_y)) %>% 
  ggplot() + 
  geom_col(aes = )

### Create Math Image Export
library(latex2exp)
png("/blue/guralnick/millerjared/SuperBlooms/figures/latex_step1_math.png", width = 400, height = 200)
plot(1:10, type = "n", axes = FALSE, xlab = "", ylab = "")
text(5, 5, TeX('
S_y \\sim Bernoulli(p) \\\\
\\log(N_y) = \\alpha_0 + \\alpha_1 \\times S_y \\\\
\\lambda_{y,d} = A_{y,d} = N_y \\times \\text{NormalPDF}(d|\\mu_y,\\sigma_y) \\\\
\\text{TrueCounts}_{y,d} \\sim \\text{Poisson}(\\lambda_{y,d}) \\\\
\\text{ObservedCounts}_{y,d} \\sim \\text{Binomial}(\\text{TrueCounts}_{y,d}, p_{det})
'), cex = 1.2)
dev.off()