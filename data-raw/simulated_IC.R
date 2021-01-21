## code to prepare `simulated_IC` dataset
library(point)
library(rlang)
library(tidyverse)

# repetition
reps <- 15

# Types of R variation
var_R <- c("symmetric", "asymmetric")

# Varying linear trends in the ionization efficiency (per mille heavy ion)
var_T <- seq(0, 500, length.out = 6)

# Varying isotope offset (delta per mille)
var_I <- seq(0, -30, length.out = 6)

# Seeds for number generation
n_tot <- prod(sapply(list(var_T, var_R, var_I), length))
var_seed <- 1:n_tot

# Cross all possible parameter combinations
simu_R <- tidyr::crossing(.sys = var_T, .type = var_R, .devR = var_I) %>%
  mutate(.seed = var_seed) %>%
  transmute(params = purrr::pmap(lst(.sys, .type, .seed, .devR), lst))

# Function to map over params
map_R <- function(fn, args) {
  exec(fn, !!! args, .reps = reps, .ion1 = "13C", .ion2 = "12C")
  }

# Execute
simu_IC <- purrr::map2_dfr(
  rep("simu_R", n_tot),
  simu_R$params,
  map_R
  )

# Save data for application paper
saveRDS(simu_IC, "/home/amandus/Documents/work/manuscripts/Schobbenetal_SIMS_carb_stds/data/sim_IC.rds")

# Limited data internally only to show functionality of eval_diag()
simu_IC <- filter(
  simu_IC,
  trend.nm == 500,
  type.nm == "asymmetric",
  force.nm %in% c(-30, -12, 0),
  spot.nm %in% 1:5
  )
usethis::use_data(simu_IC, overwrite = TRUE, compress = "xz")

