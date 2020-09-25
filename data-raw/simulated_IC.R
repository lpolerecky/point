## code to prepare `simulated_IC` dataset

# repetition
reps <- 10

# Types of R variation
var_R <- c("ideal", "symmetric", "asymmetric")

# Varying linear trends in the ionization efficiency
var_T <- seq(0, 500, length.out = 4)

# # varying isotope offset
var_I <- c(-5,-25,-45,-65)

# Seeds for number generation
tot_length <- length(var_T) * length(var_R) * length(var_I)
var_seed <- 1:tot_length

# Cross all possible parameter combinations
sim_R.ext <- tidyr::crossing(sys = var_T, type = var_R, offsetR = var_I) %>%
  mutate(seed = var_seed) %>%
  transmute(params = purrr::pmap(lst(sys, type, seed, offsetR), lst))


sim_IC <- purrr::map2_dfr(rep("sim_R", tot_length),
                          sim_R.ext$params,
                          function(fn, args) exec(fn,
                                                  !!! args,
                                                  reps = reps,
                                                  ion1 = "13C",
                                                  ion2 = "12C",
                                                  baseR = 5
                                                  )
                          )


usethis::use_data(sim_IC, overwrite = TRUE, compress = "xz")

# limited dataset only to show range in systematic of an isotopically
# ideal substance but different systematic trends
sim_IC_systematic <- sim_IC %>%
  filter(trend == "linear trend (var: 0)" |
         trend == "linear trend (var: 500)"
         )

usethis::use_data(sim_IC_systematic, overwrite = TRUE, compress = "xz")

# limited dataset only to show range in systematic and random R trends for
sim_IC_extremes <- sim_IC  %>%
  filter(trend == "linear trend (var: 0)" |
         trend == "linear trend (var: 500)",
         repetition == 1
         )

usethis::use_data(sim_IC_extremes, overwrite = TRUE)

# limited dataset only to show range in systematic of an isotopically
# ideal substance
sim_IC_ideal <- sim_IC_extremes %>%
  filter(simulation == "ideal")

usethis::use_data(sim_IC_ideal, overwrite = TRUE)

