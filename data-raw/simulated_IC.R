## code to prepare `simulated_IC` dataset

# repetition
reps <-50

# Types of R variation
var_R <- c("ideal", "constant", "gradient")

# Varying linear trends in the ionization efficieny
var_T <- seq(0, 0.5, length.out = 6)

# Seeds for number generation
tot_length <- length(var_T) * length(var_R)
var_seed <- 1:tot_length

sim_R.ext <- tidyr::crossing(sys = var_T, type = var_R) %>%
  mutate(seed = var_seed) %>%
  transmute(params = purrr::pmap(lst(sys, type, seed), lst))


sim_IC <- purrr::map2_dfr(rep("sim_R", tot_length) %>%
                               set_names(nm = paste("run", 1:tot_length)),
                             sim_R.ext$params,
                             function(fn, args) exec(fn,
                                                     !!! args,
                                                     reps = reps,
                                                     ion1 = "13C",
                                                     ion2 = "12C",
                                                     baseR = 5,
                                                     offsetR = -55
                             )
                          ) %>%
  tidyr::separate(simulation, sep = "-", c("simulation", "repetition"))

usethis::use_data(sim_IC, overwrite = TRUE, compress = "xz")

# limited dataset only to show range in systematic of an isotopically
# ideal substance but different systematic trends
sim_IC_systematic <- sim_IC %>%
  filter(trend == "linear trend (var: 0)" |
         trend == "linear trend (var: 0.5)"
         )

usethis::use_data(sim_IC_systematic, overwrite = TRUE, compress = "xz")

# limited dataset only to show range in systematic and random R trends for
sim_IC_extremes <- sim_IC  %>%
  filter(trend == "linear trend (var: 0)" |
         trend == "linear trend (var: 0.5)",
         repetition == 1
         )

usethis::use_data(sim_IC_extremes, overwrite = TRUE)

# limited dataset only to show range in systematic of an isotopically
# ideal substance
sim_IC_ideal <- sim_IC_extremes %>%
  filter(simulation == "ideal")

usethis::use_data(sim_IC_ideal, overwrite = TRUE)

