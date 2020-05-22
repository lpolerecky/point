#' Simulate ion count data
#'
#' @export
sim_R <- function(n = 3000,
                  N_range = 10 ^ 6,
                  reps = 1,
                  ion1,
                  ion2,
                  sys,
                  type,
                  baseR = NULL,
                  offsetR = NULL,
                  seed,
                  ...){

  average_n <- N_range / n
  start_n <- n


  tibble::tibble(simulation = type,
                 n = n,
                 N = as.integer(N_range),
                 R.input = R_gen(start_n,
                                 baseR,
                                 offsetR,
                                 input = "delta",
                                 type = type
                 ),
                 drift = seq(average_n * (1 - sys), average_n * (1 + sys),
                             length.out = start_n
                 ),
                 intercept = average_n
  ) %>%
    tidyr::expand_grid(., rep = c(1:reps), species = c(ion1, ion2)) %>%
    mutate(simulation = paste(.data$simulation, .data$rep, sep = "-"),
           seed = seed * rep
    ) %>%
# Convert common isotope N
    mutate(N = if_else(species == ion2, iso_conv(.data$N,
                                                 .data$R.input), .data$N
                       )
           ) %>%
# Calculate N of abundant isotope species
    group_by(.data$simulation, .data$species) %>%
    tidyr::nest() %>%
# Random variation (Number generation)
    mutate(N.sim= purrr::map(.data$data, ~N_gen(.x, N, n, seed))) %>%
    tidyr::unnest(cols = c(.data$data, .data$N.sim)) %>%
# Systematic variation
    mutate(diff = .data$drift - .data$intercept,
           diff = if_else(species == ion2,
                          as.double(iso_conv(.data$diff,
                                             .data$R.input
                          )
                          ),
                          .data$diff
           ),
           N.sim= .data$N.sim + .data$diff,
           Xt.sim = .data$N.sim,
           trend = paste0("linear trend (var: ", sys, ")")
    ) %>%
    ungroup() %>%
    select(-c(drift, intercept, diff))

}

#-------------------------------------------------------------------------------
# Random Poisson ion count generator
#-------------------------------------------------------------------------------
N_gen <- function(df, N, n, seed) {


  N <- enquo(N)
  n <- enquo(n)
  seed <- enquo(seed)

  N <- df %>% pull(!! N)
  n <- df %>% pull(!! n)
  seed <- df %>% pull(!! seed)

  set.seed(seed)

  Nsim <- as.double(rpois(n = n, lambda = N / n))

}

#-------------------------------------------------------------------------------
# Calculate common isotope count from rare isotope
#-------------------------------------------------------------------------------
iso_conv <- function(N, R.sim)  as.integer(N * (1 / R.sim))


#-------------------------------------------------------------------------------
# Create isotopic gradients and offsets
#-------------------------------------------------------------------------------
R_gen <- function(reps, baseR, offsetR, input = "delta", type) {

  baseR <- calib_R(baseR,
                   standard = "VPDB",
                   type = "composition",
                   input = input,
                   output = "R"
  )

  offsetR <- calib_R(offsetR,
                     standard = "VPDB",
                     type = "composition",
                     input = input,
                     output = "R"
  )

  if (type == "ideal") {
    R.sim <- rep(baseR, reps)
    return(R.sim)
  }

  if (type == "constant"){

    R.sim <- approx(c(1, 5 * reps /6, reps),
                    c(baseR, offsetR, offsetR),
                    n = reps ,
                    method = "constant"
    )$y
    return(R.sim)
  }

  if (type == "gradient") {


    R.sim <- approx(c(1, reps),
                    c(offsetR, baseR),
                    n = reps ,
                    method = "linear"
    )$y

    return(R.sim)
  }
}
