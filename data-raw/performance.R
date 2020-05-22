## code to prepare `performance` dataset

# Run diagnostics for conversion rates estimates over multiple iterations
ls.tb <- lst(1, 2, 3, 4, 5, 6, 7, 8)

ls.tb[[1]] <- lst(df = sim_IC_systematic, results = NULL)

# function wrapper for repeated treatment
rerun_dia_R <- function(out, input,
                        variables = set_names(colnames(ls.tb[[1]]$df)),
                        gr_by = quos(simulation, trend, repetition),
                        method = "CooksD"
) {

  variables <- variables[!variables %in% sapply(gr_by, as_name)]

  out <- diag_R(out$df,
                method = method,
                args = expr_R(Xt = "Xt.sim",
                              N = "N.sim",
                              species = "species",
                              ion1 = "13C",
                              ion2 = "12C"
                ),
                !!! gr_by,
                output = "flag"
  )

  # save augmented dataframe for next cycle
  df.aug <- out %>%
    filter(flag == "non-influential") %>%
    select(!!!gr_by, !!!variables)

  # save results with ID for executo
  results <- select(out, -c(!!!variables)) %>%
    mutate(execution = input) %>%
    distinct(!!! gr_by, .keep_all = TRUE)

  out <- lst(df = df.aug, results = results)

  return(out)

}

# Execute repeated cycles of augmentation
ls.tb <- ls.tb  %>%
  accumulate(rerun_dia_R)

# Reduce the results to a single dataframe
CooksD_results <- ls.tb %>%
  transpose() %>%
  pluck("results") %>%
  reduce(bind_rows)

usethis::use_data(CooksD_results, overwrite = TRUE, compress = "xz")

# Reduce the results to a single dataframe
tb.data <- ls.tb %>%
  purrr::transpose() %>%
  purrr::pluck("df") %>%
  purrr::map_dfr(~stat_R(.x,
                         Xt.sim,
                         N.sim,
                         species,
                         ion1 = "13C",
                         ion2 = "12C",
                         simulation,
                         trend,
                         repetition
  ),
  .id = "execution"
  )

# Stepwise differences in performance
CooksD_stepwise <- tb.data %>%
  group_by(simulation, trend, repetition) %>%
  arrange(desc(execution)) %>%
  # SSE
  transmute(execution = execution,
            diff_n = c(diff(n_R_Xt.sim , 1), NA),
            var_R  = c(diff(RS_R_Xt.sim, 1)^2, NA),
            SSA  = var_R * diff_n,
            # SST
            SST = RS_R_Xt.sim^2 * n_R_Xt.sim,
            # R2 of adjustement (variance explained)
            w2 = SSA / SST,
            # Chi squared based on R mean model (X2 with df = S-1 = 2-1)
            Chi_w2 = w2 *  n_R_Xt.sim
  )

usethis::use_data(CooksD_stepwise, overwrite = TRUE, compress = "xz")



ls.tb <- lst(1,2,3,4,5,6,7,8)

ls.tb[[1]] <- lst(df = sim_IC_systematic, results = NULL)

# execute repeated cycles of augmentation
ls.tb   <- ls.tb  %>% accumulate(rerun_dia_R,  method = "Cameca")

# reduce the results to a single dataframe
tb.data <- ls.tb %>%
  purrr::transpose() %>%
  purrr::pluck("df") %>%
  purrr::map_dfr(~stat_R(.x,
                         Xt.sim,
                         N.sim,
                         species,
                         ion1 = "13C",
                         ion2 = "12C",
                         simulation,
                         trend,
                         repetition
  ),
  .id = "execution"
  )

# stepwise differences in perforance
Cameca_stepwise <- tb.data %>%
  group_by(simulation, trend, repetition) %>%
  arrange(desc(execution)) %>%
  # SSE
  transmute(execution = execution,
            diff_n = c(diff(n_R_Xt.sim , 1), NA),
            var_R  = c(diff(RS_R_Xt.sim, 1)^2, NA),
            SSA  = var_R * diff_n,
            # SST
            SST = RS_R_Xt.sim^2 * n_R_Xt.sim,
            # R2 of adjustement (variance explained)
            w2 = SSA / SST,
            # Chi squared based on R mean model (X2 with df = S-1 = 2-1)
            Chi_w2 = w2 *  n_R_Xt.sim
  )

usethis::use_data(Cameca_stepwise, overwrite = TRUE, compress = "xz")

