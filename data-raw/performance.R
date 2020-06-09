## code to prepare `performance` dataset

#-------------------------------------------------------------------------------
# Cooks D
#-------------------------------------------------------------------------------

# Run diagnostics for conversion rates estimates over multiple iterations
max <- 8
# empty list for iteration storage
ls.tb <- rlang::rep_named(vctrs::vec_cast(1:max, character()), lst())
# data
data <- sim_IC_systematic
# method
method <- "CooksD"
# variables for results
results_vars <- set_names(colnames(data))
# grouping vars
group_vars <- quos(simulation, trend, repetition)

results_vars <- parse_exprs(results_vars[!results_vars  %in%
                                         sapply(group_vars, as_name)]
                            )


ls.tb[[1]] <- lst(df = sim_IC_systematic, results = NULL)

# function wrapper for repeated treatment
rerun_dia_R <- function(out, input, variables, gr_by, method) {

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
  purrr::accumulate(rerun_dia_R,
                    variables = results_vars,
                    gr_by = group_vars,
                    method = method
                    )

# Reduce the results to a single dataframe
CD_all <- ls.tb %>%
  purrr::transpose() %>%
  purrr::pluck("results") %>%
  purrr::reduce(bind_rows)

usethis::use_data(CD_all, overwrite = TRUE, compress = "xz")

# Reduce the results to point statistics
CD_stat <- ls.tb %>%
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

usethis::use_data(CD_stat, overwrite = TRUE, compress = "xz")

# Evaluation of stepwise differences in performance
CD_eval <- CD_stat %>%
  eval_diag(RS_R_Xt.sim, n_R_Xt.sim, execution, simulation, trend, repetition)

usethis::use_data(CD_eval, overwrite = TRUE, compress = "xz")

#-------------------------------------------------------------------------------
# Cameca
#-------------------------------------------------------------------------------

# empty list for iteration storage
ls.tb <- rlang::rep_named(vctrs::vec_cast(1:max, character()), lst())

    # method
    method <- "Cameca"

    # grouping vars
    group_vars <- quos(simulation, trend, repetition, bl)

    ls.tb[[1]] <- lst(df = sim_IC_systematic, results = NULL)

    file.remove("data/sim_IC_systematic.rda")

    # execute repeated cycles of augmentation
    ls.tb <- ls.tb  %>%
      purrr::accumulate(rerun_dia_R,
                        variables = results_vars,
                        gr_by = group_vars,
                        method = method
                        )

    # reduce the results to a single dataframe
    Cameca_R <- ls.tb %>%
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

    usethis::use_data(Cameca_R, overwrite = TRUE, compress = "xz")

    # stepwise differences in perforance
    Cameca_stepwise <- Cameca_R %>%
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

