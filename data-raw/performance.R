  ## code to prepare `performance` dataset

  #-------------------------------------------------------------------------------
  # Cooks D
  #-------------------------------------------------------------------------------

  expr_R_stat <- expr_R(Xt = "Xt.sim",
                        N = "N.sim",
                        species = "species",
                        ion1 = "13C",
                        ion2 = "12C"
                        )

  grps <- quos(simulation, trend, repetition)

  varXt <- sim_IC_systematic %>%
    mutate(varXt = readr::parse_number(trend)) %>%
    select(!!!grps, varXt)

  ls.tb <- diag_R(sim_IC_systematic,
                  method = "CooksD",
                  args = expr_R_stat,
                  reps = 6,
                  !!! grps
                  )

  # Reduce the results to point statistics
  CD_stat <- reduce_diag(ls.tb, type = "df", expr_R_stat, !!! grps)

  usethis::use_data(CD_stat, overwrite = TRUE, compress = "xz")

  # Evaluation of stepwise differences in performance
  CD_eval <- CD_stat %>%
    group_by(!!! grps) %>%
    arrange(desc(execution)) %>%
    mutate(execution = execution,
           nz = n_R_Xt.sim,
           diff_n = c(diff(n_R_Xt.sim, 1), NA),
           diff_RS =  c(diff(RS_R_Xt.sim, 1), NA),
           SSD_z =  diff_RS ^ 2 * diff_n,
           SSO_z =  RS_R_Xt.sim ^ 2 * n_R_Xt.sim,
           n2_z = SSD_z / SSO_z,
           chi_n2 = n2_z * nz
           ) %>%
    ungroup()

  usethis::use_data(CD_eval, overwrite = TRUE, compress = "xz")

  #-------------------------------------------------------------------------------
  # Cameca
  #-------------------------------------------------------------------------------

  ls.tb <- diag_R(sim_IC_systematic,
                  method = "Cameca",
                  args = expr_R_stat,
                  reps = 6,
                  !!! grps
                  )

  # reduce the results to a single dataframe
  CM_stat <- reduce_diag(ls.tb, type = "df", expr_R_stat, !!! grps)

  usethis::use_data(CM_stat, overwrite = TRUE, compress = "xz")

  # Evaluation of stepwise differences in performance
  CM_eval <- CM_stat %>%
    group_by(!!! grps) %>%
    arrange(desc(execution)) %>%
    mutate(execution = execution,
           nz = n_R_Xt.sim,
           diff_n = c(diff(n_R_Xt.sim, 1), NA),
           diff_RS =  c(diff(RS_R_Xt.sim, 1), NA),
           SSD_z =  diff_RS ^ 2 * diff_n,
           SSO_z =  RS_R_Xt.sim ^ 2 * n_R_Xt.sim,
           n2_z = SSD_z / SSO_z,
           chi_n2 = n2_z * nz
    ) %>%
    ungroup()

  usethis::use_data(CM_eval, overwrite = TRUE, compress = "xz")

