  ## code to prepare `CooksD_ideal` dataset

  CD_ideal <- filter(sim_IC, simulation == "ideal" & repetition < 4) %>%
    diag_R(df = .,
           method = "CooksD",
           args = expr_R(Xt = "Xt.sim",
                         N = "N.sim",
                         species = "species",
                         ion1 = "13C",
                         ion2 = "12C"
                         ),
          reps = 1,
          simulation,
          trend,
          repetition
          ) %>%
    purrr::transpose() %>%
    purrr::pluck("results") %>%
    purrr::pluck("2") %>%
    group_by(trend) %>%
    summarise(Chi_R2 = mean(Chi_R2),
              SE_beta = mean(SE_beta )
              # ,
              # S_Xt.sim.12C = mean(S_Xt.sim.12C)
              )

  usethis::use_data(CD_ideal, overwrite = TRUE, compress = "xz")

  # remove large simulation dataset
  file.remove("data/sim_IC.rda")
