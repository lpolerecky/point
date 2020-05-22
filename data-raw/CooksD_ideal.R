## code to prepare `CooksD_ideal` dataset

CooksD_ideal <- filter(sim_IC, simulation == "ideal" & rep < 4) %>%
  diag_R(method = "CooksD",
         args = expr_R(Xt = "Xt.sim",
                       N = "N.sim",
                       species = "species",
                       ion1 = "13C",
                       ion2 = "12C"
                       ),
        simulation,
        trend,
        repetition,
        output = "complete"
        )%>%
  group_by(trend) %>%
  summarise(Chi_R2 = mean(Chi_R2),
            SE_beta = mean(SE_beta ),
            S_Xt.sim.12C = mean(S_Xt.sim.12C)
            )

usethis::use_data(CooksD_ideal, overwrite = TRUE)
