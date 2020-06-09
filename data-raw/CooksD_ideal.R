## code to prepare `CooksD_ideal` dataset

sim <- sim_IC

CooksD_ideal <- filter(sim, simulation == "ideal" & rep < 4) %>%
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
        ) %>%
  group_by(trend) %>%
  summarise(Chi_R2 = mean(Chi_R2),
            SE_beta = mean(SE_beta ),
            S_Xt.sim.12C = mean(S_Xt.sim.12C)
            )

usethis::use_data(CooksD_ideal, overwrite = TRUE, compress = "xz")
# remove large simulation dataset
use_directory("data")
paths <- fs::path("data", sim, ext = "rda")
file.remove(paths)
