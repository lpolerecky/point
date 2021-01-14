## code to prepare nm_stat and R_stat dataset goes here

nm_stat_Xt <- tibble(
  nm = c(
    "n",
    "Ntot",
    "M",
    "S",
    "RS",
    "SeM",
    "hat_S",
    "hat_SeM",
    "chi2"
    ),
  latex = c(
    "$n$",
    "$N_{tot}$",
    "$\\bar{X}$",
    "$s_X$",
    "$\\epsilon_X \\,$ (\u2030)",
    "$s_\\bar{X}$",
    "$\\hat{s}_N$",
    "$\\hat{s}_\\bar{N}$",
    "$\\chi^{2}$"
  ),
  label = c(
    "observations",
    "total count",
    "mean",
    "standard deviation",
    "relative standard deviation",
    "standard error of the mean",
    "predicted standard deviation",
    "predicted standard error of the mean",
    "chi squared"
  )
)

usethis::use_data(nm_stat_Xt, overwrite = TRUE)

nm_stat_R <- tibble(
  nm = c(
    "n",
    "M",
    "S",
    "RS",
    "SeM",
    "RSeM",
    "hat_S",
    "hat_RS",
    "hat_SeM",
    "hat_RSeM",
    "chi2"
    ),
  latex = c(
    "$n$",
    "$\\bar{R}$",
    "$s_{R}$",
    "$\\epsilon_{R} \\,$ (\u2030)",
    "$s_{\\bar{R}}$",
    "$\\epsilon_{\\bar{R}} \\,$ (\u2030)",
    "$\\hat{s}_{R}$",
    "$\\hat{\\epsilon}_{R} \\,$ (\u2030)",
    "$\\hat{s}_{\\bar{R}}$",
    "$\\hat{\\epsilon}_{\\bar{R}} \\,$ (\u2030)",
    "$\\chi^{2}$"
    ),
  label = c(
    "observations",
    "mean",
    "standard deviation",
    "relative standard deviation",
    "standard error of the mean",
    "relative standard error of the mean",
    "predicted standard deviation",
    "predicted relative standard deviation",
    "predicted standard error of the mean",
    "predicted relative standard error of the mean",
    "chi squared"
    )
)

usethis::use_data(nm_stat_R, overwrite = TRUE)
