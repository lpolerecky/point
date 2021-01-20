## code to prepare nm_stat and R_stat dataset goes here

names_stat_Xt <- tibble(
  name = c(
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

usethis::use_data(names_stat_Xt, overwrite = TRUE)

names_stat_R <- tibble(
  name = c(
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

usethis::use_data(names_stat_R, overwrite = TRUE)
