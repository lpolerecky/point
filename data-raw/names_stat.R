## code to prepare nm_stat and R_stat dataset goes here

names_stat_X <- tibble(
  name = c(
    "n",
    "tot",
    "M",
    "S",
    "RS",
    "SeM",
    "hat_S",
    "hat_RS",
    "hat_SeM",
    "chi2"
    ),
  origin = c("t", "N", rep("X", 4), rep("N", 4)),
  label = c(
    "observations",
    "total count",
    "mean",
    "standard deviation",
    "relative standard deviation",
    "standard error of the mean",
    "predicted standard deviation",
    "predicted relative standard deviation",
    "predicted standard error of the mean",
    "chi squared"
  )
)

usethis::use_data(names_stat_X, overwrite = TRUE)

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
  ratio = "R",
  origin = c("t", rep("X", 5), rep("N", 5)),
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
