## code to prepare `DATASET` dataset goes here

names_diag <- tibble(
  name = c(
    "Cameca",
    "Rm",
    "norm_E",
    "CooksD",
    "QQ",
    "CV",
    "QSA",
    "timeseries"
    ),
  label = c(
    "Cameca (R sigma rejection)",
    "Linear R model (residual sigma rejection)",
    "Residual vs. Leverage",
    "Linear R model (Cook's D)",
    "Normal QQ plot",
    "Scale-location plot",
    "QSA test",
    "timeseries"
    ),
  xaxis = c("ionct", "ionct", "hat_Xi", "ionct", "TQ", "hat_Y", "ionct", "time"),
  yaxis = c("ionct", "ionct", "studE", "ionct", "RQ", "studE", "R", "Xct")
  )

usethis::use_data(names_diag, overwrite = TRUE)



