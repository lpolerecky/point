names_plot <- tibble::tibble(
  name = c(
    "Cameca",
    "Rm",
    "norm_E",
    "CooksD",
    "QQ",
    "CV",
    "QSA",
    "IR",
    "timeseries"
    ),
  inference = c(rep("eval_diag", 4), rep("external", 5)),
  label = c(
    "Cameca (R sigma rejection)",
    "Linear R model (residual sigma rejection)",
    "Residual vs. Leverage",
    "Linear R model (Cook's D)",
    "Normal QQ plot",
    "Scale-location plot",
    "QSA test",
    "ACF plot",
    "timeseries"
    ),
  xaxis = c("ionct", "ionct", "hat_Xi", "ionct", "TQ", "hat_Y", "ionct", "", "time"),
  yaxis = c("ionct", "ionct", "studE", "ionct", "RQ", "studE", "R", "", "Xct")
  )

usethis::use_data(names_plot, overwrite = TRUE)
