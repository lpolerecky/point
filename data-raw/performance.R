## code to prepare `performance` dataset

#-------------------------------------------------------------------------------
# variables simulated dataset
#-------------------------------------------------------------------------------

expr_R_stat <- expr_R(
  Xt = "Xt.sim",
  N = "N.sim",
  species = "species",
  ion1 = "13C",
  ion2 = "12C"
  )

grps <- quos(simulation, trend, repetition, iso_offset)

#-------------------------------------------------------------------------------
# Cooks D
#-------------------------------------------------------------------------------


ls.tb <- diag_R(
  sim_IC,
  method = "CooksD",
  args = expr_R_stat,
  reps = 1,
  !!! grps,
  output = "complete",
  plot = FALSE
  )

# Evaluation of performance
CD_eval <- eval_diag(
  ls.tb,
  expr_R_stat,
  flag,
  !!! grps
  )

usethis::use_data(CD_eval, overwrite = TRUE, compress = "xz")

#-------------------------------------------------------------------------------
# Cameca
#-------------------------------------------------------------------------------


ls.tb <- diag_R(
  sim_IC,
  method = "Cameca",
  args = expr_R_stat,
  reps = 1,
  !!! grps,
  output = "complete",
  plot = FALSE
  )

# Evaluation of performance
CM_eval <- eval_diag(
  ls.tb,
  expr_R_stat,
  flag,
  !!! grps
  )

usethis::use_data(CM_eval, overwrite = TRUE, compress = "xz")
