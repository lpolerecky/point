test_that("weighing structure can be build",{
  # Ratio method
  expect_equal(
    generate_weight(Xt.pr.12C, "Rm"), parse_expr("1/Xt.pr.12C")
    )
  expect_equal(
    generate_weight(Xt.pr.12C, "Rm", "ppt"), parse_expr("1/I(c(Xt.pr.12C)/1000)")
    )
  expect_equal(
    generate_weight(Xt.pr.12C, "Rm", "log"), parse_expr("1/log(Xt.pr.12C)")
    )
  # LME or GLS
  fml <- generate_weight(Xt.pr.12C, "GLS")
  expect_equal(
    fml, as.formula("~1/Xt.pr.12C", env = rlang::get_env(fml))
  )
  fml <- generate_weight(Xt.pr.12C, "GLS", "ppt")
  expect_equal(
    fml, as.formula("~1/I(c(Xt.pr.12C)/1000)", env = rlang::get_env(fml))
  )
  fml <- generate_weight(Xt.pr.12C, "GLS", "log")
  expect_equal(
    fml, as.formula("~1/log(Xt.pr.12C)", env = rlang::get_env(fml))
  )
})
