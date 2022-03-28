test_that("whole formula can be generated", {
  # only formulas
  # OLS
  call <- formula_parser(real_IC, Xt.pr.13C, Xt.pr.12C, execute = FALSE)
  expect_equal(
    rlang::as_label(call),
    "lm(Xt.pr.13C ~ Xt.pr.12C, data = data)"
  )
  # GLS
  call <- formula_parser(real_IC, Xt.pr.13C, Xt.pr.12C, type = "GLS", execute = FALSE)
  expect_equal(
    rlang::as_label(call),
    "nlme::gls(Xt.pr.13C ~ -1 + Xt.pr.12C, data = data, weights = ~1/Xt.pr.12C)"
  )
  # LME
  call <- formula_parser(real_IC, Xt.pr.13C, Xt.pr.12C, type = "LME", nest = nest, execute = FALSE)
  expect_equal(
    deparse(call),
   c("nlme::lme(Xt.pr.13C ~ -1 + Xt.pr.12C, random = ~-1 + Xt.pr.12C | ",
    "    execution/nest, data = data, weights = ~1/Xt.pr.12C)")
  )


  # evaluate with data

})


test_that("predictors can be transformed", {
  expect_equal(
    predictor_transformer(Xt.pr.12C, "ppt"), parse_expr("I(c(Xt.pr.12C)/1000)")
  )
  expect_equal(
    predictor_transformer(Xt.pr.12C, "log"), parse_expr("log(Xt.pr.12C)")
  )
})

test_that("fixed formula terms can be assembled", {
  # different model types
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "OLS")
  expect_equal(
    fml, as.formula("Xt.pr.13C ~ Xt.pr.12C", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "Rm")
  expect_equal(
    fml, as.formula("Xt.pr.13C ~ -1 + Xt.pr.12C", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "GLS")
  expect_equal(
    fml, as.formula("Xt.pr.13C ~ -1 + Xt.pr.12C", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "LME")
  expect_equal(
    fml, as.formula("Xt.pr.13C ~ -1 + Xt.pr.12C", env = rlang::get_env(fml))
  )
  # transformations
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "OLS", transformation = "ppt")
  expect_equal(
    fml, as.formula("Xt.pr.13C ~ I(c(Xt.pr.12C)/1000)", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "OLS", transformation = "log")
  expect_equal(
    fml, as.formula("Xt.pr.13C ~ log(Xt.pr.12C)", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "GLS", transformation = "ppt")
  expect_equal(
    fml, as.formula("Xt.pr.13C ~ -1 + I(c(Xt.pr.12C)/1000)", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "GLS", transformation = "log")
  expect_equal(
    fml, as.formula("Xt.pr.13C ~ -1 + log(Xt.pr.12C)", env = rlang::get_env(fml))
  )
  # flag alone
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "OLS", flag)
  expect_equal(
    fml, as.formula("Xt.pr.13C ~ Xt.pr.12C * flag", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "GLS", flag)
  expect_equal(
    fml, as.formula("Xt.pr.13C ~ -1 + Xt.pr.12C * flag", env = rlang::get_env(fml))
  )
  # flag and transformation
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "OLS", flag, transformation = "ppt")
  expect_equal(
    fml, as.formula("Xt.pr.13C ~ I(c(Xt.pr.12C)/1000) * flag", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "GLS", flag, transformation = "ppt")
  expect_equal(
    fml, as.formula("Xt.pr.13C ~ -1 + I(c(Xt.pr.12C)/1000) * flag", env = rlang::get_env(fml))
  )
})

test_that("random structure can be created", {
  # LME
  fml <- generate_random(Xt.pr.12C, "LME", nest)
  expect_equal(
    fml, as.formula("~-1 + Xt.pr.12C | execution/nest", env = rlang::get_env(fml))
  )
  # LME with log transform
  fml <- generate_random(Xt.pr.12C, "LME", nest , "log")
  expect_equal(
    fml, as.formula("~-1 + log(Xt.pr.12C) | execution/nest", env = rlang::get_env(fml))
  )
  # LME with ppt transform
  fml <- generate_random(Xt.pr.12C, "LME", nest , "ppt")
  expect_equal(
    fml, as.formula("~-1 + I(c(Xt.pr.12C)/1000) | execution/nest", env = rlang::get_env(fml))
  )
  # special case
  #
  # QSA
  fml <- generate_random(Xt.pr.12C, "QSA", nest)
  expect_equal(
    fml, as.formula("~Xt.pr.12C | nest", env = rlang::get_env(fml))
  )
})


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
