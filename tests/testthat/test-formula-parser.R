test_that("whole formula can be generated", {
  # only formulas
  quo <- rlang::quo
  # OLS
  call <- formula_parser(real_IC, quo(Xt.pr.13C), quo(Xt.pr.12C),
                         execute = FALSE)
  expect_equal(
    rlang::as_label(call),
    "lm(Xt.pr.13C ~ Xt.pr.12C, data = data)"
  )
  # GLS
  call <- formula_parser(real_IC, quo(Xt.pr.13C), quo(Xt.pr.12C), type = "GLS",
                         execute = FALSE)
  expect_equal(
    rlang::as_label(call),
    "nlme::gls(Xt.pr.13C ~ -1 + Xt.pr.12C, data = data, weights = ~1/Xt.pr.12C)"
  )
  # LME
  call <- formula_parser(real_IC, quo(Xt.pr.13C), quo(Xt.pr.12C), type = "LME",
                         nest = quo(nest), execute = FALSE)
  expect_equal(
    deparse(call),
   c("nlme::lme(Xt.pr.13C ~ -1 + Xt.pr.12C, random = ~-1 + Xt.pr.12C | ",
    "    execution/nest, data = data, weights = ~1/Xt.pr.12C)")
  )
  # LME with transformation
  call <- formula_parser(real_IC, quo(Xt.pr.13C), quo(Xt.pr.12C), type = "LME",
                         nest = quo(nest), transformation = "ppt",
                         execute = FALSE)
  expect_equal(
    deparse(call),
    c("nlme::lme(Xt.pr.13C ~ -1 + Xt.pr.12C, random = ~-1 + Xt.pr.12C | ",
      "    execution/nest, data = data, weights = ~1/Xt.pr.12C)")
  )
  # with flag
  call <- formula_parser(real_IC, quo(Xt.pr.13C), quo(Xt.pr.12C), quo(flag),
                         type = "Rm", execute = FALSE)
  expect_equal(
    rlang::as_label(call),
    "lm(Xt.pr.13C ~ -1 + Xt.pr.12C * flag, data = data, weights = 1/Xt.pr.12C)"
  )

  # evaluate with data
  xc <- dplyr::filter(real_IC, .data$file.nm == "2018-01-19-GLENDON_1_1")
  xc <- cov_R(xc, c("12C", "13C"), sample.nm, file.nm)
  # simple
  expect_snapshot(
    formula_parser(xc, quo(Xt.pr.13C), quo(Xt.pr.12C))
  )
  expect_snapshot(
    formula_parser(xc, quo(Xt.pr.13C), quo(Xt.pr.12C), type = "Rm")
  )
  # flag
  xc <- dplyr::mutate(xc, flag = dplyr::if_else(t.nm < 15, "divergent",
                                                "confluent"))
  expect_snapshot(
    formula_parser(xc, quo(Xt.pr.13C), quo(Xt.pr.12C), quo(flag), type = "Rm")
  )
  # GLS with transformation
  expect_snapshot(
    formula_parser(xc, quo(Xt.pr.13C), quo(Xt.pr.12C), type = "GLS",
                   transformation = "ppt")
  )
  # LME with transformation
  expect_snapshot(
    formula_parser(real_IC, quo(Xt.pr.13C), quo(Xt.pr.12C), type = "LME",
                   nest = fil.nm, transformation = "ppt")
  )
})


test_that("predictors can be transformed", {
  expect_equal(
    predictor_transformer(Xt.pr.12C, "ppt"),
    rlang::parse_expr("I(c(Xt.pr.12C)/1000)")
  )
  expect_equal(
    predictor_transformer(Xt.pr.12C, "log"),
    rlang::parse_expr("log(Xt.pr.12C)")
  )
})

test_that("fixed formula terms can be assembled", {
  # different model types
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "OLS")
  expect_equal(
    fml,
    as.formula("Xt.pr.13C ~ Xt.pr.12C", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "Rm")
  expect_equal(
    fml,
    as.formula("Xt.pr.13C ~ -1 + Xt.pr.12C", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "GLS")
  expect_equal(
    fml,
    as.formula("Xt.pr.13C ~ -1 + Xt.pr.12C", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "LME")
  expect_equal(
    fml,
    as.formula("Xt.pr.13C ~ -1 + Xt.pr.12C", env = rlang::get_env(fml))
  )
  # transformations
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "OLS", transformation = "ppt")
  expect_equal(
    fml,
    as.formula("Xt.pr.13C ~ I(c(Xt.pr.12C)/1000)", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "OLS", transformation = "log")
  expect_equal(
    fml,
    as.formula("Xt.pr.13C ~ log(Xt.pr.12C)", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "GLS", transformation = "ppt")
  expect_equal(
    fml,
    as.formula("Xt.pr.13C ~ -1 + I(c(Xt.pr.12C)/1000)", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "GLS", transformation = "log")
  expect_equal(
    fml,
    as.formula("Xt.pr.13C ~ -1 + log(Xt.pr.12C)", env = rlang::get_env(fml))
  )
  # flag alone
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "OLS", flag)
  expect_equal(
    fml,
    as.formula("Xt.pr.13C ~ Xt.pr.12C * flag", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "GLS", flag)
  expect_equal(
    fml,
    as.formula("Xt.pr.13C ~ -1 + Xt.pr.12C * flag", env = rlang::get_env(fml))
  )
  # flag and transformation
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "OLS", flag, transformation = "ppt")
  expect_equal(
    fml,
    as.formula("Xt.pr.13C ~ I(c(Xt.pr.12C)/1000) * flag", env = rlang::get_env(fml))
  )
  fml <- generate_fixed(Xt.pr.13C, Xt.pr.12C, "GLS", flag, transformation = "ppt")
  expect_equal(
    fml,
    as.formula("Xt.pr.13C ~ -1 + I(c(Xt.pr.12C)/1000) * flag", env = rlang::get_env(fml))
  )
})

test_that("random structure can be created", {
  # LME
  fml <- generate_random(Xt.pr.12C, "LME", nest)
  expect_equal(
    fml,
    as.formula("~-1 + Xt.pr.12C | execution/nest", env = rlang::get_env(fml))
  )
  # LME with log transform
  fml <- generate_random(Xt.pr.12C, "LME", nest , "log")
  expect_equal(
    fml,
    as.formula("~-1 + log(Xt.pr.12C) | execution/nest",
               env = rlang::get_env(fml))
  )
  # LME with ppt transform
  fml <- generate_random(Xt.pr.12C, "LME", nest , "ppt")
  expect_equal(
    fml,
    as.formula("~-1 + I(c(Xt.pr.12C)/1000) | execution/nest",
               env = rlang::get_env(fml))
  )
  # special case
  #
  # QSA
  fml <- generate_random(Xt.pr.12C, "QSA", nest)
  expect_equal(
    fml,
    as.formula("~Xt.pr.12C | nest", env = rlang::get_env(fml))
  )
})


test_that("weighing structure can be build",{
  # Ratio method
  expect_equal(
    generate_weight(Xt.pr.12C, "Rm"),
    rlang::parse_expr("1/Xt.pr.12C")
  )
  expect_equal(
    generate_weight(Xt.pr.12C, "Rm", "ppt"),
    rlang::parse_expr("1/I(c(Xt.pr.12C)/1000)")
  )
  expect_equal(
    generate_weight(Xt.pr.12C, "Rm", "log"),
    rlang::parse_expr("1/log(Xt.pr.12C)")
  )
  # LME or GLS
  fml <- generate_weight(Xt.pr.12C, "GLS")
  expect_equal(
    fml,
    as.formula("~1/Xt.pr.12C", env = rlang::get_env(fml))
  )
  fml <- generate_weight(Xt.pr.12C, "GLS", "ppt")
  expect_equal(
    fml,
    as.formula("~1/I(c(Xt.pr.12C)/1000)", env = rlang::get_env(fml))
  )
  fml <- generate_weight(Xt.pr.12C, "GLS", "log")
  expect_equal(
    fml,
    as.formula("~1/log(Xt.pr.12C)", env = rlang::get_env(fml))
  )
})
