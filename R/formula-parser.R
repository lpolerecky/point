#' Linear model parser
#'
#' @param data ion count data.
#' @param arg1 quosure of outcome variable.
#' @param arg2 quosure of predictor variable.
#' @param flag quosure of nominal predictor variable for outliers.
#' @param transformation transformation of the predictor ("log" or "ppt"
#'  dividing through thousand).
#' @param nest quosure of nest variable for mixed effect moddeling this is the
#'  variable for the nesting structure.
#' @param type type of model fitting (ratio method (RM), ordinary least squares
#'  (OLS), generalized least squares (GLS), special case for Quasi Simultaneous
#'  Arrival estimate (QSA) or linear mixed effects (LME).
#' @param execute either evaluate or return expression.
#'
#' @return call or evaluated call in context of data.
formula_parser <- function(data, arg1, arg2, flag = NULL, transformation = NULL,
                           nest = NULL, type = "OLS", execute = TRUE) {

  # fixed terms
  withRestarts(
    err <- tryCatch({
      fixed <- generate_fixed(!!arg1, !!arg2, type, flag = !!flag,
                              transformation)
      },
      error = function(e) {
      invokeRestart("rerun")
    }),
    rerun = function(...) fixed <<- generate_fixed(!!arg1, !!arg2, type,
                                                   flag = NULL, transformation)
  )

  # weight terms
  weights <- generate_weight(!!arg2, type, transformation)

  # random terms
  try(rand <- generate_random(!!arg2, type, !!nest, transformation), silent = TRUE)

  # Switch between types of linear model
  fml <- switch(
    type,
    OLS = rlang::call2("lm", fixed, data = rlang::expr(data)),
    QSA = rlang::call2("lme", fixed, random = rand, data = rlang::expr(data),
                       .ns = "nlme"),
    GLS = rlang::call2("gls", fixed, data = rlang::expr(data),
                       weights = weights, .ns = "nlme"),
    Rm =  rlang::call2("lm", fixed, data = rlang::expr(data), weights = weights),
    LME = rlang::call2("lme", fixed, random = rand, data = rlang::expr(data),
                       weights = weights, .ns = "nlme")
  )

  if(isTRUE(execute)) rlang::eval_tidy(fml, data) else fml
}

# The random terms
#
# arg2: predictor variable
# type: sort of model (ratio method (RM), linear mixed effect (LME),
#  generalized least square (GLS), and ordinary least square (OLS))
# transformation: predictor transform (parts per thousand (ppt) or log)
generate_random <- function(arg2, type, nest, transformation = NULL) {

  # symbols
  arg2 <- rlang::ensym(arg2)
  nest <- rlang::ensym(nest)

  # transform
  arg2 <- predictor_transformer(!! arg2, transformation)

  if (type == "QSA") {
    rhs <- paste(rlang::as_label(arg2), rlang::as_name(nest), sep = "|")
  } else if (type == "LME") {
    rhs <- paste(
      paste0("-1 +", rlang::as_label(arg2)),
      paste("execution", rlang::as_name(nest), sep = "/"),
      sep = "|"
    )
  }

  # formula
  rlang::new_formula(NULL, rlang::parse_expr(rhs))
}

# The fixed terms
#
# arg1: outcome variable
# arg2: predictor variable
# type: sort of model (ratio method (RM), linear mixed effect (LME),
#  generalized least square (GLS), and ordinary least square (OLS))
# transformation: predictor transform (parts per thousand (ppt) or log)
generate_fixed <- function(arg1, arg2, type, flag = NULL,
                           transformation = NULL) {

  # symbols
  arg1 <- rlang::ensym(arg1)
  arg2 <- rlang::ensym(arg2)
  try(flag <- rlang::ensym(flag), silent = TRUE)

  # most types can be generalized to one specific type, expect OLS
  if (type %in% c("Rm", "LME", "GLS")) {
    type <- "Rm"
  } else if (type %in% c("OLS", "QSA")) {
    type <- "OLS"
  }

  # transform
  arg2 <- predictor_transformer(!! arg2, transformation)

  # rhs of the formula
  if (type == "Rm") {
    rhs <- paste0("-1 +", rlang::as_label(arg2))
  } else if (type == "OLS") {
    rhs <- rlang::as_label(arg2)
  }

  # add flag variable if needed
  if (!is.null(flag)) {
    rhs <- paste0(rhs, "*", rlang::as_name(flag))
  }

  rlang::new_formula(arg1, rlang::parse_expr(rhs))
}

# The weighing term
#
# nlme weighing https://www.r-bloggers.com/2012/12/a-quick-note-in-weighting-with-nlme/
# use formula argument (gls(y ~ x, data = dat, weights = ~1/n))
#
# arg2: predictor variable
# type: sort of model (ratio method (RM), linear mixed effect (LME),
#  generalized least square (GLS), and ordinary least square (OLS))
# transformation: predictor transform (parts per thousand (ppt) or log)
generate_weight <- function(arg2, type, transformation = NULL) {

  # symbols
  arg2 <- rlang::ensym(arg2)

  # weighing structure for LME and GLS is equal but differs from `lm`
  if (type %in% c("LME", "GLS")) {
    type <- "LME"
  } else {
    type <- "Rm"
  }

  # transform
  arg2 <- predictor_transformer(!! arg2, transformation)

  # construct weight
  switch(
    type,
    Rm =  rlang::expr(1 / !! arg2),
    LME = rlang::new_formula(NULL, rlang::expr(1 / !! arg2))
  )
}

# predictor transformation
#
# function returns and expression
predictor_transformer <- function(arg2, transformation = NULL) {

  # symbols
  arg2 <- rlang::ensym(arg2)

  # transformations
  if (is.null(transformation)) {
    arg2
    # log transformation
  } else if (transformation == "log") {
    rlang::parse_expr(paste0("log(", as_name(arg2), ")"))
    # parts per thousand transformation
  } else if (transformation == "ppt") {
    rlang::parse_expr(paste0("I(c(", as_name(arg2), ")/ 1000)"))
  } else {
    stop("Transformation unkown!", call. = FALSE)
  }
}
