#' Linear model
#'
#' @param data ion count data
#' @param arg1 outcome variable
#' @param arg2 predictor variable
#' @param flag nominal predictor variable for outliers
#' @param trans transformation of the predictor ("log" or "ppt" dividing through
#'  thousand)
#' @param vorce for mixed effect modelling this dictates whether inter-analysis
#' isotope variance, or intra-analysis isotope variance is modelled
#' @param nest for mixed effect moddeling this is the variable for the nesting
#  structure
#' @param type type of model fitting (ratio method (RM), ordinary least squares
#' (OLS), generalized least squares (GLS), or linear mixed effects (LME)
#'
#' @return
lm_form <- function(data, arg1, arg2, flag = NULL, trans = NULL, vorce = NULL,
                    nest = NULL, type = "OLS") {

  lm_switch <- function(type, itrc = NULL, trns = NULL){

    if (type %in% c("Rm", "LME", "GLS")) type <- "Rm" else type = "OLS"
    if (is_quosure(itrc)) {
      type <- paste(type, deparse(substitute(itrc)), sep = "_")
      } else {
        if (is.null(trns)) {
          type <- paste(type, itrc, sep = "_")
        } else {
          type <- paste(type, trns, sep = "_")
        }
    }

    switch(
      type,
      Rm_flag =
        new_formula(
          quo_get_expr(arg1),
          parse_expr(paste0("-1 +", as_name(arg2), "*", as_name(flag)))
          ),
      Rm_ppt =
        new_formula(
          quo_get_expr(arg1),
          parse_expr(paste0("-1 + I(", as_name(arg2), "/ 1000)"))
          ),
      Rm_ =
        new_formula(
          quo_get_expr(arg1),
          parse_expr(paste0("-1 +", as_name(arg2)))
          ),
      Rm_log =
        new_formula(
          parse_expr(paste0("log(", as_name(arg1), ")")),
          parse_expr(paste0("-1 + log(", as_name(arg2), ")"))
          ),
      Rm_intra =
        new_formula(
          quo_get_expr(arg1),
          parse_expr(
            paste0(
              "-1 + I(",
              as_name(arg2),
              "/ 1000) + I(",
              as_name(arg2),
              "/ 1000) : execution")
            )
          ),
      OLS_ = new_formula(quo_get_expr(arg1), quo_get_expr(arg2))
        )

  }

  ran_switch <- function(clst, trns  = NULL){

    if (is_character(trns)) clst <- paste(clst, trns, sep = "_")
    switch(
      clst,
      QSA = new_formula(NULL, parse_expr(paste(as_name(arg2), as_name(nest), sep = "|"))),
      intra =
        new_formula(
          NULL,
          parse_expr(paste("-1 + execution", as_name(nest), sep = "|"))
          ),
      inter_ppt =
        new_formula(
          NULL,
          parse_expr(
            paste(
              paste0("-1 + I(c(", as_name(arg2), ")/ 1000)"),
              paste("execution", as_name(nest), sep = "/"),
              sep = "|"
              )
            )
          ),
      inter_log =
        new_formula(
          NULL,
          parse_expr(
            paste(
              paste0("-1 + log(", as_name(arg2), ")"),
              paste("execution", as_name(nest), sep = "/"),
              sep = "|")
            )
          )
    )

  }



  # Perform model in correct env
  data_env <- env(data = data)
  # Switch between types of linear model
  lm_method <- function(type) {
    switch(
      type,
      OLS = eval(call2("lm", lm_switch("OLS") , data = expr(data)), data_env),
      QSA = eval(call2("lme", lm_switch("OLS"), random = ran_switch("QSA"), data = expr(data),
                       .ns = "nlme"), data_env),
      GLS = eval(
        call2(
          "gls",
          lm_switch("GLS", trans),
          data = expr(data),
          weights = wght_switch("GLS", trans),
          .ns = "nlme"
        ),
        data_env
      ),
      Rm = eval(
        call2(
          "lm",
          lm_switch("Rm", flag),
          data = expr(data),
          weights = wght_switch("Rm")
          ),
        data_env
        ),
      LME = eval(
        call2(
          "lme",
          lm_switch("LME", vorce, trans),
          random = ran_switch(vorce, trans),
          data = expr(data),
          method = if_else(vorce == "inter", "REML", "ML"),
          weights = wght_switch("LME", trans),
          .ns = "nlme"
          ),
        data_env
        )
      )
  }

  lm_method(type)
}


# The weighing term
#
# nlme weighing https://www.r-bloggers.com/2012/12/a-quick-note-in-weighting-with-nlme/
# use formula argument (gls(y ~ x, data = dat, weights = ~1/n))
generate_weight <- function(arg2, type, transformation = NULL) {

  arg2 <- rlang::ensym(arg2)

  # weighing structure for LME and GLS is equal but differs from `lm`
  if (type %in% c("LME", "GLS")) {
    type <- "LME"
  } else {
    type <- "Rm"
  }

  # parts per thousand transformation
  if (is.null(transformation)) {
    tr <- arg2
  } else if (transformation == "log") {
    tr <- rlang::parse_expr(paste0("log(", as_name(arg2), ")"))
  } else if (transformation == "ppt") {
    tr <- rlang::parse_expr(paste0("I(c(", as_name(arg2), ")/ 1000)"))
  } else {
    stop("Transformation unkown!", .call = FALSE)
  }

  # construct weight
  switch(
    type,
    Rm =  rlang::expr(1 / !! tr),
    LME = rlang::new_formula(NULL, rlang::expr(1 / !! tr))
  )

}
