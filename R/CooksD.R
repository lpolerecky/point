#' @rdname Cameca
#'
#' @export
CV <- Rm <- norm_E <- CooksD <- QQ <- IR <- function(
  .df, .ion1, .ion2, ..., .Xt = Xt.pr, .t = t.nm,
  .output = "complete", .hyp = "none", .alpha_level = 0.05){

  # Grouping
  gr_by <- enquos(...)

  # function name
  fun_nm <- as_name(match.call()[[1]])

  # check if name has hypothesis test
  if(fun_nm != "CV" & .hyp == "bp") {
    stop("Wrong hypothesis test for this method.")
    }
  if(fun_nm != "QQ" & (.hyp == "norm" | .hyp == "ttest")) {
    stop("Wrong hypothesis test for this method.")
    }
  if(fun_nm != "IR" & .hyp == "ljung") {
    stop("Wrong hypothesis test for this method.")
    }
  if (!(fun_nm == "QQ" | fun_nm == "CV" | fun_nm == "IR") &  .hyp != "none") {
    .hyp <- "none"
    warning("No hypothesis test avalaible for this method.")
    }

  # Quoting the call (user-supplied expressions)
  Xt <- enquo(.Xt)
  t <- enquo(.t)

  # Heavy isotope
  Xt1 <- quo_updt(Xt, post = .ion1) # count rate
  # Light isotope
  Xt2 <- quo_updt(Xt, post = .ion2) # count rate

  # Execute
  df <- nest_R_lm(.df, gr_by, Xt1, Xt2, t, method = fun_nm, .hyp = .hyp, .alpha_level = .alpha_level)

  # Output
  if (fun_nm == "IR") return(unnest(select(df, -c(t, data)), cols = c(extr, flag)))
  if (.output == "flag") return(unnest(select(df, -data), cols = c(t, extr, flag)))
  if (.output == "complete") return(unnest(select(df, -t), cols = c(data, extr, flag)))
}

#-------------------------------------------------------------------------------
# Not exportet helper functions
#-------------------------------------------------------------------------------

# nest lm (args as quos)
nest_R_lm <- function(df, gr_by, Xt1, Xt2, t, method, .hyp, .alpha_level){

  tidyr::nest(df, t = !!t, data = -c(!!! gr_by)) %>%
  mutate(
    R_lm =
      purrr::map(data, ~lm_form(.x, Xt1, Xt2, type = "Rm")),
    aug =
      purrr::map(R_lm, broom::augment),
    extr =
      purrr::map(aug, ~transmute_reg(.x, Xt1, Xt2, method)),
    extr =
      purrr::map(extr, ~QQ_trans(.x, method, .hyp, .alpha_level)),
    extr =
      purrr::map(extr, ~IR_trans(.x, method, .hyp, .alpha_level)),
    extr =
      purrr::map2(aug, extr, ~bp_wrap(.x, .y, Xt2, method, .hyp, .alpha_level)),
    flag =
      purrr::map2(extr, aug, ~flag_set(.x, .y, Xt1, Xt2, method, .alpha_level))
    ) %>%
  select(-c(R_lm, aug))

}

# prefix for modelled values
prefix <- c("hat", "hat_s", "hat_e")

# augment function transform  and rename variables to standards of point
transmute_reg <- function(df, Xt1, Xt2, type) {

  # predicted variable and sigma level boundaries
  hat_args <- purrr::map(prefix, ~quo_updt(Xt1, pre = .x))
  hat_nm <- paste(prefix, as_name(Xt1), sep = "_")

  # predicted sigma level boundaries studentized residuals
  hat_args_studE_nm <- paste(prefix, "studE", sep = "_")
  hat_args_studE <- parse_exprs(hat_args_studE_nm)

  # model args
  args <- quos(
    hat_E = .data$.resid,
    !!hat_args[[1]] := .data$.fitted,
    !!hat_args[[2]] := sigma_calc(hat_E),
    studE = .data$.std.resid,
    hat_Xi = .data$.hat,
    CooksD = .data$.cooksd,
    )

  if (type == "Rm") args <- args[names(args) %in% c(hat_nm, "hat_E")]
  if (type == "norm_E") args <- args[!names(args) %in% hat_nm]
  if (type == "CooksD") args <- args[names(args) %in% c(hat_nm, "hat_E", "CooksD")]
  if (type == "CV") args <- args[names(args) %in% c(hat_nm[1], "studE")]
  if (type == "QQ"| type == "IR") args <- args["studE"]

  # Execute
  transmute(df, !!!args)
}

# create flag variable
flag_set <- function(df1, df2, Xt1, Xt2, type, .alpha_level){

  # residual sigma level boundaries
  hat_args_sigma <- purrr::map(prefix, ~quo_updt(Xt1, pre = .x))

  # CI boundaries of QQ
  hat_args_QQ <- paste(prefix, "RQ", sep = "_") %>%
    parse_exprs()

  # CI boundaries acf
  hat_args_acf <- paste(prefix, "acf", sep = "_") %>%
    parse_exprs()

  if (type == "Rm") {
   df <- flagger(
     df1,
     hat_E,
     !!hat_args_sigma[[2]],
     fct = qnorm((1 -.alpha_level / 2))
     )
   return(df)
   }
  if (type == "CooksD" | type == "norm_E") {
   df <- transmute(
     df1,
     flag = factor(if_else(CooksD < {4 / (n() - 2)}, "confluent", "divergent"))
     )
   return(df)
   }
  if (type == "QQ") {
    df <- flagger(
      df1,
      QE,
      !!hat_args_QQ[[3]],
      fct = qt((1 - .alpha_level / 2), n())
      )
    return(df)
   }
  if (type == "CV") return(flagger(df1, studE, 3.5))
  if (type == "IR") return(flagger(df1, acf, !!hat_args_acf[[3]]))

}

flagger <- function(df, value, bound, fct = 1){

    transmute(
      df,
      flag = factor(
        if_else(
          between({{value}}, - fct * {{bound}}, fct * {{bound}}),
          "confluent",
          "divergent"
          )
        )
      )

}

# quantile transformations and hypothesis tests
QQ_trans <- function(df, type, .hyp, .alpha_level) {

  if (type!= "QQ") return(df)

  # predicted variable and standard error and CI boundaries
  hat_args <- paste(prefix, "RQ", sep = "_") %>%
    parse_exprs()

  # Normality hypothesis test
  if (.hyp == "norm") {
    hyp_result <- nortest::ad.test(df$studE)$p.value
    Ha <- "Ha (non-normal)"
    H0 <- "H0 (normal)"
    }
  # t-test flag for mu0 (aka the conditional mean of residual) being zero
  if (.hyp == "ttest") {
    hyp_result <- t.test(df$studE, mu = 0)$p.value
    Ha <- "Ha (mu0 is not zero)"
    H0 <- "H0 (mu0 is zero)"
    }

  df <- transmute(
    df,
    RQ = unname(quantile(studE, probs = ppoints(n()))), #vector_probs(n()))),
    # Calculate normal (Theoretical) quantiles using mean and standard deviation
    TQ = qnorm(ppoints(n()), mean(RQ), sd(RQ)),
    QE = RQ - TQ,
    # The standard error
    !!hat_args[[1]] := mean(RQ) + sd(RQ) * TQ,
    !!hat_args[[3]] := hat_QR_se(RQ, TQ, ppoints(n()), n()),
    )

  if (.hyp != "none") {
    return(mutate(df, hyp = if_else(hyp_result < .alpha_level, Ha, H0)))
    } else {
      return(df)
      }
}

# auto-correlation and hypothesis tests
IR_trans <- function(df, type, .hyp, .alpha_level) {

  if (type!= "IR") return(df)

  # independence test
  if (.hyp == "ljung") {
    hyp_result <- stats::Box.test(df$studE, type = "Ljung-Box")$p.value
    Ha <- "Ha (dependence of residuals)"
    H0 <- "H0 (independence of residuals)"
    }

  # predicted variable and standard error and CI boundaries
  hat_args <- paste(prefix, "acf", sep = "_") %>%
    parse_exprs()

  acf <- acf(df$studE, plot = FALSE)
  si <- qnorm(.alpha_level / 2) / sqrt(length(df$studE))

  df <- tibble(
    lag = as.vector(acf$lag)[-1],
    acf := as.vector(acf$acf)[-1],
    !!hat_args[[3]] := -si
    )

  if (.hyp != "none") {
    return(mutate(df, hyp = if_else(hyp_result < .alpha_level, Ha, H0)))
    } else {
      return(df)
      }
  }

# Hetroscadasticity test (Breusch Pagan test)(level of confidence 95%;
# cut-off 0.05 for H0 rejection)
bp_wrap <- function(df1, df2, Xt2, type, .hyp, .alpha_level){

  # Breusch Pagan test
  if (type == "CV" & .hyp == "bp") {
    Chi_R2 <- custom_bp(df1, Xt2)
    Ha <- "Ha (heteroskedasticity)"
    H0 <- "H0 (homoskedasticity)"
    df2 <- mutate(
      df2,
      hyp = if_else(Chi_R2 > qchisq((1 - .alpha_level), df = 1), Ha, H0)
      )
    return(df2)
    } else {
      return(df2)
  }

  }

custom_bp <- function(df, Xt2){

  res_lm <- lm_form(df, quo(.std.resid), Xt2)
  R2 <- pull(broom::glance(res_lm), `r.squared`)
  SE_beta <- pull(broom::tidy(res_lm), std.error)[2]

  return(R2 * length(R2))
  }

sigma_calc <- function(res) sqrt(sum((res ^ 2)) / (length(res) - 1))

# # use the formula i - 0.5/ in, for i = 1,..,n
# # this is a vector of the n probabilities (theoretical cumulative distribution function CDF)
# vector_probs <- function(n){
#   ((1:unique(n)) - 0.5) / (unique(n))
# }

# standard error of quantiles model
hat_QR_se <- function(RQ, TQ, pb, n){
  (sd(RQ) / dnorm(TQ)) * sqrt((pb * (1 - pb))/ unique(n))
}
# confidence interval regression model
hat_Y_se <- function(sigma, hat_Xi){
  sigma * sqrt(hat_Xi)
}
