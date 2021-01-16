#' @rdname Cameca
#'
#' @export
CV <- Rm <- norm_E <- CooksD <- QQ <- function(
  .df, .ion1, .ion2, ..., .Xt = Xt.pr, .t = t.nm, .output = "flag",
  .hyp = "none"
  ){

  # Grouping
  gr_by <- enquos(...)

  # function name
  fun_nm <- as_name(match.call()[[1]])

  # check if name has hypothesis test
  if(fun_nm == "QQ" & .hyp == "bp") {
    stop("Wrong hypothesis test for this method.")
    }
  if(fun_nm == "CV" & (.hyp == "norm" | .hyp == "ttest")) {
    stop("Wrong hypothesis test for this method.")
    }
  if (!(fun_nm == "QQ" | fun_nm == "CV") &  .hyp != "none") {
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
  df <- nest_R_lm(.df, gr_by, Xt1, Xt2, t, method = fun_nm, .hyp = .hyp)

  # Output
  if (.output == "flag") return(unnest(select(df, -data), cols = c(t, extr, flag)))
  if (.output == "complete") return(unnest(select(df, -t), cols = c(data, extr, flag)))
  # if (.output == "nest") return(df)
}


#' Independence of residuals
#'
#' @export
IR <- function(.df, args = expr_R(NULL), time_steps, ..., plot = FALSE){

  time_steps <- enquo(time_steps)
  time_steps <- quo_updt(time_steps, as_name(args[["ion1"]]))
  gr_by <- enquos(...)

  df.acf <- Rm(.df, args = args, !!! gr_by, output = "augment") %>%
    left_join(.df, ., by = "ID") %>%
    group_by(!!! gr_by) %>%
    arrange(!! time_steps) %>%
    select(!!! gr_by, .data$.std.resid) %>%
    mutate(LB_test = stats::Box.test(.data$.std.resid, type = "Ljung-Box")$p.value,
           flag_IR = if_else(
             LB_test < 0.05,
             "Ha (dependence of residuals)",
             "H0 (independence of residuals)"
           ),
           flag_IR = as.factor(flag_IR)
           ) %>%
    tidyr::nest(data = .data$.std.resid) %>%
    mutate(ACF = purrr::map(data, acf_calc)) %>%
    select(-data) %>%
    ungroup() %>%
    tidyr::unnest(cols = c(ACF))

  if (plot) {

  df.acf %T>%
      { print(ggplot(., mapping = aes(x = lag, y = acf)) +
                geom_hline(aes(yintercept = 0)) +
                geom_segment(mapping = aes(xend = lag, yend = 0)) +
                geom_hline(aes(yintercept = ci_upper), color = "darkblue") +
                geom_hline(aes(yintercept = ci_lower), color = "darkblue") +
                facet_wrap(vars(!!! gr_by), scales = "free") +
                ggtitle("ACF plot") +
                theme_classic()
              )
        }

  } else {

    return(df.acf)
  }


}





#-------------------------------------------------------------------------------
# Not exportet helper functions
#-------------------------------------------------------------------------------

# calculate autocorrelation
acf_calc <- function(data){
  acf  <- acf(data$.std.resid, plot = FALSE)
  ci <- qnorm((1 - 0.95) / 2) / sqrt(length(data$.std.resid))
  df.acf <- tibble(
    lag = as.vector(acf$lag)[-1],
    acf = as.vector(acf$acf)[-1],
    ci_upper = ci,
    ci_lower = -ci
  )
}





# nest lm (args as quos)
nest_R_lm <- function(df, gr_by, Xt1, Xt2, t, method, .hyp){


  #group_by(df, ) %>%
  tidyr::nest(df, t = !!t, data = -c(!!! gr_by)) %>%
  mutate(
    R_lm =
      purrr::map(data, ~lm_form(.x, Xt1, Xt2, type = "Rm")),
    aug =
        purrr::map(R_lm, broom::augment),
    extr =
      purrr::map(aug, ~transmute_reg(.x, Xt1, Xt2, method)),
    extr =
      purrr::map(extr, ~QQ_trans(.x, method, .hyp = .hyp)),
    extr =
      purrr::map2(aug, extr, ~bp_wrap(.x, .y, Xt2, method, .hyp = .hyp)),
    flag =
      purrr::map2(extr, aug, ~flag_set(.x, .y, Xt1, Xt2, type = method))
    ) %>%
  select(-c(R_lm, aug)) #%>%
    #ungroup()


}


# augment function transform  and rename variables to standards of point
transmute_reg <- function(df, Xt1, Xt2, type) {

  # predicted variable and sigma level boundaries
  prefix <- c("hat", "hat_sl", "hat_su")
  hat_args <- purrr::map(prefix, ~quo_updt(Xt1, pre = .x))
  hat_nm <- paste(prefix, as_name(Xt1), sep = "_")

  # predicted sigma level boundaries studentized residuals
  hat_args_studE_nm <- paste(c("hat", "hat_sl", "hat_su"), "studE", sep = "_")
  hat_args_studE <- parse_exprs(hat_args_studE_nm)


  # model args
  args <- quos(
    !!hat_args[[1]] := .data$.fitted,
    !!hat_args[[2]] := .data$.fitted - 2 * sigma_calc(.data$.resid),
    !!hat_args[[3]] := .data$.fitted + 2 * sigma_calc(.data$.resid) ,
    studE = .data$.std.resid,
    hat_Xi = .data$.hat,
    CooksD = .data$.cooksd,
    !!hat_args_studE[[1]] := 0,
    !!hat_args_studE[[2]] := -3.5,
    !!hat_args_studE[[3]] :=  3.5
    )

  if (type == "Rm") args <- args[names(args) %in% hat_nm]
  if (type == "norm_E") args <- args[!names(args) %in% hat_nm]
  if (type == "CooksD") args <- args[names(args) %in% c(hat_nm, "CooksD")]
  if (type == "CV") args <- args[names(args) %in% c(hat_nm[1], "studE", hat_args_studE_nm)]
  if (type == "QQ") args <- args["studE"]

  transmute(df, !!!args)


}

# create flag variable
flag_set <- function(df1, df2, Xt1, Xt2, type){

  # predicted variable and sigma level boundaries
  prefix <- c("hat", "hat_sl", "hat_su")
  hat_args_sigma <- purrr::map(prefix, ~quo_updt(Xt1, pre = .x))

  # predicted variable and standard error and CI boundaries
  hat_args_CI <- paste(c("hat", "hat_cl", "hat_cu"), "RQ", sep = "_") %>%
    parse_exprs()

  # predicted sigma level boundaries studentized residuals
  hat_args_studE <- paste(c("hat", "hat_sl", "hat_su"), "studE", sep = "_") %>%
    parse_exprs()


  if (type == "Rm") {
     df <- mutate(df1, !! Xt1 := pull(df2, !! Xt1)) %>%
       transmute(flag = as.factor(if_else(between(!! Xt1, !! hat_args_sigma[[2]], !! hat_args_sigma[[3]]), "good", "bad")))
    return(df)
    }
  if (type == "CooksD" | type == "norm_E") {
    df <- transmute(df1, flag = as.factor(if_else(CooksD < {4 / (n() - 2)}, "good", "bad")))
    return(df)
    }
  if (type == "QQ") {
    df <- transmute(df1, flag = as.factor(if_else(!between(RQ, !!hat_args_CI[[2]], !!hat_args_CI[[3]]), "good", "bad")))
    return(df)
    }
  if (type == "CV") {
    df <- transmute(df1, flag = as.factor(if_else(between(studE, !!hat_args_studE[[2]], !!hat_args_studE[[3]]), "good", "bad")))
    return(df)
    }

}

# quantile transformations and hypothesis tests
QQ_trans <- function(df, type, .hyp) {

  if (type!= "QQ") return(df)

  # predicted variable and standard error and CI boundaries
  hat_args <- paste(c("hat", "hat_cl", "hat_cu"), "RQ", sep = "_") %>%
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

  # if (.hyp == "none") hyp_result <- NULL; Ha = NULL; H0 = NULL

  df <- transmute(
    df,
    RQ = unname(quantile(studE, probs = vector_probs(n()))),
    # Calculate normal (Theoretical) quantiles using mean and standard deviation
    TQ = qnorm(vector_probs(n()), mean(RQ), sd(RQ)),
    # The standard error
    !!hat_args[[1]] := mean(RQ) + sd(RQ) * TQ,
    !!hat_args[[2]] := !!hat_args[[1]] - 2 * hat_QR_se(RQ, TQ, vector_probs(n()), n()),
    !!hat_args[[3]] := !!hat_args[[1]] + 2 * hat_QR_se(RQ, TQ, vector_probs(n()), n())
    )

  if (.hyp != "none") {
    return(mutate(df, hyp = if_else(hyp_result < 0.05, Ha, H0)))
    } else {
      return(df)
      }

}

# Hetroscadasticity test (Breusch Pagan test)(level of confidence 95%;
# cut-off 0.05 for H0 rejection)
bp_wrap <- function(df1, df2, Xt2, type, .hyp){

  # Breusch Pagan test
  if (type == "CV" & .hyp == "bp") {
    Chi_R2 <- custom_bp(df1, Xt2)
    Ha <- "Ha (heteroskedasticity)"
    H0 <- "H0 (homoskedasticity)"
    # return(mutate(df2, hyp = Chi_R2))
    return(mutate(df2, hyp = if_else(Chi_R2 > qchisq(.95, df = 1), Ha, H0)))
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

# use the formula i - 0.5/ in, for i = 1,..,n
# this is a vector of the n probabilities (theoretical cumulative distribution function CDF)
vector_probs <- function(n){
  ((1:unique(n)) - 0.5) / (unique(n))
}

# standard error of quantiles model
hat_QR_se <- function(RQ, TQ, pb, n){
  (sd(RQ) / dnorm(TQ)) * sqrt((pb * (1 - pb))/ unique(n))
}
# confidence interval regression model
hat_Y_se <- function(sigma, hat_Xi){
  sigma * sqrt(hat_Xi)
}
