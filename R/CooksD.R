#' Ideal linear R model
#'
#' @export
Rm <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

# heavy isotope
  Xt1 <- quo_updt(args[["Xt"]], as_name(args[["ion1"]]))
# light isotope
  Xt2 <- quo_updt(args[["Xt"]], as_name(args[["ion2"]]))

# standard deviation of R
  S_R <- quo_updt(my_q = args[["Xt"]], x = "S_R")

# variable names of original dataset
  var_names <- c(colnames(df), aug_names)
  var_names <- var_names[!var_names %in% "ID"]

  df %>%
    group_by(!!! gr_by) %>%
    tidyr::nest() %>%
    mutate(R_lm = purrr::map(data, ~lm_form(.x, Xt1, Xt2, type = "Rm")),
           model = purrr::map(R_lm, {~broom::augment(.x) %>% select(-c(!!Xt1, !!Xt2))})) %>%
    tidyr::unnest(cols = c(data, model)) %>%
    mutate(
# Modelled Y values
            hat_Y = .data$.fitted,
            hat_min = hat_Y - 2 * sigma_calc(.data$.resid),
            hat_max = hat_Y + 2 * sigma_calc(.data$.resid),
            flag = if_else(!! Xt1 > hat_min &
                           !! Xt1 < hat_max,
                           "good",
                           "bad"
                           )
            ) %>%
    ungroup() %>%
    eval_tidy(expr = mod_out(output, var_names))

}

#' Normalization of residuals
#'
#' @export
norm_E <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

# variable names of original dataset
  var_names <- c(colnames(df), aug_names)
  var_names <- var_names[!var_names %in% "ID"]

  Rm(df, args = args, !!! gr_by, output = "augment") %>%
    left_join(df, ., by = "ID") %>%
    group_by(!!! gr_by) %>%
    mutate(
      studE = .data$.std.resid,
      hat_Xi = .data$.hat,
      CooksD = .data$.cooksd,
      CooksD_cf = purrr::map_dbl(
        !! quo_updt(my_q = args[["Xt"]], x = "n_R"),
        cooksD_cut_calc
       ),
      flag = if_else(CooksD > CooksD_cf,
                     "bad",
                     "good"
      )
      ) %>%
    ungroup() %>%
    eval_tidy(expr =  mod_out(output, var_names))
  }


#' @rdname Cameca_R
#'
#' @export
CooksD <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

  # variable names of original dataset
  var_names <- c(colnames(df), aug_names)
  var_names <- var_names[!var_names %in% "ID"]

  df_cd <- norm_E(df, args = args, !!! gr_by, output = "flag")

  Rm(df, args = args, !!! gr_by, output = "complete") %>%
    select(-flag) %>%
    left_join(df_cd, ., by = "ID") %>%
    eval_tidy(expr =  mod_out(output, var_names))
}

#' Normality test
#'
#' @export
QQ <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

  # variable names of original dataset
  var_names <- c(colnames(df), aug_names)
  var_names <- var_names[!var_names %in% "ID"]

  Rm(df, args = args, !!! gr_by, output = "augment") %>%
     left_join(df, ., by = "ID") %>%
     group_by(!!! gr_by) %>%
     mutate(prob_vc =  vector_probs(!! quo_updt(my_q = args[["Xt"]],
                                               x = "n_R"
                                               )
                                  ),
            RQ = unname(quantile(.data$.std.resid, probs = prob_vc)),
# Calculate normal (Theoretical) quantiles using mean and standard deviation from
            TQ = qnorm(prob_vc, mean(RQ), sd(RQ)),
# The standard error is calculated with,
            hat_Y = mean(RQ) + sd(RQ) * TQ,
            hat_RQ_se = hat_QR_se(RQ,
                                  TQ,
                                  prob_vc,
                                  !! quo_updt(my_q = args[["Xt"]],
                                              txt = as_name(args[["ion2"]]),
                                              x = "n"
                                              )
                                  ),
            hat_min =  hat_Y - 2 * hat_RQ_se,
            hat_max =  hat_Y + 2 * hat_RQ_se,
            flag = if_else(RQ > hat_min & RQ < hat_max, "good", "bad"),
            flag_QQ = if_else(
              (nortest::ad.test(.data$.std.resid))$p.value < 0.05,
              "Ha (non-normal)",
              "H0 (normal)"
              ),
# t-test flag for mu0 (aka the conditional mean of epsilon) being zero
            flag_CM = if_else(
              (t.test(.data$.std.resid, mu = 0))$p.value  < 0.05,
              "Ha (mu0 is not zero)",
              "H0 (mu0 is zero)"
              )
            ) %>%
            ungroup() %>%
            eval_tidy(expr =  mod_out(output, var_names))
}


#' Constant variance test
#'
#' @export
CV <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

# heavy isotope
  Xt2 <- quo_updt(args[["Xt"]], as_name(args[["ion2"]]))

# variable names of original dataset
  var_names <- c(colnames(df), aug_names)
  var_names <- var_names[!var_names %in% "ID"]

  Rm(df, args = args, !!! gr_by, output = "augment") %>%
    left_join(df, ., by = "ID") %>%
    group_by(!!! gr_by) %>%
# Hetroscadasticity test (Breusch Pagan test)(level of confidence 95%;
# cut-off 0.05 for H0 rejection),
    tidyr::nest() %>%
    mutate(
      res_lm = purrr::map(data, ~lm_form(.x, quo(.std.resid), Xt2)),
      R2 = purrr::map(res_lm, ~{broom::glance(.x)$r.squared}),
      SE_beta = purrr::map(res_lm, ~{pull(broom::tidy(.x), std.error)[2]}),
          ) %>%
    tidyr::unnest(cols = c(data, R2, SE_beta)) %>%
    mutate(Chi_R2 = R2 * !! quo_updt(my_q = args[["Xt"]],
                                      txt = as_name(args[["ion2"]]),
                                      x = "n"
                                      ),
           flag_CV = if_else(Chi_R2 >
                        qchisq(.95, df = 1),
                        "Ha (heteroskedasticity)",
                        "H0 (homoskedasticity)"
                            ),
           hat_min = -3.5,
           hat_max = 3.5,
           fitted = .data$.fitted,
           studE = .data$.std.resid,
           hat_Y = 0,
           flag = if_else(studE  > hat_min &
                            studE < hat_max,
                          "good",
                          "bad"
                          )
           ) %>%
    ungroup() %>%
    eval_tidy(expr =  mod_out(output, var_names))
}


#' Independence of residuals
#'
#' @export
IR <- function(df, args = expr_R(NULL), time_steps, ..., plot = FALSE){

  time_steps <- enquo(time_steps)
  time_steps <- quo_updt(time_steps, as_name(args[["ion1"]]))
  gr_by <- enquos(...)

  df.acf <- Rm(df, args = args, !!! gr_by, output = "augment") %>%
    left_join(df, ., by = "ID") %>%
    group_by(!!! gr_by) %>%
    arrange(!! time_steps) %>%
    select(!!! gr_by, .data$.std.resid) %>%
    mutate(LB_test = stats::Box.test(.data$.std.resid, type = "Ljung-Box")$p.value,
           flag_IR = if_else(
             LB_test < 0.05,
             "Ha (dependence of residuals)",
             "H0 (independence of residuals)"
           )
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


#' Standard error of the regression
#' @export
sigma_calc <- function(res) sqrt(sum((res ^ 2)) / (length(res) - 1))


#-------------------------------------------------------------------------------
# Not exportet helper functions
#-------------------------------------------------------------------------------
aug_names <- c(".cooksd",
               ".fitted",
               ".hat",
               ".resid",
               ".sigma",
               ".std.resid"
               )

# calculate autocorrelation
acf_calc <- function(data){
  acf  <- acf(data$.std.resid, plot = FALSE)
  ci <- qnorm((1 - 0.95) / 2) / sqrt(length(data$.std.resid))
  df.acf <- tibble(lag = as.vector(acf$lag)[-1],
                   acf = as.vector(acf$acf)[-1],
                   ci_upper = ci,
                   ci_lower = -ci
  )
}


# Switch output to flag, augmentation and complete
mod_out <- function(output, vars = NULL) {
  switch(output,
         flag = call2("select",
                      expr(.),
                      # expr(!contains("."))
                      !!! parse_exprs(paste0("-", vars))
                      ),
         complete = call2( "invisible", expr(.)),
         augment =  call2("select",
                          expr(.),
                          expr(starts_with(c("ID", ".")))
         )
  )
}


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

# Cook's D cutt-off value calculator
cooksD_cut_calc  <- function(n) 4 / (n - 2)
