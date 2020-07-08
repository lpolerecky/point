#' Ideal linear R Model
#'
#' @export
Rm <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

# variable names of original dataset
  var_names <- colnames(df)
  var_names <- var_names[!var_names %in% "ID"]

  df %>%
    group_by(!!! gr_by) %>%
    mutate(
# Hat values
            hat_Xi = stats::hat(!! quo_updt(my_q = args[["Xt"]],
                                            txt = as_name(args[["ion2"]])
                                            )
                                ),
# Modelled Y values
            hat_Y = !! quo_updt(my_q = args[["Xt"]],
                                x = "M_R"
                                ) *
                    !! quo_updt(my_q = args[["Xt"]],
                                txt = as_name(args[["ion2"]])
                                ),
# residuals
            E = !! quo_updt(my_q = args[["Xt"]],
                            txt = as_name(args[["ion1"]])) -
                hat_Y,
# Standard error of the regression
            sigma = sigma_calc(E),
            hat_min = hat_Y - 2 * sigma,
            hat_max = hat_Y + 2 * sigma,
            flag = if_else(!! quo_updt(my_q = args[["Xt"]],
                                       txt = as_name(args[["ion1"]])
                                       ) > hat_min &
                           !! quo_updt(my_q = args[["Xt"]],
                                       txt = as_name(args[["ion1"]])
                                       ) < hat_max,
                           "good",
                           "bad"
                           )
            ) %>%
    ungroup() %>%
    eval_tidy(expr = mod_out(output, var_names))

}


# Standard error of the regression
sigma_calc <- function(res) sqrt(sum((res ^ 2)) / (length(res) - 1))

# Switch output complete dataset, stats or summary stats
mod_out <- function(output, vars = NULL) {
  switch(output,
         flag = call2("select",
                      expr(.),
                      !!! parse_exprs(paste0("-", vars))
                      ),
         complete = call2( "invisible", expr(.)),
         normE = call2("select",
                       expr(.),
                       expr(.data$ID),
                       expr(.data$hat_Y),
                       expr(.data$hat_Xi),
                       expr(.data$E)
                       ),
         CooksD =call2("select",
                       expr(.),
                       expr(.data$ID),
                       expr(.data$studE),
                       expr(.data$hat_Xi),
                       expr(.data$flag)
                       ),
         Rm = call2("select",
                    expr(.),
                    expr(.data$ID),
                    expr(.data$hat_Y),
                    expr(.data$hat_min),
                    expr(.data$hat_max)
                    ),
         studE = call2("select",
                       expr(.),
                       expr(.data$ID),
                       expr(.data$studE),
                       expr(.data$hat_Y),
                       expr(.data$hat_Xi)
                      )
         )
  }


#' Normalization of residuals
#'
#' @export
norm_E <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

# variable names of original dataset
  var_names <- colnames(df)
  var_names <- var_names[!var_names %in% "ID"]

  Rm(df, args = args, !!! gr_by, output = "normE") %>%
    left_join(df, ., by = "ID") %>%
    group_by(!!! gr_by) %>%
    tidyr::nest() %>%
    mutate(
# jackknifed mean R
      M_Ri = purrr::map(
        data,
        ~jack_meanR(df = .x,
                    ion1 = !! quo_updt(my_q = args[["Xt"]],
                                       txt = as_name(args[["ion1"]])
                                       ),
                    ion2 = !! quo_updt(my_q = args[["Xt"]],
                                       txt = as_name(args[["ion2"]]))
                    )
        )
      ) %>%
    tidyr::unnest(cols = c(data, M_Ri)) %>%
    mutate(
# jackknifed modelled Y values
      hat_Yi = purrr::map2_dbl(
        M_Ri,
        !! quo_updt(my_q = args[["Xt"]],
                    txt = as_name(args[["ion2"]])),
        ~{.x * .y}
        ),
# residuals with i-th value removed
      Ei = purrr::map2_dbl(
        !! quo_updt(my_q = args[["Xt"]],
                    txt = as_name(args[["ion1"]])),
        hat_Yi,
        ~{.x - .y}
        )
      ) %>%
    tidyr::nest() %>%
    mutate(
# standard error of regression (external with i-th residual removed)
      sigma_i = purrr::map(data, ~jack_sigma(.x, Ei))
        ) %>%
    tidyr::unnest(cols = c(data, sigma_i)) %>%
    mutate(
      studE = purrr::pmap_dbl(
        lst(res = E, sigma_i = sigma_i, hat_Xi = hat_Xi),
        studE_calc
        ),
      CooksD = purrr::map2_dbl(studE, hat_Xi, cooksD_calc),
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


# calculation of studentized residuals
studE_calc <- function(res, sigma_i, hat_Xi) res / (sigma_i * sqrt(1 - hat_Xi))

# jackknifed mean R values
jack_meanR <- function(df, ion1, ion2){

  Xt_ion1 <- enquo(ion1)
  Xt_ion2 <- enquo(ion2)

  Xt_ion1 <- select(df, !! Xt_ion1) %>% pull(!! Xt_ion1)
  Xt_ion2 <- select(df, !! Xt_ion2) %>% pull(!! Xt_ion2)

  c((resample::jackknife(Xt_ion1 , mean))$replicates /
    (resample::jackknife(Xt_ion2 , mean))$replicates)

}


# jackknifed standard error of the regression
jack_sigma <- function(df, res){

  res <- enquo(res)
  res <- select(df, !! res) %>% pull(!! res)
  c((resample::jackknife(res, sigma_calc))$replicates)

}


#' @rdname Cameca_R
#'
#' @export
CooksD_R <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

# Variable names of original dataset
  var_names <- colnames(df)
  var_names <- var_names[!var_names %in% "ID"]

  df_cd <- norm_E(df, args = args, !!! gr_by, output = "CooksD") %>%
    left_join(df, ., by = "ID")

  Rm(df_cd, args = args, !!! gr_by, output = "Rm") %>%
    left_join(df_cd, ., by = "ID") %>%
    eval_tidy(expr =  mod_out(output, var_names))
  }

#' Normality test
#'
#' @export
QQ <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

# Variable names of original dataset
  var_names <- colnames(df)
  var_names <- var_names[!var_names %in% "ID"]

  norm_E(df, args = args, !!! gr_by, output = "studE") %>%
    left_join(df, ., by = "ID") %>%
    group_by(!!! gr_by) %>%
    mutate(RQ = studE) %>%
    arrange(RQ) %>%
    mutate(prob_vc =  vector_probs(!! quo_updt(my_q = args[["Xt"]],
                                               x = "n_R"
                                               )
                                  ),
           RQ = unname(quantile(studE, probs = prob_vc)),
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
             (nortest::ad.test(studE))$p.value < 0.05,
             "Ha (non-normal)",
             "H0 (normal)"
             ),
# t-test flag for mu0 (aka the conditional mean of epsilon) being zero
           flag_CM = if_else(
             (t.test(studE, mu = 0))$p.value  < 0.05,
             "Ha (mu0 is not zero)",
             "H0 (mu0 is zero)"
             )
           ) %>%
           ungroup() %>%
           eval_tidy(expr =  mod_out(output, var_names))
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

#' Constant variance test
#'
#' @export
CV <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

  # Variable names of original dataset
  var_names <- colnames(df)
  var_names <- var_names[!var_names %in% "ID"]

  norm_E(df, args = args, !!! gr_by, output = "studE") %>%
    left_join(df, ., by = "ID") %>%
    group_by(!!! gr_by) %>%
# hetroscadasticity test (Breusch Pagan test)(level of confidence 95%;
# cut-off 0.05 for H0 rejection),
    tidyr::nest() %>%
    mutate(
      res_lm = purrr::map(data, ~lm_res(.x, args = args)),
      R2 = purrr::map(res_lm, ~(summary(.x)$r.squared)),
      SE_beta = purrr::map(res_lm, ~(unname(summary(.x)$coefficients[2,2]))),
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
           fitted = hat_Y,
           hat_Y = 0,
           flag = if_else(studE > hat_min & studE < hat_max,
                          "good",
                          "bad"
                          )
           ) %>%
    ungroup() %>%
    eval_tidy(expr =  mod_out(output, var_names))
}


# lm residuals (this can be done better!! tidy formula call)
lm_res <- function(data, args){
  call_lm <- parse_expr(paste0("data$studE~data$", as_name(args[["Xt"]]), ".", as_name(args[["ion2"]])))
  eval_tidy(call2("lm", call_lm))

}

#' Indpendence of residuals
#'
#' @export
IR <- function(df, args = expr_R(NULL), time_steps, ..., plot = FALSE){

  time_steps <- enquo(time_steps)
  time_steps <- quo_updt(time_steps, as_name(args[["ion1"]]))
  gr_by <- enquos(...)


  df.acf <- norm_E(df, args = args, !!! gr_by, output = "studE") %>%
    left_join(df, ., by = "ID") %>%
    group_by(!!! gr_by) %>%
    arrange(!! time_steps) %>%
    select(!!! gr_by, studE) %>%
    mutate(LB_test = stats::Box.test(studE, type = "Ljung-Box")$p.value,
           flag_IR = if_else(
             LB_test < 0.05,
             "Ha (dependence of residuals)",
             "H0 (independence of residuals)"
           )
           ) %>%
    tidyr::nest(data = studE) %>%
    mutate(ACF = purrr::map(data, acf_calc)) %>%
    select(-data) %>%
    ungroup() %>%
    tidyr::unnest(cols = c(ACF))

  if (plot){

  df.acf %T>%
      {print(ggplot(., mapping = aes(x = lag, y = acf)) +
               geom_hline(aes(yintercept = 0)) +
               geom_segment(mapping = aes(xend = lag, yend = 0)) +
               geom_hline(aes(yintercept = ci_upper), color = "darkblue") +
               geom_hline(aes(yintercept = ci_lower), color = "darkblue") +
               facet_wrap(vars(!!! gr_by), scales = "free") +
               ggtitle("ACF plot") +
               theme_classic())}

  } else {

    return(df.acf)
  }


}


acf_calc <- function(data){
  acf  <- acf(data$studE, plot = FALSE)
  ci <- qnorm((1 - 0.95) / 2) / sqrt(length(data$studE))
  df.acf <- tibble(lag = as.vector(acf$lag)[-1],
                   acf = as.vector(acf$acf)[-1],
                   ci_upper = ci,
                   ci_lower = -ci
                   )
  }

# https://stackoverflow.com/questions/27264266/multiple-ggplots-with-magrittr-tee-operator
ggpass <- function(gg){
  print(gg)
  return(gg$data)
}
