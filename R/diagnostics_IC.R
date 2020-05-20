#' Diagnostics isotope count data
#'
#' \code{diag_R} wrapper function for diagnostics on isotope count data
#'
#' The \code{diag_R} function performs an internal call to stat_R to perform
#' diagnostics on the influence of individial measurements on the blockwise or
#' global (i.e., a complete analysis) stastics. It identifies potentially
#' influential measurements that can subsequentally be filtered from the
#' original dataset to improve the precision. See
#' \code{vignette("IC-diagnostics", package = "point")} for more information
#' on how to use the function, and possible methods for diagnostics.
#'
#' @param df A tibble containing processed ion count data.
#' @param method Character string for the type of diagnostics. "Cameca" pertains
#' to the default Cameca software setting. "CooksD" pertains to global
#' regression diagnostics based on Cook's D statistics.
#' @param args A list of quosures pertaining to the variables required for a
#' call to stat_R. The function expr_R can be used to simplify setting this
#' argument.
#' @param ... Variables for grouping.
#' @param output Character string determing whether the returned values in a
#' minimal version `"flag"` (original dataset + diagnostics) or an extended
#' version with all the intermediate steps of ion- and isotope-wise summary
#' statistics `"complete"`.
#'
#' @return A t\code{\link[tibble:tibble]{tibble}} containing either the original
#' dataset with new columns related to the diagnostics an extended version with
#' all the intermediate steps of ion- and isotope-wise summary statistics. The
#' flag variable enables convenient filtering of the original tibble for an
#' augmentation of the original dataset.
#' @export
#' @examples
#' # Use point_example() to access the examples bundled with this package
#'
#' # raw data containing 13C and 12C counts on carbonate
#' tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # Processing raw ion count data
#' tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt)
#'
#' # CAMECA style augmentatio of ion count data for isotope ratios
#' tb.aug <- diag_R(tb.pr,
#'                  method = "Cameca",
#'                  args = expr_R(Xt = "Xt.pr",
#'                                N = "N.pr",
#'                                species = "species.nm",
#'                                ion1 = "13C",
#'                                ion2 = "12C"),
#'                  file.nm,
#'                  bl.mt)
diag_R <- function(df, method = "Cameca", args = expr_R(NULL), ...,
                   output = "flag"){

  gr_by <- enquos(...)

# ID for connecting flag to original dataframe
  df <- ID_builder(df, !! args[["species"]], !!! gr_by)

# Method selection
  diag_method <- function(method){
    switch(method,
           Cameca = call2("Cameca_R", expr(.), args, !!!gr_by, output = output),
           CooksD = call2("CooksD_R", expr(.), args, !!!gr_by, output = output)
           )
  }

# # Remove zeros
#   df <- zeroCt(df,
#                !! args[["N"]],
#                !! args[["species"]],
#                as_name(args[["ion1"]]),
#                as_name(args[["ion2"]]),
#                !!! gr_by,
#                warn = FALSE
#                )

# Descriptive an predictive statistics for ion ratios
  tb.aug <- stat_R(df,
                   Xt = !! args[["Xt"]],
                   N =!! args[["N"]],
                   species = !! args[["species"]],
                   ion1 = as_name(args[["ion1"]]),
                   ion2 = as_name(args[["ion2"]]),
                   !!! gr_by,
                   output = "complete",
                   zero = FALSE
                   ) %>%
              eval_tidy(expr = diag_method(method))

# Datafile with flag values associated to diagnostics
  if (output == "flag"){
    return(left_join(df, tb.aug, by = "ID") %>% select(-.data$ID))
    }
   if (output == "complete"){
     return(tb.aug)
   }
  }

#' Family of diagnostics functions for isotope count ratios
#'
#' \code{Cameca_R} default CAMECA diagnostics
#' \code{CooksD_R} regression diagnostics based on Cook's D
#'
#' These functions perform a specific set of diagnostics to term anomolous
#' values inr raw ion count data of an isotope pair. The wrapper function
#' \code{Diag_R} is more convenient as it defines all the ion- and isotope-wise
#' statistics required for the diagnostics.
#'
#' @param df A tibble with ion count data and statistics for ion- and
#' isotope-wise statistics.
#' @param args A list of quosures pertaining to the variables required for a
#' call to stat_R. The function expr_R can be used to simplify setting this
#' argument.
#' @param ... Variables for grouping.
#' @param output Character string determing whether the returned values in a
#' minimal version `"flag"` (original dataset + diagnostics) or an extended
#' version with all the intermediate steps of ion- and isotope-wise summary
#' statistics `"complete"`.
#'
#' @return A t\code{\link[tibble:tibble]{tibble}} containing either the original
#' dataset with new columns related to the diagnostics or only the diagnostics.
#' The flag variable enables convenient filtering of the original tibble for an
#' augmentation of the original dataset.
#' @export
#' @examples
#' # Use point_example() to access the examples bundled with this package
#'
#' # raw data containing 13C and 12C counts on carbonate
#' tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # Processing raw ion count data
#' tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt)
#'
#' # Descriptive an predictive statistics for 13C/12C ratios (note output
#' # argument)
#' tb.R <- stat_R(tb.pr, Xt.pr, N.pr, species.nm, ion1 = "13C",
#'                ion2 = "12C", file.nm, bl.mt, output = "complete")
#'
#' # CAMECA style augmentatio of ion count data for isotope ratios
#' tb.aug <- Cameca_R(tb.R,
#'                    args = expr_R(Xt = "Xt.pr",
#'                                  N = "N.pr",
#'                                  species = "species.nm",
#'                                  ion1 = "13C",
#'                                  ion2 = "12C"),
#'                   file.nm,
#'                   bl.mt,
#'                   output = "flag")
Cameca_R <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

  # Switch output complete dataset, stats or summary stats
  mod_out <- function(output) {
    switch(output,
           flag = call2( "select", expr(.),
                         expr(.data$ID),
                         expr(.data$lower),
                         expr(.data$upper),
                         expr(.data$flag)
                         ),
           complete = call2( "invisible", expr(.)))
    }

  df %>%
    group_by(!!! gr_by) %>%
    mutate(lower = !! quo_updt(my_q = args[["Xt"]],
                               x = "M_R") - 2 *
             !! quo_updt(my_q = args[["Xt"]],
                         x = "S_R"),
           upper = !! quo_updt(my_q = args[["Xt"]],
                               x = "M_R") + 2 *
             !! quo_updt(my_q = args[["Xt"]],
                         x = "S_R"),
           flag = if_else(
             between(!! quo_updt(my_q = args[["Xt"]], x = "R"),
                     unique(.data$lower),  # mean - 2SD
                     unique(.data$upper)), # mean + 2SD
             "non-influential",
             "influential"
             )
           ) %>%
    ungroup() %>%
    eval_tidy(expr =  mod_out(output))
  }

#' @rdname Cameca_R
#'
#' @export
CooksD_R <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

# Switch output complete dataset, stats or summary stats
  mod_out <- function(output) {
    switch(output,
           flag = call2("select",
                        expr(.),
                        expr(.data$ID),
                        expr(.data$E),
                        expr(.data$sigma),
                        expr(.data$hat_Y),
                        expr(.data$hat_Y_min),
                        expr(.data$hat_Y_max),
                        expr(.data$hat_Xi),
                        expr(.data$studE),
                        expr(.data$CooksD),
                        expr(.data$CooksD_cf),
                        expr(.data$flag),
                        expr(.data$RQ),
                        expr(.data$TQ),
                        expr(.data$hat_RQ),
                        expr(.data$hat_RQ_min),
                        expr(.data$hat_RQ_max),
                        expr(.data$flag_QQ),
                        expr(.data$flag_CM),
                        expr(.data$R2),
                        expr(.data$SE_beta),
                        expr(.data$Chi_R2),
                        expr(.data$flag_CV),
                        expr(.data$ACF)
                         ),
           complete = call2( "invisible", expr(.)),
    )
    }

  df %>%
    group_by(!!! gr_by) %>%
    mutate(
# Hat values
           hat_Xi = stats::hat(!! quo_updt(my_q = args[["Xt"]],
                                            txt = as_name(args[["ion2"]])
                                           )
                               ),
# modelled Y values
           hat_Y = !! quo_updt(my_q = args[["Xt"]],
                               x = "M_R"
                               ) *
                   !! quo_updt(my_q = args[["Xt"]],
                               txt = as_name(args[["ion2"]])
                               ),

# residuals
           E = !! quo_updt(my_q = args[["Xt"]],
                           txt = as_name(args[["ion1"]])
                           ) - hat_Y,
# Standard error of the regression
           sigma = sigma_calc(E),
# 95 CI of the regression
           hat_Y_min = hat_Y - 2 * hat_Y_se(sigma, hat_Xi),
           hat_Y_max = hat_Y + 2 * hat_Y_se(sigma, hat_Xi),
# Mean Square Error
           MSE = sum(E ^ 2) / !! quo_updt(my_q = args[["Xt"]],
                                          x = "n_R"
                                          )
            ) %>%
    tidyr::nest() %>%
    mutate(
# jackknifed mean R
           M_Ri = purrr::map(data,
                             ~jack_meanR(df = .x,
                                         ion1 = !! quo_updt(my_q = args[["Xt"]],
                                                            txt = as_name(args[["ion1"]])),
                                         ion2 = !! quo_updt(my_q = args[["Xt"]],
                                                            txt = as_name(args[["ion2"]]))
                                        )
                            )

           ) %>%
    tidyr::unnest(cols = c(data, M_Ri)) %>%
    mutate(
# jackknifed modelled Y values
           hat_Yi = purrr::map2_dbl(M_Ri,
                                    !! quo_updt(my_q = args[["Xt"]],
                                                txt = as_name(args[["ion2"]])),
                                    ~{.x * .y}),
# residuals with i-th value removed
           Ei = purrr::map2_dbl(!! quo_updt(my_q = args[["Xt"]],
                                            txt = as_name(args[["ion1"]])),
                                hat_Yi,
                                ~{.x - .y})
          ) %>%
    tidyr::nest() %>%
    mutate(
# standard error of regression (external with i-th residual removed)
          sigma_i = purrr::map(data, ~jack_sigma(.x, Ei))
          ) %>%
    tidyr::unnest(cols = c(data, sigma_i)) %>%
    mutate(
           studE = purrr::pmap_dbl(lst(res = E,
                                       sigma_i = sigma_i,
                                       hat_Xi = hat_Xi
                                       ),
                                   studE_calc
                                   ),
           CooksD = purrr::map2_dbl(studE, hat_Xi, cooksD_calc),
           CooksD_cf = purrr::map_dbl(!! quo_updt(my_q = args[["Xt"]],
                                                  x = "n_R"
                                                  ),
                                      cooksD_cut_calc
                                      ),
           flag = if_else(CooksD > CooksD_cf, "influential", "non-influential")
          ) %>%
# normality test
    # mutate(RQ = studE) %>%
    # arrange(RQ) %>%
# use the formula i - 0.5/ in, for i = 1,..,n for probs
# this is a vector of the n probabilities ( theoretical cumulative distribution function CDF)
   mutate(prob_vc =  vector_probs(!! quo_updt(my_q = args[["Xt"]],
                                              txt = as_name(args[["ion2"]]),
                                              x = "n"
                                              )
                                   ),
          RQ = unname(quantile(studE, probs = prob_vc)),
# calculate normal (Theoretical) quantiles using mean and standard deviation from
          TQ = qnorm(prob_vc, mean(RQ), sd(RQ)),
# the standard error is calculated with,
          hat_RQ = mean(RQ) + sd(RQ) * TQ,
          hat_RQ_se = hat_QR_se(RQ,
                                TQ,
                                prob_vc,
                                !! quo_updt(my_q = args[["Xt"]],
                                            txt = as_name(args[["ion2"]]),
                                            x = "n"
                                            )
                                ),
          hat_RQ_min =  hat_RQ - 2 * hat_RQ_se,
          hat_RQ_max =  hat_RQ + 2 * hat_RQ_se,
          flag_QQ = if_else(
            (nortest::ad.test(studE))$p.value < 0.05,
            "Ha (non-normal)",
            "H0 (normal)"
            ),
# t-test flag for mu0 (aka the conditional mean of epsilon) being zero
          flag_CM = if_else((t.test(studE, mu = 0))$p.value  < 0.05,
                            "Ha (mu0 is not zero)",
                            "H0 (mu0 is zero)"
                            )
          ) %>% #,
# hetroscadasticity test (Breusch Pagan test)(level of confidence 95%;
# cut-off 0.05 for H0 rejection),
    tidyr::nest() %>%
    mutate(
          res_lm = purrr::map(data, ~lm_res(.x, args =args)),
          # res_lm = !!quo(lm(!!quo_updt(args[["Xt"]],
          #                                      txt = as_name(args[["ion2"]]),
          #                                      x = "studE",
          #                                      sepfun = "~")
          #                           )
          #                 ),
          # R2 = summary(res_lm)$r.squared,
          R2 = purrr::map(res_lm, ~(summary(.x)$r.squared)),
          SE_beta = purrr::map(res_lm, ~(unname(summary(.x)$coefficients[2,2]))),
# acf test
          ACF = purrr::map(data, acf_calc)
          # Chi_R2 = R2  * !! quo_updt(my_q = args[["Xt"]],
          #                               txt = as_name(args[["ion2"]]),
          #                               x = "n"
          #                               ),
          # flag_CV = if_else(Chi_R2 >
          #                   qchisq(.95, df = 1),
          #                   "heteroskedasticity",
          #                   "homoskedasticity"
          #                   )
            ) %>%
    tidyr::unnest(cols = c(data,
                           R2,
                           SE_beta#,
                           # Chi_R2,
                           # flag_CV
                           )
                  ) %>%
    mutate(Chi_R2 = R2  * !! quo_updt(my_q = args[["Xt"]],
                                      txt = as_name(args[["ion2"]]),
                                      x = "n"
                                      ),
           flag_CV = if_else(Chi_R2 >
                        qchisq(.95, df = 1),
                        "Ha (heteroskedasticity)",
                        "H0 (homoskedasticity)"
                        )) %>%
    ungroup() %>%
    eval_tidy(expr =  mod_out(output))
}



#' Create stat_R call quosure
#'
#' \code{expr_R} function to generate an quosure that mimics a stat_R call
#' for subsequent usage in dia_R
#'
#' The \code{diag_R} function performs an internal call to stat_R to perform diagnostics
#' on the influence of individial measurements on the blockwise or
#' global (i.e., a complete analysis) stastics. This function provides a
#' convenient way to enter this call as an quosure into the argument
#' \code{args} of \code{diag_R}.
#'
#' @param Xt A character string constituting the ion count rate.
#' @param N A character string constituting the ion counts.
#' @param species A character string constituting the species analysed.
#' @param ion1 A character string constituting the heavy isotope ("13C").
#' @param ion2 A character string constituting the light isotope ("12C").
#'
#' @return A list containing the input arguments as a
#' \code{\link[rlang:quosure]{quosure}}.
#' @export
#' @examples
#' expr_R(Xt = "Xt.pr", N = "N.pr", species = "species.nm", ion1 = "13C",
#'        ion2 = "12C")
expr_R <- function(Xt, N, species, ion1, ion2){

  as_quosures(lst(Xt = parse_expr(Xt),
                  N = parse_expr(N),
                  species = parse_expr(species),
                  ion1 = ion1,
                  ion2 = ion2),
              env = caller_env())

}

#-------------------------------------------------------------------------------
# Not exportet helper functions
#-------------------------------------------------------------------------------

# lm residuals

lm_res <- function(data, args){
  call_lm <- parse_expr(paste0("data$studE~data$", as_name(args[["Xt"]]), ".", as_name(args[["ion2"]])))
  eval_tidy(call2("lm", call_lm))

}

# use the formula i - 0.5/ in, for i = 1,..,n
# this is a vector of the n probabilities ( theoretical cumulative distribution function CDF)
vector_probs <- function(n){
  ((1:unique(n)) - 0.5) / (unique(n))
}

# standard error of quantiles model
hat_QR_se <- function(RQ, TQ, pb, n){
  (sd(RQ) / dnorm(TQ)) * sqrt((pb * (1 - pb))/ unique(n))
}
# confidence interval regression model
hat_Y_se <- function(sigma, hat){
  sigma * sqrt(hat)
}

acf_calc <- function(data){
  acf  <- acf(data$studE, plot = FALSE)
  ci <- qnorm((1 - 0.95) / 2) / sqrt(length(data$studE))
  df.acf <- tibble(lag =  as.vector(acf$lag)[-1],
                   acf =  as.vector(acf$acf)[-1],
                   ci_upper = ci,
                   ci_lower = -ci
                 ) }



# Standard error of regression calculator
#' @export
sigma_calc <- function(res) sqrt(sum((res ^ 2)) / (length(res) - 2))

#' Studentized residuals calculator
#' @export
studE_calc <- function(res, sigma_i, hat_Xi) res / (sigma_i * sqrt(1 - hat_Xi))

#' Cook's D calulator
#' @export
cooksD_calc  <- function(studE, hat_Xi) (studE ^ 2 / 2 ) * (hat_Xi / (1- hat_Xi))

#' Cook's D cutt-off value calculator
#' @export
cooksD_cut_calc  <- function(n) 4 / (n - 2)

#' Jackknife regression coefficient
#' @export
jack_meanR <- function(df, ion1, ion2){

  Xt_ion1 <- enquo(ion1)
  Xt_ion2 <- enquo(ion2)

  Xt_ion1 <- select(df, !! Xt_ion1) %>% pull(!! Xt_ion1)
  Xt_ion2 <- select(df, !! Xt_ion2) %>% pull(!! Xt_ion2)


  R_i <- c((resample::jackknife(Xt_ion1 , mean))$replicates /
           (resample::jackknife(Xt_ion2 , mean))$replicates)

}

#' Jackknife standard error of regression
#' @export
jack_sigma <- function(df, res){

  res <- enquo(res)

  res <- select(df, !! res) %>% pull(!! res)

  sigma_i <- c((resample::jackknife(res, sigma_calc))$replicates)

}

#' @export
sim_R <- function(n = 3000, N_range = 10 ^ 6, reps = 1, ion1, ion2, sys, type, baseR = NULL, offsetR = NULL, seed, ...){

  average_n <- N_range / n
  start_n <- n


  tibble::tibble(simulation = type,
                       n = n,
                       N = as.integer(N_range),
                       R.input = R_gen(start_n,
                                       baseR,
                                       offsetR,
                                       input = "delta",
                                       type = type
                                       ),
                       drift = seq(average_n * (1 - sys), average_n * (1 + sys),
                                   length.out = start_n
                                   ),
                       intercept = average_n
                       ) %>%
      tidyr::expand_grid(., rep = c(1:reps), species = c(ion1, ion2)) %>%
      mutate(simulation = paste(.data$simulation, .data$rep, sep = "-"),
             seed = seed * rep
             ) %>%
# convert common isotope N
      mutate(N = if_else(species == ion2, iso_conv(.data$N, .data$R.input), .data$N)) %>%
# Calculate N of abundant isotope species
      group_by(.data$simulation, .data$species) %>%
      tidyr::nest() %>%
# random variation (Number generation)
      mutate(N.sim= purrr::map(.data$data, ~N_gen(.x, N, n, seed))) %>%
      tidyr::unnest(cols = c(.data$data, .data$N.sim)) %>%
# systematic variation
      mutate(diff = .data$drift - .data$intercept,
             diff = if_else(species == ion2,
                            as.double(iso_conv(.data$diff,
                                               .data$R.input
                                               )
                                      ),
                            .data$diff
                            ),
             N.sim= .data$N.sim + .data$diff,
             Xt.sim = .data$N.sim,
             trend = paste0("linear trend (var: ", sys, ")")
             ) %>%
      ungroup() %>%
      select(-c(drift, intercept, diff))

}

#
# random Poisson ion count generator
N_gen <- function(df, N, n, seed) {


  N <- enquo(N)
  n <- enquo(n)
  seed <- enquo(seed)

  N <- df %>% pull(!! N)
  n <- df %>% pull(!! n)
  seed <- df %>% pull(!! seed)

  set.seed(seed)

  Nsim <- as.double(rpois(n = n, lambda = N / n))

  }

# calculate common isotope count from rare isotope
iso_conv <- function(N, R.sim)  as.integer(N * (1 / R.sim))

R_gen <- function(reps, baseR, offsetR, input = "delta", type) {



  baseR <- calib_R(baseR,
                   standard = "VPDB",
                   type = "composition",
                   input = input,
                   output = "R"
                   )

  offsetR <- calib_R(offsetR,
                    standard = "VPDB",
                    type = "composition",
                    input = input,
                    output = "R"
                    )

  if (type == "ideal") {
    R.sim <- rep(baseR, reps)
    return(R.sim)
    }

  if (type == "constant"){

    R.sim <- approx(c(1, 5 * reps /6, reps),
                    c(baseR, offsetR, offsetR),
                    n = reps ,
                    method = "constant"
                    )$y
    return(R.sim)
    }

  if (type == "gradient") {


    R.sim <- approx(c(1, reps),
                    c(offsetR, baseR),
                    n = reps ,
                    method = "linear"
                    )$y

   return(R.sim)
    }
  }




