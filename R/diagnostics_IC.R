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
#' @param reps The number of iterations of diagnostic-augmentation cycles.
#' @param ... Variables for grouping.
#' @param plot Logical indicating whether plot is produced, if FALSE a list of
#' tibbles with length `reps` is produced.
#' @param Character string determing whether the returned a \code{"static"} plot
#' or an \code{"interactive"} plot with \code{\link[plotly:ggplotly]{ggplotly}}.
#'
#' @return A list containing several \code{\link[tibble:tibble]{tibble}} the
#' elements named df contain the original dataset and augmented versions from
#' element number 2 onward. The elements named results contain the flag
#' variable for outlier detection as well as associated statistics of the
#' selected procedure.
#' @export
#' @examples
#' # Modelled ion count dataset
#' # Create expression for mimicking R_stat() call
#' expr_R_stat <- expr_R(Xt = "Xt.sim",
#'                       N = "N.sim",
#'                       species = "species",
#'                       ion1 = "13C",
#'                       ion2 = "12C"
#'                       )
#'
#' # Cook's D style diagnostic-augmentation of ion count data for
#' # isotope ratios; 3 repeats
#' diag_R(sim_IC_extremes,
#'        method = "CooksD",
#'        args = expr_R_stat,
#'        reps = 3,
#'        simulation,
#'        trend
#'        )
diag_R <- function(df,
                   method = "Cameca",
                   args = expr_R(NULL),
                   reps = 1,
                   ...,
                   plot = TRUE,
                   plot_type = "interactive"
                   ){

  gr_by <- enquos(...)

  max <- reps + 1

# empty list for iteration storage
  ls.tb <- rlang::rep_named(as.character(1:max), lst())

# Remove zeros
  if (method != "Cameca"){

    df <- zeroCt(df,
                 !! args[["N"]],
                 !! args[["species"]],
                 as_name(args[["ion1"]]),
                 as_name(args[["ion2"]]),
                 !!! gr_by,
                 warn = FALSE
                 )
  }

# ID for connecting flag to original dataframe
  df <- ID_builder(df, !! args[["species"]], !!! gr_by)
# set initial dataset
  df <- filter(df,
               !! args[["species"]] == !! args[["ion1"]] |
               !! args[["species"]] == !! args[["ion2"]]
               )
  ls.tb[[1]] <- lst(df = df, results = NULL)

# Execute repeated cycles of augmentation
  ls.tb <- ls.tb  %>%
    purrr::accumulate(rerun_diag_R,
                      method = method,
                      args = args,
                      !!! gr_by
                      )
  if (plot) {

    ls.tb %>%
      plot_diag_R(args = args, !!! gr_by, plot_title = method, type = plot_type)

    } else {

      return(ls.tb)

      }

}

#' @export
rerun_diag_R <- function(out,
                         input,
                         method = "Cameca",
                         args = expr_R(NULL),
                         ...
                         ){

  gr_by <- enquos(...)

  # variables to be saved
  results_vars <- set_names(colnames(out$df))
  deselect_vars <- results_vars  %in% append(sapply(gr_by, as_name), "ID")

  variables <- parse_exprs(results_vars[!deselect_vars])

  # if called this way then output is set fixed (more fleible use with diag_R_exec)
  out <- diag_R_exec(out$df,
                     method = method,
                     args = args,
                     !!! gr_by,
                     output = "flag"
                     )

  # save augmented dataframe for next cycle
  df.aug <- filter(out, flag == "good") %>%
    select(.data$ID, !!!gr_by, !!!variables)

  # save results; flag ad statistics
  results <- select(out, -c(!!!variables))


  out <- lst(df = df.aug, results = results)

  return(out)

}

#' @export
diag_R_exec <- function(df,
                        method = "Cameca",
                        args = expr_R(NULL),
                        ...,
                        output
                        ){

  gr_by <- enquos(...)

# Method selection
  diag_method <- function(method){
    switch(method,
           Cameca = call2("Cameca_R", expr(.), args, !!! gr_by, output = output),
           CooksD = call2("CooksD_R", expr(.), args, !!! gr_by, output = output),
           norm_E = call2("norm_E", expr(.), args, !!! gr_by, output = output),
           Rm = call2("Rm", expr(.), args, !!! gr_by, output = output),
           QQ = call2("QQ", expr(.), args, !!! gr_by, output = output)

    )
    }

  # if (method == "CooksD" | method == "Cameca"){
# Descriptive an predictive statistics for ion ratios
    tb.aug <- stat_R(df,
                     Xt = !! args[["Xt"]],
                     N =!! args[["N"]],
                     species = !! args[["species"]],
                     ion1 = as_name(args[["ion1"]]),
                     ion2 = as_name(args[["ion2"]]),
                     !!! gr_by,
                     output = "complete"
                     ) %>%
               eval_tidy(expr = diag_method(method))
  # }


# Datafile with flag values associated to diagnostics
  if (output == "flag") return(left_join(df, tb.aug, by = "ID"))
  if (output == "complete") return(tb.aug)

  }

#' Evaluate effect size of diagnostics
#'
#' \code{eval_diag} function for the evaluation of effect size of diagnostics
#'
#' @param ls_df A list of tibbles as generated by \code{diag_R}.
#' @param args A list of quosures pertaining to the variables required for a
#' call to stat_R. The function expr_R can be used to simplify setting this
#' argument.
#' @param df2 A dataframe with the ionization tend of the major ion in per
#' mille.
#' @param RSXt A variable constituting the relative standard deviation of the
#' major ion in the dataframe enetered for the argument \code{df2}.
#' @param n A variable indicating the number of observation in the list of
#' dataframes entered for the argument \code{ls_df}
#' @param ... Variables for grouping.
#' @param output A character string for output as summary statistics ("sum") and
#' statistics with the original data ("complete").
#'
#' @export
#' @examples
#'
#' # estimate strength of ionization trend as variation in major ion in per mille
#' tb.var12C <- predict_var(tb.pr, Xt.pr, N.pr, t.rw, species.nm, "12C", file.nm)
#'
#' # evaluation of diagnostics
#' tb.eval <- eval_diag(ls.diag,
#'                      args = expr_R_stat,
#'                      tb.var12C,
#'                      RS_Xt.pr.12C,
#'                      n_R_Xt.pr,
#'                      file.nm,
#'                      bl.mt,
#'                      output = "sum"
#'                      )
eval_diag <- function(ls_df,
                      args,
                      df2,
                      RSXt,
                      n,
                      ...,
                      output = "sum"
                      ){

  RSR <- quo_updt(args[["Xt"]], x = "RS_R")
  RSXt <- enquo(RSXt)
  n <- enquo(n)
  gr_by <- enquos(...)

# connecting the
  gr_vc <- sapply(gr_by, as_name)
# connecting the single ion variable
  uni_gr <- gr_vc[gr_vc %in% colnames(df2)]

  df  <- reduce_diag(ls_df,
                     type = "df",
                     args = args,
                     !!! gr_by
                     )

  ls.trans <- lst(
    quo(!! n),
    quo(c(diff(!! n, 1), NA)),
    quo(c(diff(!! RSR, 1), NA)),
    quo(diff_RS ^ 2 * diff_n),
    quo((!! RSR) ^ 2 * !! n),
    quo(SSD / SSA),
    quo(n2 * nz)
    )

# Names of parameters for diagnostic evaluation
  ls.names <- c("nz", "diff_n", "diff_RS", "SSD", "SSA", "n2", "chi_n2")

# Set names
  ls.trans <- set_names(ls.trans, nm = ls.names)


# Switch output complete dataset, stats or summary stats
  mod_cal <- function(type) {
    switch(type,
           complete = call2( "mutate", expr(.), !!! ls.trans),
           sum = call2( "transmute",  expr(.), !!! ls.trans)
    )
    }

  ls.preserve <-lst(
    quo(sum(n2, na.rm = TRUE)),
    quo(sum(chi_n2, na.rm = TRUE))
    )

# Names to be preserved
  ls.names2 <- c("n2", "chi_n2")

# Set names
  ls.preserve <- set_names(ls.preserve, nm = ls.names2)

# Switch to determine what extend of the data is preserved
  filter_cal <- function(type) {
    switch(type,
           complete = call2("select",
                            expr(.),
                            expr(-c(ls.names[!ls.names %in% ls.names2]))
                            ),
           sum = call2("summarise",
                       expr(.),
                       !!! ls.preserve
                       )
           )
    }


# Evaluate expressions and calls
    df_eval <- df %>%
      group_by(!!! gr_by) %>%
      arrange(desc(execution)) %>%
      eval_tidy(expr = mod_cal(output)) %>%
      eval_tidy(expr = filter_cal(output))

    if (output == "sum") {
      df <- lst(filter(df, execution == 1), df_eval, df2) %>%
        purrr::reduce2(lst(gr_vc, uni_gr), left_join) %>%
        mutate(crit_val = purrr::map_dbl(!! RSXt, crit_size),
               flag = if_else(.data$chi_n2 > crit_val, "variable", "non-variable")
               )
      return(df)
    }

    if (output == "complete") {
      df <- left_join(df_eval, df2, by = uni_gr) %>%
        mutate(crit_val = purrr::map_dbl(!! RSXt, crit_size),
               flag = if_else(.data$chi_n2 > crit_val, "variable", "non-variable")
               )
      return(df)
    }


  }

#' Reduce diagnostis
#'
#' @export
reduce_diag <- function(ls, type = "df", args = expr_R(NULL), ..., output = "sum"){

  gr_by <- quos(...)

  type_reduction <- function(type){
    switch(type,
           results = call2("bind_rows", expr(.), .id = "execution"),
           df = call2("map_dfr",
                      expr(.),
                      expr(~stat_R(df = .x,
                                   Xt = !! args[["Xt"]],
                                   N = !! args[["N"]],
                                   species = !! args[["species"]],
                                   ion1 = as_name(args[["ion1"]]),
                                   ion2 = as_name(args[["ion2"]]),
                                   !!! gr_by,
                                   output = output
                      )
                      ),
                      .id = "execution",
                      .ns = "purrr"
           )
    )
  }

  # Reduce the results to a single dataframe
  ls %>%
    purrr::transpose() %>%
    purrr::pluck(type) %>%
    eval_tidy(expr = type_reduction(type))

}


# CooksD_R <- function(df, args = expr_R(NULL), ..., output){
#
#   gr_by <- enquos(...)
#
# # Switch output complete dataset, stats or summary stats
#   mod_out <- function(output) {
#     switch(output,
#            flag = call2("select",
#                         expr(.),
#                         expr(.data$ID),
#                         expr(.data$E),
#                         expr(.data$sigma),
#                         expr(.data$hat_Y),
#                         expr(.data$hat_R2),
#                         expr(.data$hat_Y_min),
#                         expr(.data$hat_Y_max),
#                         expr(.data$hat_Xi),
#                         expr(.data$studE),
#                         expr(.data$CooksD),
#                         expr(.data$CooksD_cf),
#                         expr(.data$flag),
#                         expr(.data$RQ),
#                         expr(.data$TQ),
#                         expr(.data$hat_RQ),
#                         expr(.data$hat_RQ_min),
#                         expr(.data$hat_RQ_max),
#                         expr(.data$flag_QQ),
#                         expr(.data$flag_CM),
#                         expr(.data$R2),
#                         expr(.data$SE_beta),
#                         expr(.data$Chi_R2),
#                         expr(.data$flag_CV),
#                         expr(.data$ACF),
#                         expr(.data$flag_IR)
#                          ),
#            complete = call2( "invisible", expr(.)),
#     )
#   }
#
# # update ion names to match stat_R output in case of space seperation
#   args[["ion1"]] <- str_replace_all(as_name(args[["ion1"]]), " ", "")
#   args[["ion2"]] <- str_replace_all(as_name(args[["ion2"]]), " ", "")
#
#
#   df %>%
#     group_by(!!! gr_by) %>%
#     mutate(
# # Hat values
#            hat_Xi = stats::hat(!! quo_updt(my_q = args[["Xt"]],
#                                             txt = as_name(args[["ion2"]])
#                                            )
#                                ),
# # modelled Y values
#            hat_Y = !! quo_updt(my_q = args[["Xt"]],
#                                x = "M_R"
#                                ) *
#                    !! quo_updt(my_q = args[["Xt"]],
#                                txt = as_name(args[["ion2"]])
#                                ),
#            # hat_R2 = (sd(hat_Y) ^ 2 *
#            #           !! quo_updt(my_q = args[["Xt"]],
#            #                         x = "n_R"
#            #                       )
#            #           ) / ((!! quo_updt(my_q = args[["Xt"]],
#            #                             txt = as_name(args[["ion2"]]),
#            #                             x = "S"
#            #                             )
#            #                  ) ^ 2 * !! quo_updt(my_q = args[["Xt"]],
#            #                                    x = "n_R"
#            #                                     )
#            #                ),
#            ESS = (sd(hat_Y) ^ 2 * !! quo_updt(my_q = args[["Xt"]],
#                                               x = "n_R"
#                                               )
#                   ),
#            # hat_Chi2 =  hat_R2 * !! quo_updt(my_q = args[["Xt"]],
#            #                                  x = "n_R"
#            #                                  ),
#
# # residuals
#            E = !! quo_updt(my_q = args[["Xt"]],
#                            txt = as_name(args[["ion1"]])
#                            ) - hat_Y,
# # Standard error of the regression
#            sigma = sigma_calc(E),
# # 95 CI of the regression
#            hat_Y_min = hat_Y - 2 * hat_Y_se(sigma, hat_Xi),
#            hat_Y_max = hat_Y + 2 * hat_Y_se(sigma, hat_Xi),
# # sum of Square residuals
#            SSR = sum(E ^ 2),
# # coeffecient of determination
#            hat_R2 =  1 - (SSR/(ESS + SSR))
#
#             ) %>%
#     tidyr::nest() %>%
#     mutate(
# # jackknifed mean R
#            M_Ri = purrr::map(data,
#                              ~jack_meanR(df = .x,
#                                          ion1 = !! quo_updt(my_q = args[["Xt"]],
#                                                             txt = as_name(args[["ion1"]])),
#                                          ion2 = !! quo_updt(my_q = args[["Xt"]],
#                                                             txt = as_name(args[["ion2"]]))
#                                         )
#                             )
#
#            ) %>%
#     tidyr::unnest(cols = c(data, M_Ri)) %>%
#     mutate(
# # jackknifed modelled Y values
#            hat_Yi = purrr::map2_dbl(M_Ri,
#                                     !! quo_updt(my_q = args[["Xt"]],
#                                                 txt = as_name(args[["ion2"]])),
#                                     ~{.x * .y}),
# # residuals with i-th value removed
#            Ei = purrr::map2_dbl(!! quo_updt(my_q = args[["Xt"]],
#                                             txt = as_name(args[["ion1"]])),
#                                 hat_Yi,
#                                 ~{.x - .y})
#           ) %>%
#     tidyr::nest() %>%
#     mutate(
# # standard error of regression (external with i-th residual removed)
#           sigma_i = purrr::map(data, ~jack_sigma(.x, Ei))
#           ) %>%
#     tidyr::unnest(cols = c(data, sigma_i)) %>%
#     mutate(
#            studE = purrr::pmap_dbl(lst(res = E,
#                                        sigma_i = sigma_i,
#                                        hat_Xi = hat_Xi
#                                        ),
#                                    studE_calc
#                                    ),
#            diag_R(tb.pr,
#                   method = "normE",
#                   args = expr_R_stat,
#                   reps = 3,
#                   file.nm
#            )
#            ) %>%
# # normality test
#     # mutate(RQ = studE) %>%
#     # arrange(RQ) %>%
# # use the formula i - 0.5/ in, for i = 1,..,n for probs
# # this is a vector of the n probabilities ( theoretical cumulative distribution function CDF)
#    mutate(prob_vc =  vector_probs(!! quo_updt(my_q = args[["Xt"]],
#                                               x = "n_R"
#                                               )
#                                    ),
#           RQ = unname(quantile(studE, probs = prob_vc)),
# # calculate normal (Theoretical) quantiles using mean and standard deviation from
#           TQ = qnorm(prob_vc, mean(RQ), sd(RQ)),
# # the standard error is calculated with,
#           hat_RQ = mean(RQ) + sd(RQ) * TQ,
#           hat_RQ_se = hat_QR_se(RQ,
#                                 TQ,
#                                 prob_vc,
#                                 !! quo_updt(my_q = args[["Xt"]],
#                                             txt = as_name(args[["ion2"]]),
#                                             x = "n"
#                                             )
#                                 ),
#           hat_RQ_min =  hat_RQ - 2 * hat_RQ_se,
#           hat_RQ_max =  hat_RQ + 2 * hat_RQ_se,
#           flag_QQ = if_else(
#             (nortest::ad.test(studE))$p.value < 0.05,
#             "Ha (non-normal)",
#             "H0 (normal)"
#             ),
# # t-test flag for mu0 (aka the conditional mean of epsilon) being zero
#           flag_CM = if_else((t.test(studE, mu = 0))$p.value  < 0.05,
#                             "Ha (mu0 is not zero)",
#                             "H0 (mu0 is zero)"
#                             )
#           ) %>% #,
# # hetroscadasticity test (Breusch Pagan test)(level of confidence 95%;
# # cut-off 0.05 for H0 rejection),
#     tidyr::nest() %>%
#     mutate(
#           res_lm = purrr::map(data, ~lm_res(.x, args = args)),
#           # res_lm = !!quo(lm(!!quo_updt(args[["Xt"]],
#           #                                      txt = as_name(args[["ion2"]]),
#           #                                      x = "studE",
#           #                                      sepfun = "~")
#           #                           )
#           #                 ),
#           # R2 = summary(res_lm)$r.squared,
#           R2 = purrr::map(res_lm, ~(summary(.x)$r.squared)),
#           SE_beta = purrr::map(res_lm, ~(unname(summary(.x)$coefficients[2,2]))),
# # acf test
#           ACF = purrr::map(data, acf_calc)
#
#
#           # Chi_R2 = R2  * !! quo_updt(my_q = args[["Xt"]],
#           #                               txt = as_name(args[["ion2"]]),
#           #                               x = "n"
#           #                               ),
#           # flag_CV = if_else(Chi_R2 >
#           #                   qchisq(.95, df = 1),
#           #                   "heteroskedasticity",
#           #                   "homoskedasticity"
#           #                   )
#             ) %>%
#     tidyr::unnest(cols = c(data,
#                            R2,
#                            SE_beta#,
#                            # Chi_R2,
#                            # flag_CV
#                            )
#                   ) %>%
#     mutate(Chi_R2 = R2  * !! quo_updt(my_q = args[["Xt"]],
#                                       txt = as_name(args[["ion2"]]),
#                                       x = "n"
#                                       ),
#            flag_CV = if_else(Chi_R2 >
#                         qchisq(.95, df = 1),
#                         "Ha (heteroskedasticity)",
#                         "H0 (homoskedasticity)"
#                         ),
#            LB_test = stats::Box.test(studE, type = "Ljung-Box")$p.value,
#            flag_IR = if_else(LB_test < 0.05,
#                              "Ha (dependence of residuals)",
#                              "H0 (independence of residuals)"
#                              )
#            ) %>%
#     ungroup() %>%
#     eval_tidy(expr =  mod_out(output))
# }


#' Create stat_R call quosure
#'
#' \code{expr_R} and \code{expr_cor} functions generate an quosure that mimics a
#' stat_R call or cor_IC call for subsequent usage in dia_R
#'
#' The \code{diag_R} function performs an internal call to stat_R or cor_IC call
#' to perform diagnostics on the influence of individial measurements on the
#' blockwise or global (i.e., a complete analysis) stastics. This function
#' provides a convenient way to enter this call as an quosure into the argument
#' \code{args} of \code{diag_R}.
#'
#' @param Xt A character string constituting the ion count rate.
#' @param N A character string constituting the ion counts.
#' @param species A character string constituting the species analysed.
#' @param ion1 A character string constituting the heavy isotope ("13C").
#' @param ion2 A character string constituting the light isotope ("12C").
#' @param t A character string constituting the time increment of measurement.
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
                  ion2 = ion2
                  ),
              env = caller_env()
              )

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
hat_Y_se <- function(sigma, hat_Xi){
  sigma * sqrt(hat_Xi)
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
cooksD_calc  <- function(studE, hat_Xi) (studE ^ 2 / 2 ) * (hat_Xi / (1 - hat_Xi))

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

crit_size <- function(RS_Xt.ion2){

  null_dist <- null_dist %>%
    mutate(min.val = abs(RS_Xt.ion2 - parse_number(trend))) %>%
    filter(min.val == min(.data$min.val)) %>%
    pull(y_star)

}

