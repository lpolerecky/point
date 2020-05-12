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

# Remove zeros
  df <- zeroCt(df,
               !! args[["N"]],
               !! args[["species"]],
               as_name(args[["ion1"]]),
               as_name(args[["ion2"]]),
               !!! gr_by,
               warn = FALSE
               )

# Descriptive an predictive statistics for ion ratios
  tb.aug <- stat_R(df,
                   Xt = !! args[["Xt"]],
                   N =!! args[["N"]],
                   species = !! args[["species"]],
                   ion1 = as_name(args[["ion1"]]),
                   ion2 = as_name(args[["ion2"]]),
                   !!! gr_by,
                   output = "complete",
                   zero = TRUE
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
             "good",
             "bad"
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
                        expr(.data$hat_Xi),
                        expr(.data$studE),
                        expr(.data$CooksD),
                        expr(.data$CooksD_cf),
                        expr(.data$flag_CD),
                        expr(.data$RQ),
                        expr(.data$TQ),
                        expr(.data$hat_RQ),
                        expr(.data$hat_RQ_min),
                        expr(.data$hat_RQ_max),
                        # expr(.data$flag_QQ)
                         ),
           complete = call2( "invisible", expr(.)),
    )
    }

  df %>%
    group_by(!!! gr_by) %>%
    mutate(
# Sum of squared deviation from the mean of x (SSX)
           SSX = sum((!! quo_updt(my_q = args[["Xt"]],
                                  txt = as_name(args[["ion2"]])) -
                      !! quo_updt(my_q = args[["Xt"]],
                                  txt = as_name(args[["ion2"]]),
                                  x = "M")) ^ 2),
# Hat values
           hat_Xi = purrr::pmap_dbl(
             list(Xi = !! quo_updt(my_q = args[["Xt"]],
                                 txt = as_name(args[["ion2"]])),
                  M_X = !! quo_updt(my_q = args[["Xt"]],
                                    txt = as_name(args[["ion2"]]),
                                    x = "M"),
                  SSX = SSX,
                  n = !! quo_updt(my_q = args[["Xt"]],
                                  txt = as_name(args[["ion2"]]),
                                  x = "n")),
                                    hat_calc),
# modelled Y values
           hat_Y = purrr::map2_dbl(!! quo_updt(my_q = args[["Xt"]],
                                               x = "M_R"),
                                   !! quo_updt(my_q = args[["Xt"]],
                                               txt = as_name(args[["ion2"]])),
                                   ~{.x * .y}),
# residuals
           E = purrr::map2_dbl(!! quo_updt(my_q = args[["Xt"]],
                                           txt = as_name(args[["ion1"]])),
                               hat_Y,
                               ~{.x - .y}),
# Standard error of the regression
           sigma = sigma_calc(E),
# Mean Square Error
           MSE = sum(E ^ 2) / !! quo_updt(my_q = args[["Xt"]],
                                          x = "n_R")) %>%
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
           flag_CD = if_else(CooksD > CooksD_cf, "influential", "non-influential")
          ) %>%
# normality test
    mutate(RQ = studE) %>%
    arrange(RQ) %>%
# use the formula i - 0.5/ in, for i = 1,..,n
# this is a vector of the n probabilities ( theoretical cumulative distribution function CDF)
    mutate(prob_vc = vector_probs(!! quo_updt(my_q = args[["Xt"]],
                                              txt = as_name(args[["ion2"]]),
                                              x = "n"
                                              )
                                  ),
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
          hat_RQ_max =  hat_RQ + 2 * hat_RQ_se
          # flag_QQ = purrr::pmap_chr(lst(a = RQ,
          #                               b = hat_RQ_min,
          #                               c = hat_RQ_max
          #                               ),
          #                          ~function(a = a, b = b, c = c){
          #                            if_else(between(a, b, c),
          #                                   "normal",
          #                                   "non-normal"
          #                                   )
          #                           }
          #                          )
          ) %>%
    ungroup() %>%
    eval_tidy(expr =  mod_out(output))
}







#
# mutate(RQ = hat_E_star)  %>%
#   # use the formula i - 0.5/ in, for i = 1,..,n
#   # this is a vector of the n probabilities ( theoretical cumulative distribution function CDF)
#   mutate(probs = ((1:unique(!!n.Xt2)) - 0.5) / (unique(!!n.Xt2)),
#          # calculate normal (Theoretical) quantiles using mean and standard deviation from
#          TQ = qnorm(probs, mean(RQ), sd(RQ)),
#          # the standard error is calculated with,
#          hat_RQ = mean(RQ) + sd(RQ) * TQ,
#          hat_RQ_se = (sd(RQ) / dnorm(TQ)) * sqrt((probs * (1 - probs))/ unique(!!n.Xt2)),
#          hat_RQ_min =  hat_RQ - 2 * hat_RQ_se,
#          hat_RQ_max =  hat_RQ + 2 * hat_RQ_se)  %>%
#   rowwise() %>%
#   mutate(RQ_lab = if_else(between(RQ, hat_RQ_min, hat_RQ_max),
#                           "normal", "non-normal")) %>%
#   ungroup()


# use the formula i - 0.5/ in, for i = 1,..,n
# this is a vector of the n probabilities ( theoretical cumulative distribution function CDF)
vector_probs <- function(n){
  ((1:unique(n)) - 0.5) / (unique(n))
}

# standard error of quantiles model
hat_QR_se <- function(RQ, TQ, pb, n){
  (sd(RQ) / dnorm(TQ)) * sqrt((pb * (1 - pb))/ unique(n))
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

# Hat value calculcator
hat_calc <- function(Xi, M_X, SSX, n) ((Xi - M_X)^2 / SSX) + (1 / n)

# Standard error of regression calculator
#' @export
sigma_calc <- function(res) sqrt(sum(res ^ 2) / (length(res) - 2))

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

