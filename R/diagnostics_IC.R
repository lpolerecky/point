#' Diagnostics raw ion count data
#'
#' \code{diag_R} function for propagation of uncertainty for single ions
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
#' @param method Character string for the type of diagnostics. Currently only
#' "standard" is supported, which pertains to the default Cameca software
#' setting.
#' @param args A list of quosures pertaining to the variables required for a
#' call to stat_R. The function expr_R can be used to simplify setting this
#' argument.
#' @param ... Variables for grouping.
#'
#' @return A t\code{\link[tibble:tibble]{tibble}} containing the original
#' dataset and new columns related to the diagnostics. The flag variable enables
#' convenient filtering of the original tibble for an augmentation of the
#' original dataset.
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
           CooksD = call2("CooksD_R", expr(.), args, !!!gr_by, output = output))
  }

# Remove zeros
  df <- zeroCt(df,
               !! args[["N"]],
               !! args[["species"]],
               as_name(args[["ion1"]]),
               as_name(args[["ion2"]]),
               !!! gr_by,
               warn = FALSE)

# Descriptive an predictive statistics for ion ratios
  tb.aug <- stat_R(df,
                   Xt = !! args[["Xt"]],
                   N =!! args[["N"]],
                   species = !! args[["species"]],
                   ion1 = as_name(args[["ion1"]]),
                   ion2 = as_name(args[["ion2"]]),
                   !!! gr_by,
                   output = "complete",
                   zero = TRUE) %>%
              eval_tidy(expr = diag_method(method))

# Datafile with flag values associated to diagnostics
  if (output == "flag"){

    return(left_join(df, tb.aug, by = "ID") %>% select(-.data$ID))

    }

  if (output == "complete"){

      return(tb.aug)
    }
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


#' Cook's D
#' @export
CooksD_R <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

# Switch output complete dataset, stats or summary stats
  mod_out <- function(output) {
    switch(output,
           flag = call2( "select", expr(.), expr(.data$ID),
                                            expr(.data$E),
                                            expr(.data$hat_Y),
                                            expr(.data$hat_Xi),
                                            expr(.data$studE),
                                            expr(.data$CooksD),
                                            expr(.data$CooksD_cf),
                                            expr(.data$CooksD_lab)),
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
                                       hat_Xi = hat_Xi),
                                   studE_calc),
           CooksD = purrr::map2_dbl(studE, hat_Xi, cookD_calc),
           CooksD_cf = purrr::map_dbl(!! quo_updt(my_q = args[["Xt"]],
                                                  x = "n_R"),
                                      cookD_cut_calc),
           CooksD_lab = if_else(CooksD > CooksD_cf, "influential", "non-influential")
          ) %>%
    ungroup() %>%
    eval_tidy(expr =  mod_out(output))
}

#' Hat value calculcator
#' @export
hat_calc <- function(Xi, M_X, SSX, n) ((Xi - M_X)^2 / SSX) + (1 / n)

#' Standard error of regression calculator
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

  sigma_i <- c((resample::jackknife(res, sigma_calc ))$replicates)

}

#' Cameca style diagnostics
#' @export
Cameca_R <- function(df, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

# Switch output complete dataset, stats or summary stats
  mod_out <- function(output) {
    switch(output,
           flag = call2( "select", expr(.), expr(.data$ID),
                                            expr(.data$lower),
                                            expr(.data$upper),
                                            expr(.data$flag)),
           complete = call2( "invisible", expr(.)),
    )
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
                         between(!! quo_updt(my_q = args[["Xt"]],x = "R"),
                         unique(.data$lower),  # mean - 2SD
                         unique(.data$upper)), # mean + 2SD
                         "good",
                         "bad")
          ) %>%
    ungroup() %>%
    eval_tidy(expr =  mod_out(output))
    # select(.data$ID, .data$lower, .data$upper, .data$flag)
 }
