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
#'                  method = "standard",
#'                  args = expr_R(Xt = "Xt.pr",
#'                                N = "N.pr",
#'                                species = "species.nm",
#'                                ion1 = "13C",
#'                                ion2 = "12C"),
#'                  file.nm,
#'                  bl.mt)
diag_R <- function(df, method = "standard", args = expr_R(NULL), ...){

  gr_by <- enquos(...)

# ID for connecting flag to original dataframe
  df <- ID_builder(df, !! args[["species"]], !!! gr_by)

# Descriptive an predictive statistics for 13C/12C ratios
  tb.aug <- stat_R(df,
                     !! args[["Xt"]],
                     !! args[["N"]],
                     species = !! args[["species"]],
                     ion1 = as_name(args[["ion1"]]),
                     ion2 = as_name(args[["ion2"]]),
                     !!! gr_by,
                     output = "complete") %>%
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
                                    between(!! quo_updt(my_q = args[["Xt"]],
                                                        x = "R"),
                                            unique(.data$lower),  # mean - 2SD
                                            unique(.data$upper)), # mean + 2SD
                                     "good",
                                     "bad"
                                    )
                    ) %>%
              ungroup() %>%
              select(.data$ID, .data$flag)

# Cameca style augmented datafile
  tb.aug <- left_join(df, tb.aug, by = "ID")

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
