#' Diagnostics raw ion count data
#'
#' \code{diag_R} function for propagation of uncertainty for single ions
#'
#' Ion count data consisting of time-incremented integer values are
#'
#' @param df A tibble containing processed ion count data
#' @param method Character string for the type of diagnostics
#' @param args A list of expression pertaining to the variables required for a
#' call to stat_R. The function expr_R can be used to simplify this.
#' @param ... Variables for grouping
#'
#' @return A t\code{\link[tibble:tibble]{tibble}} containing the original
#' dataset and new columns related to the diagnostics. Flag variables enable
#' convenient filtering of the tibble for augmentation of the dataset.
#'
#' @examples
#' # Use point_example() to access the examples bundled with this package in the
#' # inst/extdata directory.
#'
#' # raw data containing 13C and 12C counts on carbonate
#' tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # processing raw ion count data
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
#'                  file.nm, bl.mt)
#'
#' @export
diag_R <- function(df, method = "standard", args = expr_R(NULL), ...){

  gr_by <- enquos(...)

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
                                            unique(lower),  # mean - 2SD
                                            unique(upper)), # mean + 2SD
                                     "good",
                                     "bad"
                                    )
                    ) %>%
              ungroup() %>%
              select(ID, flag)

# Cameca style augmented datafile
  tb.aug <- left_join(df, tb.aug, by = "ID")

}


# Function to create expr for the variables for stat_R call for subsequent usage in dia_R
expr_R <- function(Xt, N, species, ion1, ion2){

  as_quosures(lst(Xt = rlang::parse_expr(Xt),
                  N = rlang::parse_expr(N),
                  species = rlang::parse_expr(species),
                  ion1 = ion1,
                  ion2 = ion2),
              env = rlang::caller_env())

  }



