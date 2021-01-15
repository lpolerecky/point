#' Family of diagnostics functions for isotope count ratios
#'
#' \code{Cameca} default CAMECA diagnostics
#' \code{CooksD} regression diagnostics based on Cook's D
#'
#' These functions perform a specific set of diagnostics to term anomalous
#' values in raw ion count data of an isotope pair. The wrapper function
#' \code{Diag_R} is more convenient as it defines all the ion- and isotope-wise
#' statistics required for the diagnostics.
#'
#' @param .df A tibble with ion count data and statistics for ion- and
#' isotope-wise statistics.
#' @param .args A list of quosures pertaining to the variables required for a
#' call to stat_R. The function expr_R can be used to simplify setting this
#' argument.
#' @param ... Variables for grouping.
#' @param .output Character string determining whether the returned values in a
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
#' tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # Processing raw ion count data
#' tb_pr <- cor_IC(tb_rw)
#'
#' # Descriptive an predictive statistics for 13C/12C ratios (note .output
#' # argument and remove zero count analysis)
#' tb_R <- stat_R(tb_pr, "13C", "12C", file.nm, sample.nm, .output = "complete", .zero = TRUE)
#'
#' # Mimic stat_R call
#' args <-  expr_R("Xt.pr", "N.pr", "species.nm", "13C", "12C")
#'
#' # CAMECA style augmentation of ion count data for isotope ratios
#' tb_dia <- Cameca(tb_R, file.nm, sample.nm, .args = args, .output = "flag")
Cameca <- function(.df, ..., .args = expr_R(NULL), .output){

  # Grouping
  gr_by <- enquos(...)

  # New quosures
  # Single ion (light)
  Xt2 <- quo_updt(.args[["Xt"]], post = as_name(.args[["ion2"]]))

  # R
  R <- quo_updt(.args[["Xt"]], pre = "R")
  M_R <- quo_updt(.args[["Xt"]], pre = "M_R")
  S_R <- quo_updt(.args[["Xt"]], pre = "S_R")

  group_by(.df, !!! gr_by) %>%
    mutate(
      lower = !! M_R - 2 * !! S_R ,
      upper = !! M_R + 2 * !! S_R,
      hat_Y = !! M_R * !! Xt2,
      flag = if_else(
        between(!! R,
        unique(.data$lower),  # mean - 2SD
        unique(.data$upper)), # mean + 2SD
        "good",
        "bad"
        ),
      flag = as.factor(flag)
      ) %>%
    ungroup()



}
