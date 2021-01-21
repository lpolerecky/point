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
#' @param .df A tibble containing processed ion count data.
#' @param .ion1 A character string constituting the heavy isotope ("13C").
#' @param .ion2 A character string constituting the light isotope ("12C").
#' @param ... Variables for grouping.
#' @param .method Character string for the method for diagnostics (default =
#' \code{"CooksD"}, see details).
#' @param .reps Numeric setting the number of repeated iterations of outlier
#' detection (default = 1).
#' @param .Xt A variable constituting the ion count rate (defaults to
#' variables generated with \code{read_IC()})
#' @param .N A variable constituting the ion counts (defaults to variables
#' generated with \code{read_IC()}.).
#' @param .output Character string determining whether the returned values in a
#' minimal version `"flag"` (original dataset + diagnostics) or an extended
#' version with all the intermediate steps of ion- and isotope-wise summary
#' statistics `"complete"`.
#' @param .hyp Hypothesis test. Only usable in combination with a selection of
#' methods (see above argument \code{.method}) and details.
#' @param .alpha_level. The significance level of the hypothesis test and
#' rejection level for outliers.
#'
#' @return A \code{\link[tibble:tibble]{tibble}} containing either the original
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
#' tb_R <- stat_R(tb_pr, "13C", "12C", file.nm, sample.nm, .output = "complete",
#'                .zero = TRUE)
#'
#' # CAMECA style augmentation of ion count data for isotope ratios
#' tb_dia <- Cameca(tb_R,"13C", "12C", file.nm, .output = "flag")
Cameca <- function(.df, .ion1, .ion2, ..., .Xt = Xt.pr, .t = t.nm,
                   .output = "complete", .hyp = "none", .alpha_level = 0.05){

  # Grouping
  gr_by <- enquos(...)

  # Quoting the call (user-supplied expressions)
  Xt <- enquo(.Xt)
  t <- enquo(.t)

  # Heavy isotope
  Xt1 <- quo_updt(Xt, post = .ion1) # count rate
  # Light isotope
  Xt2 <- quo_updt(Xt, post = .ion2) # count rate

  # R
  R <- quo_updt(Xt, pre = "R")
  M_R <- quo_updt(Xt, pre = "M_R")

  # Fitted isotope R (+ heavy isotope) and variance (sigma)
  R <- quo_updt(Xt, pre = "R")
  hat_R <- quo_updt(R, pre = "hat", post = .ion1)
  hat_s_R <- quo_updt(R, pre = "hat_s", post = .ion1)
  hat_Xt1 <- quo_updt(Xt, pre = "hat", post = .ion1)
  hat_s_Xt1 <- quo_updt(Xt, pre = "hat_s", post = .ion1)

  group_by(.df, !!! gr_by) %>%
    mutate(
      !! hat_s_R := sd(!!R),
      !! hat_R := !! M_R,
      !! hat_Xt1 := !! hat_R * !!Xt2,
      !! hat_s_Xt1 := !! hat_s_R * !!Xt2,
      flag = if_else(
        between(
          !! R,
          !! hat_R + qnorm((.alpha_level / 2)) * !!hat_s_R, # mean - 2SD
          !! hat_R + qnorm(1 - (.alpha_level / 2)) * !!hat_s_R  # mean + 2SD
          ),
        "confluent",
        "divergent"
        ),
      flag = as.factor(flag)
      ) %>%
    ungroup()

}
