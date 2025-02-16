#' Family of diagnostics functions for isotope count ratios
#'
#' \code{Cameca} default CAMECA diagnostics
#' \code{CooksD} regression diagnostics based on Cook's D
#'
#' These functions perform a specific set of diagnostics to term anomalous
#' values in raw ion count data of an isotope pair. The wrapper function
#' \code{diag_R} is more convenient as it defines all the ion- and isotope-wise
#' statistics required for the diagnostics.
#'
#' @param .IC A tibble containing processed ion count data.
#' @param .ion1 A character string constituting the heavy isotope ("13C").
#' @param .ion2 A character string constituting the light isotope ("12C").
#' @param ... Variables for grouping.
#' @param .X A variable constituting the ion count rate (defaults to
#' variables generated with \code{read_IC()})
#' @param .N A variable constituting the ion counts (defaults to variables
#' generated with \code{read_IC()}.).
#' @param .species A variable constituting the species analysed (defaults to
#' variables generated with \code{read_IC()}.).
#' @param .t A variable constituting the time of the analyses (defaults to
#' variables generated with \code{read_IC()}.).
#' @param .output Character string determining whether the returned values in a
#' minimal version `"flag"` (original dataset + diagnostics) or an extended
#' version with all the intermediate steps of ion- and isotope-wise summary
#' statistics `"complete"`.
#' @param .alpha_level The significance level of the hypothesis test and
#' rejection level for outliers.
#' @param .hyp Hypothesis test appropriate for the selected method.
#' @param .mc_cores Number of workers for parallel execution (Does not work on
#'  Windows).
#'
#' @return A \code{tibble\link[tibble:tibble]{tibble}()} containing either the
#' original dataset with new columns related to the diagnostics or only the
#' diagnostics. The flag variable enables convenient filtering of the original
#' tibble for an augmentation of the original dataset.
#' @export
#' @examples
#' # Descriptive an predictive statistics for 13C/12C ratios (note .output
#' # argument and remove zero count analysis)
#' tb_R <- stat_R(real_IC, "13C", "12C", file.nm, sample.nm,
#'                .output = "complete", .zero = TRUE)
#'
#' # CAMECA style augmentation of ion count data for isotope ratios
#' Cameca(tb_R, "13C", "12C", file.nm, .X = Xt.pr, .N = N.pr,
#'        .species = species.nm, .t = t.nm, .output = "flag")
Cameca <- function(.IC, .ion1, .ion2, ..., .X = NULL, .N = NULL, .species = NULL,
                   .t = NULL, .output = "complete", .alpha_level = 0.05,
                   .hyp = "none", .mc_cores = 1) {

  # Grouping
  gr_by <- enquos(...)

  # Quoting the call (user-supplied expressions)
  args <- enquos(.X = .X, .N = .N, .species = .species, .t = .t)

  # R quosures
  args <- rlang::list2(
    !!! args,
    X1 = quo_updt(args[[".X"]], post = .ion1), # count rate rare isotope
    X2 = quo_updt(args[[".X"]], post = .ion2), # count rate common isotope
    !!! arg_builder(args, "R"),
    R  = quo_updt(args[[".X"]], pre = "R")
  )

  # new quosures
  args <- rlang::list2(
    !!! args,
    # Sigma cut-off bound
    hat_s_R = quo_updt(args[["R"]], pre = "hat_s"),
    # Predicted rare isotope count rate
    hat_X1 = quo_updt(args[["X1"]], pre = "hat")
  )

  # quantiles for cut-off bound
  fct_min <- qnorm((.alpha_level / 2))
  fct_max <- qnorm(1 - (.alpha_level / 2))

  dplyr::group_by(.IC, !!! gr_by) |>
    dplyr::mutate(
      !! args[["hat_s_R"]] := sd(!! args[["R"]]),
      !! args[["hat_X1"]] := !! args[["M_R"]] * !! args[["X2"]],
      flag = dplyr::if_else(
        dplyr::between(
          !! args[["R"]],
          unique(!! args[["M_R"]] + fct_min * !! args[["hat_s_R"]]),# mean - 2SD
          unique(!! args[["M_R"]] + fct_max * !! args[["hat_s_R"]]) # mean + 2SD
        ),
        "confluent",
        "divergent"
      ),
      flag = as.factor(.data$flag)
      ) |>
    dplyr::ungroup()
}
