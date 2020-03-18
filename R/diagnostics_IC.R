#' Diagnostics raw ion count data
#'
#' \code{diag_R} function for propagation of uncertainty for single ions
#'
#' Ion count data consisting of time-incremented integer values are
#'
#' @param df A tibble containing processed ion count data
#' @param method Character string for the type of diagnostics
#'
#' @return A tibble containing diagnostics
#'
#' @examples
#' # Use system.file() to access the examples bundled with this package in the
#' # inst/extdata directory. The examples directories are named:
#' # 2020-01-17-TREASURE and "2018-01-19-GLENDON"
#'
#'
#' tb.dia <- diag_R(tb.pr, method = "standard", Moreargs = quos(Xt = Xt.pr,
#'                  N = N.pr, species = species.nm, ion1 = "13C", ion2 = "12C",
#'                  .named = TRUE), file.nm)
#'
#' @export
diag_R <- function(df, method = "standard", Moreargs = quos(NULL), ...){

  args <- Moreargs
  gr_by <- enquos(...)

# descriptive an predictive statistics for 13C/12C ratios
  tb.aug <- zeroCt(df, !!args[["N"]], !!args[["species"]], as_name(args[["ion1"]]), as_name(args[["ion2"]]), !!! gr_by)  %>%
              stat_R(. , !!args[["Xt"]], !!args[["N"]], species = !!args[["species"]], ion1 = as_name(args[["ion1"]]), ion2 = as_name(args[["ion2"]]),
                     !!! gr_by, bl, latex = FALSE, output = "complete") %>%
              group_by(!!! gr_by, bl) %>%
              mutate(lower = !!quo_updt(my_q = args[["Xt"]] , x = "M_R") - 2 * !!quo_updt(my_q = args[["Xt"]] , x = "S_R"),
                     upper = !!quo_updt(my_q = args[["Xt"]] , x = "M_R") + 2 * !!quo_updt(my_q = args[["Xt"]] , x = "S_R"),
                     flag = if_else(
                                     between(!!quo_updt(my_q = args[["Xt"]] , x = "R"),
                                             unique(lower),  # mean - 2SD
                                             unique(upper)), # mean + 2SD
                                     "good",
                                     "bad")
                                    ) %>%
# SD blocks as in CAMECA data files (total R variance)
              group_by(!!! gr_by) %>%
              mutate(SD_bl = sd(!!quo_updt(my_q = args[["Xt"]] , x = "R")) / mean(!!quo_updt(my_q = args[["Xt"]] , x = "R")) * 100 )  %>%
              ungroup() %>%
              select(ID, SD_bl , flag)

# Cameca style augmented datafile
  tb.aug <- left_join(df, tb.aug

                      , by = "ID") # %>% select(-ID)

}
