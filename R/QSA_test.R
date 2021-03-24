#' QSA test
#'
#' \code{QSA_test} function to test for QSA.
#'
#' The accuracy of pulsed ion counting is influenced by systematic errors which
#' depend on the ion counting system. Quasi simultaneous arrival is one of
#' those potential errors that can also impact isotope ratios.
#'
#' @param .IC A tibble containing processed ion count data.
#' @param .ion1 A character string constituting the heavy isotope ("13C").
#' @param .ion2 A character string constituting the light isotope ("12C").
#' @param ... Variables for grouping.
#' @param .nest A variable hat identifies a series of analyses to calculate
#' the significance of QSA.
#' @param .X A variable constituting the ion count rate (defaults to
#' variables generated with \code{read_IC()}.)
#' @param .N A variable constituting the ion counts (defaults to variables
#' generated with \code{read_IC()}.).
#' @param .species A variable constituting the species analysed (defaults to
#' variables generated with \code{read_IC()}.).
#' @param .t A variable constituting the time of the analyses (defaults to
#' variables generated with \code{read_IC()}.).
#' @param .plot Currently not supported.
#'
#' @return A \code{tibble::\link[tibble:tibble]{tibble}()} containing the
#' original dataset and adds the variables: \code{beta}, \code{t_QSA}, and
#' \code{p_QSA} that summarise the results of an linear model fitted by OLS
#' (respectively; the slope and the associated student's t test statistic and
#' p value) on the ion count rates of the common isotope (as predictor) and the
#' isotope ratio (as dependent variable). The p value is for \code{beta} being
#' different from zero.
#' .
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
#' # QSA test
#' tb_QSA <- QSA_test(tb_pr, "13C", "12C", file.nm)
QSA_test <- function(.IC, .ion1, .ion2, ..., .nest = NULL, .X = Xt.pr,
                     .N = N.pr, .species = species.nm, .t = t.nm, .plot = TRUE){

  # Quoting the call (user-supplied expressions)
  args <- enquos(.X = .X, .N = .N, .species = .species, .t = .t, .nest  = .nest)

  # Grouping
  gr_by <- enquos(...)
  nest <- enquo(.nest)

  # Updated quosures
  X1 <- quo_updt(args[[".X"]], post = .ion1) # count rate
  X2 <- quo_updt(args[[".X"]], post = .ion2) # count rate
  R_X <- quo_updt(args[[".X"]], pre = "R")

  IC <- zeroCt(.IC, .ion1, .ion2, !!! gr_by, .N = !! args[[".N"]],
               .species = !! args[[".species"]]) %>%
    cov_R(c(.ion1, .ion2), !!! gr_by, .species = !! args[[".species"]],
          .t = !! args[[".t"]]) %>%
    mutate(!!R_X  := !! X1  / !! X2 )
  return(IC)

  df_lm <- tidyr::nest(df, data = -c(!!! gr_by)) %>%
    mutate(
      lm_out =
        purrr::map(
          .data$data,
          purrr::possibly(mlm_QSA, NA),
          .X1 = R_X,
          .X2 = X2
          )
      ) %>%
    tidyr::unnest_wider(.data$lm_out) %>%
    tidyr::unnest(cols = .data$data)

  if (is_symbol(get_expr(args[[".nest"]]))) {
    # Groups for nested data
    nest_gr <- gr_by[!sapply(gr_by, as_name) %in% as_name(args[[".nest"]])]

    df_mlm <- tidyr::nest(IC, data = -c(!!! nest_gr)) %>%
      mutate(
        mlm_out =
          purrr::map(
            .data$data,
            purrr::possibly(mlm_QSA, NA),
            .X1 = R_X,
            .X2 = X2,
            .group = nest
          )
        ) %>%
      select(-.data$data) %>%
      tidyr::unnest_wider(.data$mlm_out)

    ls_mlm <- lst(df_lm, df_mlm)
    # Prepare output
    df <- purrr::reduce(ls_mlm, left_join, by = sapply(nest_gr, as_name))
    return(df)
    }
  return(df_lm)
}


#-------------------------------------------------------------------------------
# lm and mlm model QSA
#-------------------------------------------------------------------------------

mlm_QSA <- function(.IC, .X1, .X2, .group = NULL) {

  # lm  model
  if (is.null(get_expr(.group))) {
    lm_QSA <- lm_form(.IC, .X1, .X2)
    td <- broom::tidy(lm_QSA) %>% mutate(effect = "fixed")
    min <-  min(fitted(lm_QSA))
    max <-  max(fitted(lm_QSA))
    label <- .X2
  }
  # mlm  model
  if (!is.null(get_expr(.group))) {
    mlm_QSA <- lm_form(.IC, .X1, .X2, nest = .group, type = "QSA")
    td <- broom.mixed::tidy(mlm_QSA)
    min <-  min(fitted(mlm_QSA))
    max <-  max(fitted(mlm_QSA))
    label <- .group
  }

  ls_QSA <- lst(
    # model params
    "alpha_{{label}}" :=
      pull(
        filter(td, .data$effect == "fixed" & .data$term == "(Intercept)"),
        .data$estimate
        ),
    "beta_{{label}}" :=
      pull(
        filter(td, .data$effect == "fixed" & .data$term == as_name(.X2)),
        .data$estimate
        ),
    "t_{{label}}" :=
      pull(
        filter(td, .data$effect == "fixed" & .data$term == as_name(.X2)),
        .data$statistic
        ),
    "p_{{label}}" :=
      pull(
        filter(td, .data$effect == "fixed" & .data$term == as_name(.X2)),
        .data$`p.value`
        ),
    # modelled delta value
    "delta_{{label}}" := (min / max - 1) * 1e3
  )
  return(ls_QSA)
}



