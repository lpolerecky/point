#' QSA test
#'
#' \code{QSA_test} function to test for QSA.
#'
#' The accuracy of pulsed ion counting is influenced by systematic errors which
#' depend on the ion counting system. Quasi simultaneous arrival is one of
#' those potential errors that can also impact isotope ratios.
#'
#' @param .df A tibble containing processed ion count data.
#' @param .ion1 A character string constituting the heavy isotope ("13C").
#' @param .ion2 A character string constituting the light isotope ("12C").
#' @param ... Variables for grouping.
#' @param .Xt A variable constituting the ion count rate (defaults to
#' variables generated with \code{read_IC()}.)
#' @param .N A variable constituting the ion counts (defaults to variables
#' generated with \code{read_IC()}.).
#' @param .species A variable constituting the species analysed (defaults to
#' variables generated with \code{read_IC()}.).
#' @param .t A variable constituting the time of the analyses (defaults to
#' variables generated with \code{read_IC()}.).
#' @param .plot Currently not supported.
#'
#' @return A \code{\link[tibble:tibble]{tibble}} containing the original dataset
#' and adds the variables: \code{beta}, \code{t_QSA}, and \code{p_QSA} that
#' summarise the results of an linear model fitted by OLS (respectively; the
#' slope and the associated student's t test statistic and  p value) on the
#' ion count rates of the common isotope (as predictor) and the isotope ratio
#' (as dependent variable). The p value is for \code{beta} being different from
#' zero.
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
QSA_test <- function(.df, .ion1, .ion2, ..., .nest = NULL, .Xt = Xt.pr,
                     .N = N.pr, .species = species.nm, .t = t.nm, plot = TRUE){

  # # Quoting the call (user-supplied expressions)
  Xt <- enquo(.Xt)
  N <- enquo(.N)
  species <- enquo(.species)
  t <- enquo(.t)
  gr_by <- enquos(...)
  nest <- enquo(.nest)
  Xt1 <- quo_updt(Xt, post = .ion1) # count rate
  Xt2 <- quo_updt(Xt, post = .ion2) # count rate
  R_Xt <- quo_updt(Xt, pre = "R")

  df <- zeroCt(.df, .ion1 , .ion2, !!! gr_by, .N = !!N, .species = !!species,
               .t = !!t) %>%
    cov_R(c(.ion1, .ion2), !!! gr_by, .species = !!species, .t = !!t) %>%
    mutate("R_{{Xt}}" := !! Xt1  / !! Xt2 )

  df_lm <- tidyr::nest(df, data = -c(!!! gr_by)) %>%
    mutate(
      lm_out =
        purrr::map(
          data,
          purrr::possibly(mlm_QSA, NA),
          .Xt1 = R_Xt,
          .Xt2 = Xt2
        )
      ) %>%
    tidyr::unnest_wider(lm_out) %>%
    tidyr::unnest(cols = data)

  if (is_symbol(get_expr(nest))) {
    # Groups for nested data
    nest_gr <- gr_by[!sapply(gr_by, as_name) %in% as_name(nest)]

    df_mlm <- tidyr::nest(df, data = -c(!!! nest_gr)) %>%
      mutate(
        mlm_out =
          purrr::map(
            data,
            purrr::possibly(mlm_QSA, NA),
            .Xt1 = R_Xt,
            .Xt2 = Xt2,
            .group = nest
          )
        ) %>%
      select(-data) %>%
      tidyr::unnest_wider(mlm_out)

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

mlm_QSA <- function(.df, .Xt1, .Xt2, .group = NULL) {

  # lm  model
  if (is.null(get_expr(.group))) {
    lm_QSA <- lm_form(.df, .Xt1, .Xt2)
    td <- broom::tidy(lm_QSA) %>% mutate(effect = "fixed")
    min <-  min(fitted(lm_QSA))
    max <-  max(fitted(lm_QSA))
    label <- .Xt2
  }
  # mlm  model
  if (!is.null(get_expr(.group))) {
    mlm_QSA <- lm_form(.df, .Xt1, .Xt2, nest = .group, type = "QSA")
    td <- broom.mixed::tidy(mlm_QSA)
    min <-  min(fitted(mlm_QSA))
    max <-  max(fitted(mlm_QSA))
    label <- .group
  }

  ls_QSA <- lst(
    # model params
    "alpha_{{label}}" := pull(filter(td, effect == "fixed" & term == "(Intercept)"), estimate),
    "beta_{{label}}" := pull(filter(td, effect == "fixed" & term == as_name(.Xt2)), estimate),
    "t_{{label}}" := pull(filter(td, effect == "fixed" & term == as_name(.Xt2)), statistic),
    "p_{{label}}" := pull(filter(td, effect == "fixed" & term == as_name(.Xt2)), `p.value`),
    # modelled delta value
    "delta_{{label}}" := (min / max - 1) * 1e3
  )

  return(  ls_QSA)

}



