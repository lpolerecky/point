#' Chemical species names for plots and tables
#'
#' \code{ion_labeller} and \code{R_labeller} converts a character string
#' containing chemical species names in a latex string or expression. The
#' \code{stat_labeller} function generates labels for statistics in tables
#' (latex) or on plots (expressions).
#'
#' This functions converts chemical species names of the form, e.g. `"12C"`,
#' `"13C2"`, `"12C 14N"`, or `"12C-14N"` to a character string which can be
#' parsed in Latex to species names with appropriate superscripts on th left for
#' mass and subscripts for the index on the right.
#'
#' @param ion A chemical species character string
#' @param ion1 A chemical species character string (for isotope ratios the rare
#' isotope).
#' @param ion2 A chemical species character string (for isotope ratios the
#' common isotope).
#' @param var A character string for the variable, either \code{"X"} for single
#' ions or \code{"R"} for ion ratios.
#' @param org A character string for the origin of a derived variable, e.g., R.
#' @param stat A character string for the statistic following convention of e.g.
#' \code{point::names_stat_X}.
#' @param value The numeric value for the statistic result.
#' @param label Character string indicating whether the output should be
#' \code{"latex"} or an expression (\code{"expr"}).
#'
#' @return A character string parsable in Latex or expression for usage in
#' plots.
#' @export
#' @examples
#' # plot some ion count data
#' library(ggplot2)
#' ggplot() +
#'  geom_blank() +
#'  ylab(ion_labeller("12C2-40Ca", "expr")) +
#'  xlab(R_labeller("12C2-40Ca", "13C2-40Ca", "expr"))
#'
ion_labeller <- function(ion, label = "latex") {

  # Regex for element names
  el_reg <- stringr::regex("
    ((?<=[:digit:])[:upper:]{1}?)(?<=[:upper:])[:lower:]{1}? # 2 letter element
    |(?<=[:digit:])[:upper:]{1}? # 1 letter element
    ", comments = TRUE)
  # Regex mass number
  ms_reg <- stringr::regex("
    [:digit:]{1,3}?(?=[:upper:]) # mass number in front of element
    ", comments = TRUE)
  # Regex element index
  st_reg <- stringr::regex("
    [:digit:]*?(?=[:punct:]) # index followed by punctuation
    |[:digit:]*?(?![:graph:]) # index followed by anything but graphics
    ", comments = TRUE)

  suppressWarnings(
    ls_chr <- purrr::map(
      tibble::lst(a = ms_reg, b = el_reg, c = st_reg),
      ~stringr::str_extract_all(ion, .x),
      ) %>%
      purrr::flatten() %>%
      purrr::transpose()
  )

  if (label == "latex" | label == "webtex") {
    paste_ion <- function(ls) {
      if (length(ls[["c"]]) != 0) {
        if (!stringr::str_detect(ls[["c"]], "[[:digit:]]")) ls[["c"]] <- NULL
      }
      paste0("${}^{", ls[["a"]], "}\\mathrm{", ls[["b"]], "}_{", ls[["c"]],
             "}$")
    }
    ls_ion <- purrr::map(ls_chr, paste_ion)
    purrr::reduce(ls_ion, paste0)
    } else if (label == "expr") {
    ls_ion <- purrr::map(ls_chr, ~ substitute("" ^ a * b[c] , env = .x))
    purrr::reduce(ls_ion, ~substitute(a * b, env = list(a = .x, b = .y)))
  }
}
#' @rdname ion_labeller
#'
#' @export
R_labeller <- function(ion1, ion2, label = "latex") {
  ls_R <- purrr::map(c(ion1, ion2), ion_labeller, label)
  if (label == "latex" | label == "webtex") {
    paste(ls_R[[1]], ls_R[[2]], sep = "") %>%
      stringr::str_replace("\\$\\$", "/")
  } else if (label == "expr") {
    substitute(a *"/ "* b, env = list(a = ls_R[[1]], b = ls_R[[2]]))
  }
}
#' @rdname ion_labeller
#'
#' @export
stat_labeller <- function(var, org, stat, value, label = "latex"){

  if (label == "latex" | label == "webtex") {
    if (stat == "n") return("$n$")
    if (stat == "tot") return("$N_{tot}$")

    if (stat == "F") return(paste0("$F_{", var,"}$"))
    if (stringr::str_detect(stat, "M")) {
      stat_chr <- paste0("\\bar{", var,"}")
      # external precision
      if (stringr::str_detect(org, "^M$|^Ntot$")) {
        stat_chr <- paste0("\\bar{", stat_chr,"}")
      }
      # mlm external precision
      if (stringr::str_detect(org, "^M_R$")) {
        stat_chr <- paste0("\\hat{", stat_chr,"}")
        if (stringr::str_detect(stat, "hat")) {
          stat <- stringr::str_replace(stat, "hat", "") # remove hat above stat
        }
      }
    } else {
      stat_chr <- var
      # external precision
      if (stringr::str_detect(org, "^M$|^Ntot$")) {
        stat_chr <- paste0("\\bar{", stat_chr,"}")
      }
      # mlm external precision
      if (stringr::str_detect(org, "^M_R$")) {
        stat_chr <- paste0("\\bar{", stat_chr,"}")
      }
    }

    # p values, AIC
    if (stat == "p") return(paste0("$p_{", stat_chr,"}$"))
    if (stat == "dAIC") return(paste0("$\\Delta AIC_{", stat_chr, "}$"))
    if (stat == "chi2") return(paste0("$\\chi^{2}_{",  stat_chr,"}$"))
    # Prefix of relative and normal standard deviations and errors
    if (stringr::str_detect(stat, "S")) {
      stat_chr <- paste0("_{", stat_chr,"}")
      if (stringr::str_detect(stat, "R")) {
        sd_prefix <- "\\epsilon"
      } else {
        sd_prefix <- "s"
      }
      # Hat of predicted standard deviations and errors
      if (stringr::str_detect(stat, "hat")) {
        stat_chr <- paste0("\\hat{", sd_prefix ,"}", stat_chr)
      } else {
        stat_chr <- paste0(sd_prefix, stat_chr)
      }
    }
    #  Per mille signs  of relative standard deviations and errors
    if (stringr::str_detect(stat, "R")) {
      if (label == "webtex") return(paste0("$", stat_chr, "$ (\u2030)"))
      if (label == "latex") {
        return(paste0("$", stat_chr, "$ (\\text{\\textperthousand})"))
      }
    } else {
      return(paste0("$", stat_chr, "$"))
    }
  }
  # plot labels as expressions
  if (label == "expr") {
    var <- rlang::parse_expr(var)
    if (stat == "n") {
      return(substitute(n == ~ a, list(a = sprintf(fmt = "%.0f", value))))
    }
    if (stat == "tot") {
      return(substitute(N[tot] == ~ a, list(a = sprintf(fmt = "%.0f", value))))
    }
    if (stat == "chi2") {
      return(substitute(chi^2 == ~ a, list(a = sprintf(fmt = "%.1f", value))))
    }
    if (stringr::str_detect(stat, "R")) {
      sd_prefix <- expression(epsilon)
      value <- paste(sprintf(fmt = "%.1f", value), "(\u2030)")
    } else {
      sd_prefix <- expression(s)
      value <- sprintf(fmt = "%.3f", value)
    }
    if (stringr::str_detect(stat, "hat")) {
      sd_prefix <- substitute(hat(a), list(a = sd_prefix))
    }
    if (stringr::str_detect(stat, "M")) {
        if (nchar(stat) == 1) {
          return(substitute(bar(a) == ~ b, list(a = var, b = value)))
        } else {
          sym_chr <- substitute(bar(a), list(a = var))
        }
      } else {
        sym_chr <- var
      }
    if (stringr::str_detect(stat, "S")){
      stat_chr <- substitute(
        a[b] == ~ c,
        list(a = sd_prefix, b = sym_chr, c = value)
      )
      return(stat_chr)
    }
  }
}

#' Function for co-variate conversion of isotope systems
#'
#' \code{cov_R} create a wide format tibble for an isotope pair
#'
#' This functions converts the long format data frame to a wide format tibble
#' for an isotope pair based on an unique identifier for the time of
#' measurement. The data can be linked by a combination of three variable unique
#' for the analyses; the file name, the chemical species name and the time
#' increment of measurement. Pay attention when using in combination with
#' `zeroCt()`.
#'
#' @param .IC A tibble containing processed ion count data.
#' @param .ion A character string or vector constituting ion names.
#' @param ... Variables for grouping.
#' @param .species A variable constituting the species analysed.
#' @param .t A variable constituting the time of the analyses.
#' @param .preserve A logical whether to preserve ID variable.
#'
#' @return A \code{tibble::\link[tibble:tibble]{tibble}} in wide format
#' @export
#' @examples
#'
#' # raw data containing 13C and 12C counts on carbonate
#' tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # Processing raw ion count data
#' tb_pr <- cor_IC(tb_rw)
#'
#' # wide format
#' cov_R(tb_pr, c("13C", "12C"), file.nm)
#'
cov_R <- function(.IC, .ion, ..., .species = NULL, .t = NULL,
                  .preserve = FALSE){

  # Quoting the call (user-supplied expressions)
  args <- inject_args(
    .IC,
    enquos(.t = .t, .species = .species),
    type = c("processed", "group")
  )

  gr_by <- enquos(...)

  # Observation per group
  obs_gr <- dplyr::group_size(dplyr::group_by(.IC, !!! gr_by, !! args[[".species"]]))
  # Distinct observation for time increments
  obs_t <- dplyr::n_distinct(dplyr::pull(.IC, !! args[[".t"]]))
  # Check if t steps is consistent with grouping otherwise create new ID
  if (any(obs_gr %% obs_t != 0)) {
    .IC <- dplyr::group_by(.IC, !!! gr_by, !! args[[".species"]]) %>%
      dplyr::mutate(!! args[[".t"]] := dplyr::row_number()) %>%
      dplyr::ungroup()
  }
  # Remove white space in ion names and add underscore for polyatomic species
  IC <- dplyr::mutate(
    .IC,
    !! args[[".species"]] := ion_trim(!! args[[".species"]])
  ) %>%
    dplyr::filter(!! args[[".species"]] %in% sapply(.ion, ion_trim))

  # Wide format
  tidyr::pivot_wider(
    IC,
    names_from = !! args[[".species"]],
    values_from = -c(!!! gr_by, !! args[[".t"]], !! args[[".species"]]),
    names_sep = "."
  )
}
#' Remove analytical runs with zero counts
#'
#' \code{zeroCt} removes analytical runs for isotope ratios that contain zero
#' counts.
#'
#' This functions removes analytical runs with zero counts for calculating
#' isotope ratios.
#'
#' @param .IC A tibble containing processed ion count data.
#' @param .ion1 A character string constituting the rare isotope ("13C").
#' @param .ion2 A character string constituting the common isotope ("12C").
#' @param ... Variables for grouping.
#' @param .N A variable constituting the ion counts.
#' @param .species A variable constituting the species analysed.
#' @param .warn A logical indicating whether to produce a warning.
#'
#' @return A \code{\link[tibble:tibble]{tibble}} containing the single ion count
#' dataset for the specified ion ratio. The grouping variable specifies on
#' which level the original dataset is subsetted.
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
#' # Remove analyses with zero counts
#' tb_pr <- zeroCt(tb_pr, "13C", "12C", file.nm)
zeroCt <- function(.IC, .ion1, .ion2, ..., .N = NULL, .species = NULL,
                   .warn = TRUE){

  # Quoting the call (user-supplied expressions)
  args <- inject_args(
    .IC,
    enquos(.N = .N, .species = .species),
    type = c("processed", "group")
  )
  gr_by <- enquos(...)

  # Remove white space in ion names and add underscore for polyatomic species
  .ion1 <- ion_trim(.ion1)
  .ion2 <- ion_trim(.ion2)
  IC <- dplyr::mutate(
    .IC,
    !! args[[".species"]] := ion_trim(!! args[[".species"]])
  ) %>%
    dplyr::filter(
      !! args[[".species"]] == .ion1 | !! args[[".species"]] == .ion2
    )

  if (isTRUE(.warn)) {
    if (any(dplyr::pull(IC, !! args[[".N"]]) == 0)) {
      warning("Zero counts present and removed.", call. = FALSE)
    }
  }

  ls_0 <- dplyr::filter(IC, !! args[[".N"]] == 0) %>%
    dplyr::select(!!! gr_by)
  IC <- dplyr::anti_join(IC, ls_0, by = sapply(gr_by, as_name))
  if(nrow(IC) == 0) {
    warning(
      "No more data left after removing zero count analysis.",
      call. = FALSE
    )
  }
  # Return
  IC
}

#-------------------------------------------------------------------------------
# Not exportet
#-------------------------------------------------------------------------------
# trim white space from ion character strings
ion_trim <-function(ion) {
  stringr::str_remove_all(stringr::str_trim(ion), "\\s")
}
