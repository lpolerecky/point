#' Chemical species names for plots and tables
#'
#' \code{ion_labeller} and \code{R_labeller} converts a character string
#' containing chemical species names in a latex string or expression.
#'
#' This functions converts chemical species names of the form, e.g. `"12C"`,
#' `"13C2"`, `"12C 14N"`, or `"12C-14N"` to a character string which can be
#' parsed in Latex to species names with appropriate superscripts on th left for
#' mass and subscripts for stochiometry on the right.
#'
#' @param ion A chemical species character string
#' @param ion1 A chemical species character string (for isotope ratios the heavy
#' isotope).
#' @param ion2 A chemical species character string (for isotope ratios the light
#' isotope).
#' @param label Character string indicating whether the output should be latex
#' or an expression.
#'
#' @return A character string parsable in Latex or expression for usage in
#' plots.
#' @export
#' @examples
#'
#' # plot some ion count data
#' library(ggplot2)
#' ggplot() +
#'    geom_blank() +
#'    ylab(ion_labeller("12C2-40Ca", "expr")) +
#'    xlab(R_labeller("12C2-40Ca", "13C2-40Ca", "expr"))
#'
ion_labeller <- function(ion, label = "latex") {

# Regex for element names
  el_reg <- "((?<=[[:digit:]])[[:upper:]]{1})(?<=[[:upper:]])[[:lower:]]{1}|(?<=[[:digit:]])[[:upper:]]{1}"
# Regex mass number
  ms_reg <- "[[:digit:]]{1,2}(?=[[:upper:]])"
# Regex stochiometry
  st_reg <- "(?<=[[:alpha:]])([[:digit:]]|\\s|[[:punct:]])"

  suppressWarnings(
  ls_chr <- purrr::map(
    lst(a = ms_reg, b = el_reg, c = st_reg),
    ~stringr::str_extract_all(ion, .x),
    ) %>%
    purrr::flatten() %>%
    purrr::transpose()
  )

  if (label == "latex") {

    ls_ion <- purrr::map(
      ls_chr,
      ~paste0(
        "$\\phantom{,}^{",
        .x[["a"]],
        "}$",
        .x[["b"]],
        "$_{",
        .x[["c"]],
        "}$"
        )
      )

    lb <- purrr::reduce(ls_ion, paste, sep = "--")
    return(lb)
  }
  if (label == "expr") {

    ls_ion <- purrr::map(ls_chr, ~substitute(""^a*b[c] , env = .x))

    lb <- purrr::reduce(
      ls_ion,
      function(x, y) substitute(a *"-"* b, env = lst(a = x, b = y))
      )
    return(lb)
    }
}
#' @rdname ion_labeller
#'
#' @export
R_labeller <- function(ion1, ion2, label = "latex"){

  ls_R <- purrr::map(c(ion1, ion2), ion_labeller, label)
  if (label == "latex") {
    lb <- paste(ls_R[[1]], ls_R[[2]], sep = "/")
    return(lb)
    }
  if (label == "expr") {
    lb <- substitute(a *"/ "* b, env = lst(a = ls_R[[1]], b = ls_R[[2]]))
    return(lb)
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
#' @param .df A tibble containing processed ion count data.
#' @param .ion1 A character string constituting the heavy isotope ("13C").
#' @param .ion2 A character string constituting the light isotope ("12C").
#' @param ... Variables for grouping.
#' @param .species A variable constituting the species analysed.
#' @param .t A variable constituting the time of the analyses.
#' @param .preserve A logical whether to preserve ID variable.
#'
#' @return A \code{\link[tibble:tibble]{tibble}} in wide format
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
#' cov_R(tb_pr, "13C", "12C", file.nm)
#'
cov_R <- function(.df, .ion1, .ion2, ..., .species = species.nm, .t = t.nm,
                  .preserve = FALSE){

  t <- enquo(.t)
  species <- enquo(.species)
  gr_by <- enquos(...)

  # Filtering
  df <- filter(.df, !! species == .ion1 | !! species == .ion2) %>%
    mutate(!! species := str_replace_all(!! species, " ", ""))

  # Wide format
  pivot_wider(
    df,
    c(!!!gr_by, !!t, !!species),
    names_from = !! species,
    values_from = -c(!!!gr_by, !!t, !!species),
    # values_fn = length,
    names_sep = "."
    )

  # if (.preserve) return(df) else return(select(df, -.data$ID))
}

#' Remove analytical runs with zero counts
#'
#' \code{zeroCt} removes analytical runs for isotope ratios that contain zero
#' counts.
#'
#' This functions removes analytical runs with zero counts for calculating
#' isotope ratios.
#'
#' @param .df A tibble containing processed ion count data.
#' @param .ion1 A character string constituting the heavy isotope ("13C").
#' @param .ion2 A character string constituting the light isotope ("12C").
#' @param ... Variables for grouping.
#' @param .N A variable constituting the ion counts.
#' @param .species A variable constituting the species analysed.
#' @param .t A variable constituting the time of the analyses.
#' @param warn A logical indicating whether to produce a warning.
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
#' # Vectors of isotope ratios
#' ion1 <-  c("13C", "12C 13C", "13C 14N", "12C 14N", "12C")
#' ion2 <-  c("12C", "12C2", "12C 14N", "40Ca 16O", "40Ca 16O")
#'
#' tb_pr <- purrr::map2(ion1, ion2, ~zeroCt(tb_pr, .x, .y, file.nm))
zeroCt <- function(.df, .ion1, .ion2, ..., .N = N.pr, .species = species.nm,
                   .t = t.nm, .warn = TRUE){

  N <- enquo(.N)
  t <- enquo(.t)
  species <- enquo(.species)
  gr_by <- enquos(...)

  df <- filter(.df, !!species == .ion1 | !!species == .ion2)

  if (.warn){
    if (any(select(df, !!N) %>% pull(!!N) == 0)) {
      warning("Zero counts present and removed")
    }
  }

  ls_0 <- filter(df, !! N == 0) %>%
    select(!!! gr_by)

  anti_join(df, ls_0, by = sapply(gr_by, as_name))
}

#-------------------------------------------------------------------------------
# Not exportet
#-------------------------------------------------------------------------------


# Function for building IDs
ID_builder <- function(.df, ..., .t = t.nm, .species = species.nm){

  time <- enquo(.t)
  species <- enquo(.species)
  gr_by <- enquos(...)

  tidyr::unite(
    .df,
    col = "ID",
    c(!!! gr_by, !! species, !! time),
    sep = "/",
    remove = FALSE
    )

  }

