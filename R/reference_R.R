#' The reference materials used for common isotope systems after the
#' \href{https://wwwrcamnl.wr.usgs.gov/isoig/isopubs/itchch2.html}{USGS})
#'
#' A dataset containing the reference materials and their assigned values for
#' calibrating isotope systems. These values are used by
#' \code{point::\link[point:calib_R]{calib_R}()}.
#'
#'
#' @format A data frame with 10 rows and 4 variables:
#' \describe{
#'   \item{isotope}{name of the rare isotope}
#'   \item{ratio}{isotope ratio}
#'   \item{reference}{name of reference material}
#'   \item{value}{value of refrence material}
#' }
#' @source <https://wwwrcamnl.wr.usgs.gov/isoig/isopubs/itchch2.html>
"reference_R"
