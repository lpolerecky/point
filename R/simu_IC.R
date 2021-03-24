#' Simulated ion count dataset
#'
#' A dataset containing simulated ion count ratios with no intra-isotope
#' variance ("ideal"), and intra-isotope variability following a asymmetric
#' offset and a symmetric gradient. The simulated dataset has a 60 per mille
#' offset in 13C/12C (on the VPDB-scale) in the case of a heterogeneous
#' isotopic composition. This was done with the function
#' \code{point::\link[point:simu_R]{simu_R}()}.
#'
#'
#' @format A data frame with 54.000 rows and 11 variables:
#' \describe{
#'   \item{type.nm}{name of the scenario for isotopic variance: ideal, asymmetric and symmetric}
#'   \item{trend.nm}{ionization efficiency trend}
#'   \item{base.nm}{main isotopic composition}
#'   \item{force.nm}{isotope anomoly}
#'   \item{t.nm}{time increments}
#'   \item{bl.nm}{count block number}
#'   \item{n.rw}{number of observation per analysis}
#'   \item{spot.nm}{hypothetical spot number, i.e. the individual analysis}
#'   \item{species.nm}{ion species name}
#'   \item{N.pr}{processed ion count}
#'   \item{Xt.pr}{processed ion count rate}
#' }
#' @source \code{point::\link[point:simu_R]{simu_R}()}
"simu_IC"
