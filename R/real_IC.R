#' Real ion count dataset
#'
#' The dataset is generated with the \emph{Cameca NanoSIMS 50L} at the
#' Department of Earth Sciences at Utrecht University. The material is an
#' in-house reference (a belemnite rostra) is used to check
#' repeatability/external reproducibility. Ion detection for 7
#' individual species was performed solely with electron multipliers (EM) with
#' the main purpose of producing stable carbon isotope ratios (^13^C/^12^C).
#'
#' @format A data frame with 54.000 rows and 11 variables:
#' \describe{
#'   \item{file.nm}{name of the original file}
#'   \item{t.nm}{time increments (seconds)}
#'   \item{species.nm}{ion species name}
#'   \item{sample.nm}{sample substrate name}
#'   \item{bl.nm}{count block number}
#'   \item{Xt.pr}{processed ion count rate}
#'   \item{N.pr}{processed ion count}
#' }
#'
"real_IC"
