#' Reference Cameca R statistics
#'
#' The ion count ratio statistics as performed by the Cameca software. The ion
#' counts have not been corrected for systematic offsets by EM detection
#' devices.
#'
#' @format A data frame with 54.000 rows and 11 variables:
#' \describe{
#'   \item{file.nm}{name of the original file}
#'   \item{correction block}{systematic corrections (named blocks in Cameca
#'   output)}
#'   \item{Ratio#}{ion mass ratio according to number of detection device}
#'   \item{Ratios}{Ion count ratio (same as R)}
#'   \item{Err_Poisson(\%)}{Predicted standard error of the mean (in percent)}
#'   \item{Err_mean(\%)}{Standard error of the mean (in percent)}
#'   \item{Khi2}{Reduced chi-squared}
#'   \item{SD_Block(\%)}{Block-wise standard deviation of R}
#' }
#'
"cameca_stat_R"
