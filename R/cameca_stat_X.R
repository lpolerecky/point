#' Reference Cameca X statistics
#'
#' The single ion count statistics as performed by the Cameca software. The ion
#' counts have not been corrected for systematic offsets by EM detection
#' devices
#'
#' @format A data frame with 54.000 rows and 11 variables:
#' \describe{
#'   \item{file.nm}{name of the original file}
#'   \item{correction block}{systematic corrections (named blocks in Cameca
#'   output)}
#'   \item{Mass#}{ion mass according to number of detection device}
#'   \item{Cumulated count}{Total ion count (same as N~tot~)}
#' }
#'
"cameca_stat_X"
