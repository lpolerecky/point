#' Process raw ion count data
#'
#' \code{cor_IC} function for processing ion count data.
#'
#' The accuracy of pulsed ion counting is influenced by systematic errors which
#' depend on the ion counting system. Deadtime and EM yield are two
#' prominent effects for the electron multiplier systems. The deadtime refer to
#' the timewindow when the system does not register counts; this occurs when
#' incident ions strike the EM in a small enough time window in which the EM
#' channel is electronically paralysed. The EM yield is the ratio between the
#' number of output pulses counted after the EM  discriminator threshold and
#' the number of ions arriving at the EM. The latter can be gauged with the peak
#' height distribution (PHD) which is the probability for an EM output to have a
#' certain voltage amplitude.
#'
#' @param df A tibble containing raw ion count data
#' @param N A variable constituting the ion counts
#' @param t A variable constituting the time increments
#' @param Det A character string or variable for the detection system ("EM" or
#' "FC")
#' @param deadtime A numeric value for the deadtime of the EM system
#' @param thr_PHD A numeric value for the disrcriminator threshold of the EM
#' system
#'
#' @return A \code{\link[tibble:tibble]{tibble}} containing the original dataset
#' and the processed ion count data
#'
#' @examples
#' # Use point_example() to access the examples bundled with this package in the
#' # inst/extdata directory.
#'
#' # raw data containing 13C and 12C counts on carbonate
#' tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # processing raw ion count data
#' tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt, deadtime = 44, thr_PHD = 50)
#'
#' @export
cor_IC <-function(df, N, t, Det, deadtime = 44, thr_PHD = 50){

  stopifnot(tibble::is_tibble(df))
  stopifnot(is.numeric(deadtime))

  N <- enquo(N)
  t <- enquo(t)
  Det <- enquo(Det)

  tb.pr <- df %>%
# Time increments
             mutate(dt = min(!!t),
# Count rates
                    Xt.rw = !!N / dt) %>%
# Deadtime correction on count rates
             mutate(Xt.pr = if_else(!!Det == "EM",
                                    cor_DT(Xt.rw, deadtime),
                                    Xt.rw ),
# Deadtime correction on counts
                    N.pr = if_else(!!Det == "EM",
                                   as.integer(cor_DT(Xt.rw, deadtime) * dt),
                                   N.rw ),
# Yield correction on count rates
                    Xt.pr = if_else(!!Det == "EM",
                                    cor_yield(Xt.pr, mean_PHD, SD_PHD, thr_PHD),
                                    Xt.pr),
# Yield correction on counts
                    N.pr = if_else(!!Det == "EM",
                                   as.integer(cor_yield(Xt.pr,
                                                        mean_PHD,
                                                        SD_PHD,
                                                        thr_PHD)
                                   * dt),
                                   N.pr))

    return(tb.pr %>% select(-dt))
    }

#' Correct ion detection bias
#'
#' \code{cor_yield} function to correct counting bias related to EM Yield.
#' \code{cor_DT} function to correct counting bias related to deadtime
#'
#' The accuracy of pulsed ion counting is influenced by systematic errors which
#' depend on the ion counting system. Deadtime and EM yield are two
#' prominent effects for the electron multiplier systems. The deadtime refer to
#' the timewindow when the system does not register counts; this occurs when
#' incident ions strike the EM in a small enough time window in which the EM
#' channel is electronically paralysed. The EM yield is the ratio between the
#' number of output pulses counted after the EM  discriminator threshold and
#' the number of ions arriving at the EM. The latter can be gauged with the peak
#' height distribution (PHD) which is the probability for an EM output to have a
#' certain voltage amplitude.
#'
#' @param Xt A numeric vector containing raw ion count data
#' @param mean_PHD A numeric vector coontaining the mean PHD value
#' @param SD_PHD A numeric vector coontaining the standard deviation of the
#' PHD
#' @param thr_PHD A numeric value for the disrcriminator threshold of the EM
#' system
#' @param deadtime A numeric value for the deadtime of the EM system
#'
#' #' @return A numeric vector with the corrected count rates
#'
#' @examples
#' # Original count rate of a chemical species
#' x <- 30000
#'
#' # Corrected count rate for EM Yield with a threshold of 50 V
#' cor_yield(x, mean_PHD = 210, SD_PHD = 60, thr_PHD = 50)
#'
#' # Corrected count rate for a deadtime of 44 ns
#' cor_DT(x, 44)
#'
#' @export
cor_yield <- function(Xt, mean_PHD, SD_PHD, thr_PHD){

# Lambda parameter
  lambda <- (2 * mean_PHD^2) / (SD_PHD^2 + mean_PHD)

# Probability parameter
  prob <- (SD_PHD^2 - mean_PHD) / (SD_PHD^2 + mean_PHD)
  prob <- if_else(prob < 0 | prob >= 1, NA_real_, prob, NA_real_)

  l.params <- lst(a = lambda, b = prob,  c = rep(thr_PHD, length(b)), d = Xt)
  f.polya <- function(a, b, c, d){

    if (!(is.na(a) | is.na(b))){

      Y <- polyaAeppli::pPolyaAeppli(c,
                                     lambda = a,
                                     prob = b,
                                     lower.tail = FALSE)
      Xt.pr <- d / Y

      return(Xt.pr)

    }else{

      Xt.pr <- d

      return(Xt.pr)

    }
  }

  Xt.pr <- purrr::pmap_dbl(l.params, f.polya)

  return(Xt.pr)

  }

#' @rdname  cor_yield
#'
#' @export
cor_DT <- function(Xt, deadtime) {

  Xt.pr <- Xt / (1 - (Xt * deadtime * 10^-9))

  return(Xt.pr)

  }
