#' Process raw ion count data
#'
#' \code{cor_IC} function for processing ion count data
#'
#' The accuracy of pulsed ion counting is influenced by systematic errors which
#' depending on the ion counting system. Dead time and EM yield are two
#' prominent effects for the electron multiplier systems. The dead time where
#' the system does not register counts; this occurs when to incident ions strike
#' the EM in a small enough time window in which the EM channel is
#' electronically paralysed. The EM yield is the ratio between the number of
#' output pulses counted after the EM  discriminator threshold and the number of
#' ions arriving at the EM. The latter can be gauged with the peak height
#' distribution (PHD) which is the probability for an EM output to have a
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
#' @return A \code{\link[tibble:tibble]{tibble}} containing processed
#' ion count data and metadata
#'
#' @examples
#' # Use point_example() to access the examples bundled with this package in the
#' # inst/extdata directory.
#'
#' # raw data containing 13C and 12C counts on carbonate
#' tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # processing raw ion count data
#' tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt)
#'
#' @export
cor_IC <-function(df, N, t, Det, deadtime = 44, thr_PHD = 50){

  stopifnot(tibble::is_tibble(df))
  stopifnot(is.numeric(deadtime))

  N <- enquo(N)
  t <- enquo(t)
  Det <- enquo(Det)

  df <- df %>%
# Time increments
          mutate(dt = min(!!t),
# Count rates
                 Xt.rw = !!N / dt) %>%
# Deadtime correction on count rates
          mutate(Xt.pr = if_else(!!Det == "EM",
                                 Xt.rw / (1 - (Xt.rw * deadtime * 10^-9)),
                                 Xt.rw ),
# Deadtime correction on counts
                 N.pr = if_else(!!Det == "EM",
                                (Xt.rw / (1 - (Xt.rw * deadtime * 10^-9))) * dt,
                                 Xt.rw )) %>%
# Yield correction on count rates
            mutate(Y = if_else(!!Det == "EM",
                               cor_yield(mean_PHD, SD_PHD, thr_PHD),
                               1),
                   Xt.pr = if_else(is.na(Y), Xt.pr,  Xt.pr / Y),
                   N.pr = if_else(is.na(Y), N.pr , Xt.pr * dt))

}

# Yield correction with Polya-Aeppli density probability approximation for PHD
cor_yield <- function(mean_PHD, SD_PHD, thr_PHD){

# Lambda parameter
  lambda <- (2 * mean_PHD^2) / (SD_PHD^2 + mean_PHD)

# Probability parameter
  prob <- (SD_PHD^2 - mean_PHD) / (SD_PHD^2 + mean_PHD)
  prob <- if_else(prob < 0 | prob >= 1, NA_real_, prob, NA_real_)

  l.params <- lst(a = lambda, b = prob,  c = rep(thr_PHD, length(b)))
  f.polya <- function(a, b, c){

   if (!(is.na(a) | is.na(b))){

       Y <- polyaAeppli::pPolyaAeppli(c,
                                     lambda = a,
                                     prob = b,
                                     lower.tail = FALSE)
     }else{

       Y <-  NA_real_

     }


  }

   Y <- purrr::pmap_dbl(l.params, f.polya)
   return(Y)
}

# Deadtime correction
cor_DT <- function(Xt.rw, deadtime) {

  Xt.pr <- Xt.rw / (1 - (Xt.rw * deadtime * 10^-9))

  return(Xt.pr)

  }


