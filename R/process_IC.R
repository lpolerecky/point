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
#' @param .IC A tibble containing raw ion count data.
#' @param ... Currently not supported.
#' @param .N A variable constituting the ion counts.
#' @param .t A variable constituting the time increments.
#' @param .bl_t A variable or numeric value for the blanking time
#' @param .det Variable or character string or variable for the detection
#' system ("EM" or "FC").
#' @param .deadtime A numeric value for the deadtime of the EM system.
#' @param .thr_PHD A numeric value for the discriminator threshold of the EM.
#' system
#' @param .mean_PHD A value or numeric valie of the mean PHD value
#' @param .SD_PHD A value or numeric valie of standard deviation of the PHD
#' value
#' @param .hide A logical indicating whether only processed data should be
#' returned. If \code{TRUE} The raw data is contained as an attribute named
#' \code{"rawdata"}.
#'
#' @return A \code{tibble::\link[tibble:tibble]{tibble}()} containing the original
#' dataset and adds the variables: \code{Xt.rw}, ion count rates uncorrected for
#' detection device-specific biases; \code{Xt.pr}, ion count rates corrected for
#' detection device-specific biases; and \code{N.pr}, counts corrected for
#' detection device-specific biases.
#' @export
#' @examples
#' # Use point_example() to access the examples bundled with this package
#'
#' # raw data containing 13C and 12C counts on carbonate
#' tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # Processing raw ion count data
#' cor_IC(tb_rw)
cor_IC <-function(.IC, ..., .N = N.rw, .t = t.nm, .bl_t = tc.mt,
                  .det = det_type.mt, .deadtime = NULL, .thr_PHD = NULL,
                  .mean_PHD = mean_PHD.mt, .SD_PHD = SD_PHD.mt,
                  .hide = TRUE){

  stopifnot(tibble::is_tibble(.IC))

# Quoting the call (user-supplied expressions)
  args <- enquos(N = .N, t = .t, bl_t = .bl_t, det = .det, M_PHD = .mean_PHD,
                 SD_PHD = .SD_PHD)

# New quosure
  N.pr <- quo_updt(args[["N"]], post = "pr", update_post = TRUE)

# Unfold metadata
  if (ncol(select(.IC, ends_with(".mt"))) == 0) IC <- unfold(.IC) else IC <- .IC

# Manual setting of blank time
  if (is.numeric(get_expr(args[["bl_t"]]))) {
    IC <- mutate(IC, tc.mt = !! args[["bl_t"]])
    args[["bl_t"]] <- quo(tc.mt)
  }
# Manual setting of detector type
  if (is_character(get_expr(det))) {
    IC <- mutate(IC, det_type.mt = !! det)
    det <- quo(det_type.mt)
  }

# Time increments and count rate (time between measurements minus blanking time)
  IC <- mutate(
    IC,
    dt.rw = min(!! args[["t"]]) - !! args[["bl_t"]],
    Xt.pr = !! args[["N"]] / .data$dt.rw,
    !! N.pr := !! args[["N"]]
    )

# PHD correction (Yield) on counts and count rates
  if (!is.null(.thr_PHD)) {

# Manual setting of PHD params (mean and sd)
    if (is.numeric(get_expr(args[["M_PHD"]])) &
        is.numeric(get_expr(args[["SD_PHD"]]))) {
      IC <- mutate(
        IC,
        mean_PHD.mt = !! args[["M_PHD"]],
        SD_PHD.mt = !! args[["SD_PHD"]]
        )
      M_PHD <- quo(mean_PHD.mt)
      SD_PHD <- quo(SD_PHD.mt)
      }

    IC <- tidyr::nest(IC, data = -c(!! args[["M_PHD"]], !! args[["SD_PHD"]])) %>%
      mutate(
        Y.mt =
          purrr::map2_dbl(
            !! args[["M_PHD"]],
            !! args[["SD_PHD"]],
            purrr::possibly(cor_yield, NA_real_),
            x = NULL,
            thr_PHD = .thr_PHD,
            output = "Y"
            )
        ) %>%
      tidyr::unnest(cols = .data$data) %>%
      mutate(
        Xt.pr = .data$Xt.pr / .data$Y.mt,
        !! N.pr := .data$Xt.pr * .data$dt.rw
        )
    }

# Deadtime correction on counts and count rates
  if (!is.null(.deadtime)) {
    IC  <- mutate(
      IC,
      Xt.pr =
        if_else(
          !! args[["det"]] == "EM",
          cor_DT(.data$Xt.pr, .deadtime),
          .data$Xt.pr
          ),
      !! N.pr :=
        if_else(
          !! args[["det"]] == "EM",
          .data$Xt.pr * .data$dt.rw,
          !! N.pr
          )
      )
  }

# Output
  if (.hide) {
    IC <- fold(IC, type = c(".mt", ".rw"))
    return(IC)
    }
  return(IC)
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
#' @param x A numeric vector containing raw ion count data.
#' @param mean_PHD A numeric vector coontaining the mean PHD value.
#' @param SD_PHD A numeric vector coontaining the standard deviation of the
#' PHD.
#' @param thr_PHD A numeric value for the disrcriminator threshold of the EM
#' system.
#' @param deadtime A numeric value for the deadtime of the EM system.
#' @param output Character string indicating whether to return corrected count
#' rates (\code{"ct"}) or yield value (\code{"Y"}).
#'
#' @return A numeric vector with the corrected count rates.
#' @export
#' @examples
#' # Original count rate of a chemical species
#' x <- 30000
#'
#' # Corrected count rate for EM Yield with a threshold of 50 V
#' cor_yield(x, mean_PHD = 210, SD_PHD = 60, thr_PHD = 50)
#'
#' # Corrected count rate for a deadtime of 44 ns
#' cor_DT(x, 44)
cor_yield <- function(x = NULL, mean_PHD, SD_PHD, thr_PHD, output = "ct"){

# # Stop execution if threshold = 0 an return Xt
  if (thr_PHD == 0) {
    warning("PHD discrimnator threshold 0; count rate or a yield of 1 is returned.",
            call. = FALSE)
    }

# Lambda parameter
  lambda <- (2 * mean_PHD^2) / (SD_PHD^2 + mean_PHD)

# Probability parameter
  prob <- (SD_PHD^2 - mean_PHD) / (SD_PHD^2 + mean_PHD)
  prob <- if_else(prob < 0 | prob >= 1, NA_real_, prob)

  Y <- polyaAeppli::pPolyaAeppli(thr_PHD, lambda, prob, lower.tail = FALSE)

  if (output == "Y") return(Y)
  if (output == "ct") return(x / Y)
  }

#' @rdname  cor_yield
#'
#' @export
cor_DT <- function(x, deadtime) {
  # Stop execution if threshold = 0 an return x
  if (deadtime == 0) {
    return(x)
    } else {
      x <- x / (1 - (x * deadtime * 10^-9))
      return(x)
    }
  }

#-------------------------------------------------------------------------------
# Supporting functions (NOT EXPORTET)
#-------------------------------------------------------------------------------

# Function which updates quosures for subsequent tidy evaluation
quo_updt <- function(my_q, pre = NULL, post = NULL, update_post = FALSE){

  # Get expressions
  old_expr <- get_expr(my_q)
  # Get text
  old_chr <- expr_text(old_expr)

  # Update
  if (update_post & stringr::str_detect(old_chr , "\\.") ){
    old_chr <- stringr::str_split(old_chr, "\\.")[[1]][1]
  }

  # Separators
  if (is.null(pre) & is.null(post)) {
    warning("Quosure not updated")
    return(my_q)
  }
  if (!is.null(pre) & is.null(post)) {
    new_chr <- paste0(pre, "_", old_chr)
  }
  if (is.null(pre) & !is.null(post)) {
    new_chr <- paste0(old_chr, ".", post)
  }
  if (!is.null(pre) & !is.null(post)) {
    new_chr <- paste0(pre, "_", old_chr, ".", post)
  }

  # New expression from character (remove whitespace)
  new_expr <- parse_expr(stringr::str_replace_all(new_chr, " ", ""))
  # Update old quosure
  set_expr(my_q, new_expr)

}
