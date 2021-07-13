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
#' (in milliseconds).
#' @param .det Variable or character string or variable for the detection
#' system ("EM" or "FC").
#' @param .deadtime A numeric value for the deadtime of the EM system with
#' units nanoseconds.
#' @param .thr_PHD A numeric value for the discriminator threshold of the EM.
#' system.
#' @param .M_PHD A variable or numeric value of the mean PHD value.
#' @param .SD_PHD A variable or numeric value of standard deviation of
#' the PHD value.
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
cor_IC <-function(.IC, ..., .N = NULL, .t = NULL, .bl_t = NULL,
                  .det = NULL, .deadtime = NULL, .thr_PHD = NULL,
                  .M_PHD = NULL, .SD_PHD =  NULL,
                  .hide = TRUE){

  # Checks
  stat_validator(.IC)

  # Unfold metadata
  if (ncol(select(.IC, ends_with(".mt"))) == 0) .IC <- unfold(.IC)

  # Quoting the call (user-supplied expressions)
  args <- enquos(.N = .N, .t = .t, .bl_t = .bl_t, .det = .det,
                 .M_PHD = .M_PHD, .SD_PHD = .SD_PHD)

  # Add custom arguments to tibble
  .IC <- supply_args(args, .IC)

  # Inject missing args
  args <- inject_args(.IC, args, type = c("raw", "group", "meta"))

  # Argument check
  argument_check(.IC, args, "raw")

  # New quosure
  N.pr <- quo_updt(args[[".N"]], post = "pr", update_post = TRUE)

  # Time increments (time between measurements minus blanking time)
  .IC <- mutate(
    .IC,
    dt.rw = min(!! args[[".t"]]) - !! args[[".bl_t"]] / 1e3,
    Xt.pr = !! args[[".N"]] / .data$dt.rw,
    !! N.pr := !! args[[".N"]]
    )

  # PHD correction (Yield) on counts and count rates
  if (cor_check(.thr_PHD, args, "PHD")) {
    # check if detector type is known
    det_check(args[[".det"]])
    .IC <- tidyr::nest(
      .IC,
      data = -c(!!args[[".det"]], !!args[[".M_PHD"]], !!args[[".SD_PHD"]])
      ) %>%
      mutate(
        Y.mt =
          if_else(
            !!args[[".det"]] == "EM",
            purrr::map2_dbl(
              !!args[[".M_PHD"]],
              !!args[[".SD_PHD"]],
              purrr::possibly(cor_yield, NA_real_),
              x = NULL,
              thr_PHD = .thr_PHD,
              output = "Y"
              ),
            1
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
    # check if detector type is known
    det_check(args[[".det"]])
    .IC  <- mutate(
      .IC,
      Xt.pr =
        if_else(
          !! args[[".det"]] == "EM",
          cor_DT(.data$Xt.pr, .deadtime),
          .data$Xt.pr
          ),
      !! N.pr := .data$Xt.pr * .data$dt.rw
      )
  }

  # Output
  if (.hide) .IC <- fold(.IC, type = c(".mt", ".rw"))
  .IC
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
#' @param mean_PHD A numeric vector containing the mean PHD value.
#' @param SD_PHD A numeric vector containing the standard deviation of the
#' PHD.
#' @param thr_PHD A numeric value for the disrcriminator threshold of the EM
#' system.
#' @param deadtime A numeric value for the dead-time of the EM system
#' with units nanoseconds.
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
    warning(
      "PHD discrimnator threshold 0; count rate or a yield of 1 is returned.",
      call. = FALSE
      )
    }

# Lambda parameter
  lambda <- (2 * mean_PHD^2) / (SD_PHD^2 + mean_PHD)

# Probability parameter
  prob <- (SD_PHD^2 - mean_PHD) / (SD_PHD^2 + mean_PHD)

  if (prob < 0 | prob >= 1) {
    stop("Unable to calculate probability.", call. = FALSE)
  } else {
    prob
  }

  Y <- polyaAeppli::pPolyaAeppli(thr_PHD, lambda, prob, lower.tail = FALSE)

  if (output == "Y") Y else if (output == "ct") x / Y
  }

#' @rdname  cor_yield
#'
#' @export
cor_DT <- function(x, deadtime)  x / (1 - (x * deadtime * 10^-9))

#-------------------------------------------------------------------------------
# Supporting functions (arguments)
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

# if variables are not supplied, default to standard variables
inject_args <- function(IC, args, type, check = TRUE) {
  # execute extraction procedure
  purrr::imap(args, ~extract_defaults(.x, .y, point::names_cameca, type, IC))
}

# supply custom variables
supply_args <- function(args, IC, inverse = FALSE) {
  vc_nms <- c(.bl_t = "tc.mt",  .det = "det_type.mt",  .M_PHD = "M_PHD.mt",
              .SD_PHD = "SD_PHD.mt")
  # not NULL
  arg_lgl <- purrr::map_lgl(args, ~{!is_symbol(quo_get_expr(.x))})
  args <- args[arg_lgl]
  # not a symbol
  arg_lgl <- purrr::map_lgl(args, ~{!is.null(quo_get_expr(.x))})
  args <- args[arg_lgl]
  # one of the reference vector
  args <- args[names(args) %in% names(vc_nms)]
  names(args) <- vc_nms[names(vc_nms) %in% names(args)]
  # only valid expression
  mutate(IC, !!!args)
}

# extract default IC variables
extract_defaults <- function(x, y, vars, type, IC) {
  if (!is_symbol(quo_get_expr(x))) {
    # reference table
    vars <- filter(
      vars,
      .data$argument == y,
      .data$use %in% type,
      .data$point %in% colnames(IC)
      )
    if (nrow(vars) != 0) {
      pull(vars, .data$point) %>%
        # caller environment would be similar to enquo()
        parse_quo(env = caller_env())
    } else {
      quo(NULL)
    }
  } else {
    x
  }
}
#-------------------------------------------------------------------------------
# Supporting functions (checks)
#-------------------------------------------------------------------------------

det_check <- function(det) {
  # custom detector type
  if (is.null(quo_get_expr(det))) {
      stop(
        "Supply a variable or character string for the detector type to `.det`",
        call. = FALSE
      )
    }
}

cor_check <- function(parm, args, type) {
  if (is.null(parm)) return(FALSE)
  ref <- list()
  ref$PHD <- c(".M_PHD", ".SD_PHD")
  arg_lgl <- !purrr::map_lgl(args[ref[[type]]], ~is.null(quo_get_expr((.x))))
  if (any(arg_lgl)) {
    if (all(arg_lgl)) {
      TRUE
    } else {
      stop(
        purrr::lift_dl(paste)(list2(
          type,
          "correction requires",
          !!! ref[[type]]
        )),
        call. = FALSE
      )
      }
    } else {
      FALSE
    }
}

# check if arguments exist in data frame
argument_check <- function(IC, args, type) {
  # not NULL
  arg_lgl <- purrr::map_lgl(args, ~{is.null(quo_get_expr(.x))})
  pos_args <- args[!arg_lgl]
  if (any(!sapply(pos_args, as_name) %in% colnames(IC))) {
    stop("Tibble does not contain the supplied variables!", call. = FALSE)
  }
  neg_args <- args[arg_lgl]
  # default names for missing argument (but not meta data)
  vars <- filter(
    point::names_cameca, .data$use %in% type,
    .data$argument %in% names(neg_args)
    ) %>%
    pull(.data$point)
  if (any(!vars %in% colnames(IC))) {
    stop("Tibble does not contain the default variables!", call. = FALSE)
  }
}
