#' Diagnostics isotope count data
#'
#' \code{diag_R} wrapper function for diagnostics on isotope count data
#'
#' The \code{diag_R} function performs an internal call to stat_R to perform
#' diagnostics on the influence of individual measurements on the block-wise or
#' global (i.e., a complete analysis) statistics. It identifies potentially
#' influential measurements that indicate heterogeneity in the analytical
#' substrate or measurement. See
#' \code{vignette("IC-diagnostics", package = "point")} for more information on
#' how to use the function, and possible methods for outlier detection. Options
#' are \code{"CooksD"} (default), \code{"Cameca"}, \code{"Rm"}, \code{"norm_E"},
#' \code{"CV"}, \code{"IR"}, and \code{"QQ"}, see the
#' \code{vignette("IC-diagnostics", package = "point")} for examples and
#' \code{point::\link[point:names_diag]{names_diag}}. The
#' argument \code{.output} can be used to toggle between \code{"complete"};
#' returning \code{stat_R()} and \code{stat_X()} statistics, diagnostics, and
#' inference test results, \code{"augmented"}; returning the augmented IC after
#' removing outliers, \code{"diagnostic"}; for only outlier detection results;
#' \code{"diagnostic"}; for statistics and outlier detection, or
#' \code{"inference"}; returns only inference test statistics results
#' (default = inference).
#'
#' @param .IC A tibble containing processed ion count data.
#' @param .ion1 A character string constituting the rare isotope ("13C").
#' @param .ion2 A character string constituting the common isotope ("12C").
#' @param ... Variables for grouping.
#' @param .nest A variable hat identifies a series of analyses to calculate
#' the significance of inter-isotope variability.
#' @param .method Character string for the method for diagnostics (default =
#' \code{"CooksD"}, see details).
#' @param .reps Numeric setting the number of repeated iterations of outlier
#' detection (default = 1).
#' @param .X A variable constituting the ion count rate (defaults to
#' variables generated with \code{read_IC()})
#' @param .N A variable constituting the ion counts (defaults to variables
#' generated with \code{read_IC()}.).
#' @param .species A variable constituting the species analysed (defaults to
#' variables generated with \code{read_IC()}).
#' @param .t A variable constituting the time of the analyses (defaults to
#' variables generated with \code{read_IC()}).
#' @param .output Can be set to \code{"complete"} which returns \code{stat_R()}
#' and \code{stat_X()} statistics, diagnostics, and inference test results
#' following the selected method (see above argument \code{.method});
#' \code{"augmented"} for the augmented IC data after diagnostics;
#' \code{"diagnostic"} returns \code{stat_R()} and \code{stat_X()} statistics
#' and outlier detection; \code{"outlier"} for outlier detection;
#' \code{"inference"} for only inference test statistics results
#' (default = \code{"inference"}).
#' @param .label For printing nice latex labels use \code{"latex"} (default =
#' \code{NULL}).
#' @param .meta Logical whether to preserve the metadata as an attribute
#' (defaults to TRUE).
#' @param .alpha_level Significance level of hypothesis test.
#' @param .hyp Hypothesis test appropriate for the selected method
#' (default = \code{"none"}).
#' @param .plot Logical indicating whether plot is generated.
#' @param .plot_type Character string determining whether the returned plot is
#' \code{"static"} \code{ggplot2::\link[ggplot2:ggplot2]{ggplot2}()}(currently
#' only supported option).
#' @param .plot_stat Adds a statistic label to the plot (e.g. . \code{"M"}), see
#' \code{point::nm_stat_R} for the full selection of statistics available.
#' @param .plot_iso A character string (e.g. \code{"VPDB"}) for the delta
#' conversion of R (see \code{?calib_R()} for options).
#' @param .plot_outlier_labs A character vector of length two for the colourbar
#' text for outliers (default = c("divergent", "confluent")).
#'
#' @return A \code{ggplot2::\link[ggplot2:ggplot]{ggplot}()} is returned
#' (if \code{.plot = TRUE}) along with a
#' \code{tibble::\link[tibble:tibble]{tibble}()} which can contain statistics
#' diagnostics, hypothesis test results associated with the chosen method and
#' depending on the argument \code{.output}.
#'
#' @export
#' @examples
#' # Modelled ion count dataset
#' # Cook's D style diagnostic-augmentation of ion count data for
#' # isotope ratios
#' diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .plot = TRUE)
#'
diag_R <- function(.IC, .ion1, .ion2, ..., .nest = NULL, .method = "CooksD",
                   .reps = 1, .X = Xt.pr, .N = N.pr, .species = species.nm,
                   .t = t.nm, .output = "inference", .label = NULL,
                   .meta = TRUE, .alpha_level = 0.05, .hyp = "none",
                   .plot = FALSE, .plot_type = "static", .plot_stat = NULL,
                   .plot_iso = FALSE,
                   .plot_outlier_labs = c("divergent", "confluent")){

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)

  # stat_R variables
  args <- enquos(.X = .X, .N = .N, .species = .species, .t = .t)
  # diagnostics variables
  diag_args <- list2(.method = .method, .alpha_level = .alpha_level,
                     .hyp = .hyp)
  # plot variables
  plot_args <- list2(.method = .method, .plot_type = .plot_type,
                     .labels = .plot_stat, .alpha_level = .alpha_level,
                     .plot_outlier_labs = .plot_outlier_labs)
  # warnings
  if (.plot & .output == "augmented") {
     warning("If `plot` is `TRUE`, `.output` cannot be `augmented`",
             call. = FALSE)
  }

  # Metadata
  if(.meta) meta <- unfold(.IC, merge = FALSE)

  # Set alternative hypothesis tests
  if (.output == "inference" | .output == "complete") {
    if (.method == "QQ") diag_args[[".hyp"]] <- "norm"
    if (.method == "IR") diag_args[[".hyp"]] <- "ljung"
    if (.method == "CV") diag_args[[".hyp"]] <- "bp"
  }

  # Repetitions
  max <- .reps + 1
  # Empty list for iteration storage
  ls_tb <- rep_named(as.character(1:max), list2())
  ls_tb[[1]] <- list2(IC = .IC, results = NULL)

  # Repeated cycles of augmentation call
  if (.method != "IR") {
    data_env <- env(data = ls_tb)
    diag_call <- call2(.fn = "accumulate", .x = ls_tb,
                       .f = rerun_diag_R, .ion1 = .ion1, .ion2 = .ion2,
                       !!! gr_by, !!! args, !!! diag_args, .output = .output,
                       .ns = "purrr")
  }
  # Single cycle for auto-correlation (IR) call
  if (.method == "IR") {
    data_env <- env(data = .IC)
    diag_call <- call2("diag_R_exec", .IC, .ion1 = .ion1, .ion2 = .ion2,
                       !!! gr_by, !!! args, !!! diag_args)
    }
  # Execute
  IC <- eval(diag_call, data_env)

  # Collapse repeated analysis
  if (.method != "IR") IC <- reduce_diag(IC, .output)

  # Return metadata
  if (.meta & !is.null(meta)) IC <- fold(IC, type = ".mt",  meta = meta)

  # Plot data
  if (.plot) {
    data_env <- env(data = IC)
    if (.method != "IR") {
      plot_call <- call2("gg_IC", IC, .ion1 = .ion1, .ion2 = .ion2, !!! gr_by,
                          !!! args, !!! plot_args)
    }
    if (.method == "IR") {
      plot_call <- call2("gg_IR", IC, .lag = parse_quo("lag", env = data_env),
                         .acf = parse_quo("acf", env = data_env),
                         .flag = parse_quo("flag", env = data_env), !!! gr_by,
                         .sd = parse_quo("e_acf", env = data_env))
    }
    print(eval(plot_call, data_env))
  }

  # Inferences
  if (.output == "inference" | .output == "complete") {
    if (.method %in%
        filter(point::names_diag, .data$inference == "eval_diag")$name) {
      data_env <- env(data = IC)
      eval_call <- call2("eval_diag", IC, .ion1 = .ion1, .ion2 = .ion2,
                          !!! gr_by, .nest = enquo(.nest), !!! args,
                         .flag = parse_quo("flag", env = data_env),
                         .output = .output, .label = .label)
      return(eval(eval_call, data_env))
      }
    if (.method %in%
        filter(point::names_diag, .data$inference == "external")$name) {
      distinct(IC, !!!gr_by, .data$hyp)
      }
    } else {
      return(IC)
      }
}


rerun_diag_R <- function(out, input, .ion1, .ion2, ..., .method, .X = Xt.pr,
                         .N = N.pr, .species = species.nm, .t = t.nm,
                         .output, .hyp, .alpha_level){

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)

  # Remove white space in ion names
  .ion1 <- ion_trim(.ion1)
  .ion2 <- ion_trim(.ion2)

  # stat_R variables
  args <- enquos(.X = .X, .N = .N, .species = .species, .t = .t)
  # diagnostics variables
  diag_args <- list2(.method = .method, .hyp = .hyp,
                     .alpha_level = .alpha_level)

  # Rare isotope
  X1 <- quo_updt(args[[".X"]], post = .ion1) # count rate
  N1 <- quo_updt(args[[".N"]], post = .ion1) # counts

  # Common isotope
  X2 <- quo_updt(args[[".X"]], post = .ion2) # count rate
  N2  <- quo_updt(args[[".N"]], post = .ion2) # counts

  # Execute
  out <- diag_R_exec(
    out$IC,
    .ion1,
    .ion2,
    !!! gr_by,
    .method = .method,
    .X = !!args[[".X"]],
    .N = !!args[[".N"]],
    .species = !!args[[".species"]],
    .t = !!args[[".t"]],
    .hyp = .hyp,
    .alpha_level = .alpha_level
    )

  # Save augmented data frame for next cycle
  aug <- select(
    filter(out, .data$flag == "confluent"),
    !!! gr_by, !! args[[".t"]], !! X1, !! X2, !! N1, !! N2
    ) %>%
    tidyr::pivot_longer(
      -c(!!!gr_by, !!args[[".t"]]),
      names_to = c(".value", as_name(args[[".species"]])),
      names_sep = "\\.(?=[[:digit:]]{1,2})"
      )

  # Filter all stat variables out, except Poisson prediction for rare isotope
  if (.output == "outlier") {
    except <- quo_updt(args[[".N"]], pre = "hat_S", post = .ion1) %>% as_name()
    out <- select(out, -any_of(all_args(args, .ion1, .ion2, except)))
    }
  # Output
  out <- list2(IC = aug, results = out)
  return(out)
}

diag_R_exec <- function(.IC, .ion1, .ion2, ..., .method, .X = Xt.pr, .N = N.pr,
                        .species = species.nm, .t = t.nm, .hyp, .alpha_level){

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)

  # stat_R variables
  args <- enquos(.X = .X, .N = .N, .species = .species, .t = .t)

  # Diagnostic call check
  if(!.method %in% point::names_diag$name) {
    stop("Method does not exist", call. = FALSE)
  }

  diag_call <- call2(.method, expr(.), .ion1 = .ion1, .ion2 = .ion2,
                     !!! gr_by, .X = expr(!! args[[".X"]]),
                     .t = expr(!! args[[".t"]]), .output = "complete",
                     .hyp = .hyp, .alpha_level = .alpha_level)

  # Descriptive an predictive calls for single ions and ion ratios stats
  X_call <- call2("stat_X", expr(.IC), !!! gr_by, !!! args,
                  .output = "complete")
  R_call <- call2("stat_R", expr(.), .ion1 = .ion1, .ion2 = .ion2, !!! gr_by,
                  !!! args, .output = "complete", .zero = TRUE)

  # Execute
  data_env <- env(data = .IC)
  tb_R <- eval(X_call, data_env) %>%
    eval_tidy(expr = R_call) %>%
    eval_tidy(expr = diag_call)
  }


# Reduce the results to a single data frame
reduce_diag <- function(ls, output){
  purrr::transpose(ls) %>%
    purrr::pluck(ifelse(output == "augmented", "IC", "results")) %>%
    bind_rows(.id = "execution")  %>%
    mutate(
      execution =
        as.numeric(.data$execution) - ifelse(output == "augmented", 0, 1)
      )
}

# Select arguments
all_args <- function(args, ion1, ion2, except = NULL, chr = TRUE) {
  args_Xt <- purrr::flatten(
    purrr::map2(
      c(ion1, ion2),
      seq_along(c(ion1, ion2)),
      ~arg_builder(args, "X", .x, .y)
      )
    )
  args_R <- list2(
    !!! arg_builder(args, "R"),
    R = quo_updt(args[[".X"]], pre = "R"),
    ratio = parse_quo("ratio.nm", env = quo_get_env(args[[".X"]]))
    )
  args <- append(args_Xt, args_R)
  if (isTRUE(chr)) {
    vc_args <- sapply(args, as_name)
    return(vc_args[!vc_args %in% except])
    } else {
      return(args)
      }
}
