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
#' \code{point::\link[point:names_plot]{names_plot}}. The
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
#'  the significance of inter-isotope variability.
#' @param .method Character string for the method for diagnostics (default =
#'  \code{"CooksD"}, see details).
#' @param .reps Numeric setting the number of repeated iterations of outlier
#'  detection (default = 1).
#' @param .X A variable constituting the ion count rate (defaults to
#'  variables generated with \code{read_IC()})
#' @param .N A variable constituting the ion counts (defaults to variables
#'  generated with \code{read_IC()}.).
#' @param .species A variable constituting the species analysed (defaults to
#'  variables generated with \code{read_IC()}).
#' @param .t A variable constituting the time of the analyses (defaults to
#'  variables generated with \code{read_IC()}).
#' @param .output Can be set to \code{"complete"} which returns \code{stat_R()}
#'  and \code{stat_X()} statistics, diagnostics, and inference test results
#'  following the selected method (see above argument \code{.method});
#'  \code{"augmented"} for the augmented IC data after diagnostics;
#'  \code{"diagnostic"} returns \code{stat_R()} and \code{stat_X()} statistics
#'  and outlier detection; \code{"outlier"} for outlier detection;
#'  \code{"inference"} for only inference test statistics results
#'  (default = \code{"inference"}).
#' @param .label For printing nice latex labels use \code{"latex"} (default =
#'  \code{NULL}).
#' @param .meta Logical whether to preserve the metadata as an attribute
#'  (defaults to TRUE).
#' @param .alpha_level Significance level of hypothesis test.
#' @param .hyp Hypothesis test appropriate for the selected method
#'  (default = \code{"none"}).
#' @param .plot Logical indicating whether plot is generated.
#' @param .plot_type Character string determining whether the returned plot is
#'  \code{"static"} \code{ggplot2::\link[ggplot2:ggplot2]{ggplot2}()}(currently
#'  only supported option).
#' @param .plot_stat Adds a statistic label to the plot (e.g. . \code{"M"}), see
#'  \code{point::nm_stat_R} for the full selection of statistics available.
#' @param .plot_iso A character string (e.g. \code{"VPDB"}) for the delta
#'  conversion of R (see \code{?calib_R()} for options).
#' @param .plot_outlier_labs A character vector of length two for the colourbar
#'  text for outliers (default = c("divergent", "confluent")).
#' @param .mc_cores Number of workers for parallel execution (Does not work on
#'   Windows).
#'
#' @return A \code{ggplot2::\link[ggplot2:ggplot]{ggplot}()} is returned
#'  (if \code{.plot = TRUE}) along with a
#'  \code{tibble::\link[tibble:tibble]{tibble}()} which can contain statistics
#'  diagnostics, hypothesis test results associated with the chosen method and
#'  depending on the argument \code{.output}.
#'
#' @export
#' @examples
#' # Modelled ion count dataset
#' # Cook's D style diagnostic-augmentation of ion count data for
#' # isotope ratios
#' diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .plot = TRUE)
#'
diag_R <- function(.IC, .ion1, .ion2, ..., .nest = NULL, .method = "CooksD",
                   .reps = 1, .X = NULL, .N = NULL, .species = NULL,
                   .t = NULL, .output = "inference", .label = "none",
                   .meta = FALSE, .alpha_level = 0.05, .hyp = "none",
                   .plot = FALSE, .plot_type = "static", .plot_stat = NULL,
                   .plot_iso = FALSE,
                   .plot_outlier_labs = c("divergent", "confluent"),
                   .mc_cores = 1){

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)

  # stat_R variables
  args <- inject_args(
    .IC,
    enquos(.X = .X, .N = .N, .species = .species, .t = .t),
    type = c("processed", "group")
  )

  # Argument check
  argument_check(.IC, args, "processed")

  # diagnostics variables
  diag_args <- rlang::list2(
    .method = .method,
    .alpha_level = .alpha_level,
    .hyp = .hyp
  )

  # plot variables
  plot_args <- rlang::list2(
    .method = .method,
    .plot_type = .plot_type,
    .plot_stat = .plot_stat,
    .alpha_level = .alpha_level,
    .plot_outlier_labs = .plot_outlier_labs
  )

  # warnings
  if (.plot & .output == "augmented") {
     warning("If `plot` is `TRUE`, `.output` cannot be `augmented`.",
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
  ls_tb <- rlang::rep_named(as.character(1:max), rlang::list2())
  ls_tb[[1]] <- rlang::list2(IC = .IC, results = NULL)

  # Repeated cycles of augmentation call
  if (.method != "IR") {
    diag_call <- rlang::call2(
      .fn = "accumulate",
      .x = ls_tb,
      .f = rerun_diag_R,
      .ion1 = .ion1,
      .ion2 = .ion2,
      !!! gr_by,
      .args = args,
      .diag_args = diag_args,
      .output = .output,
      .mc_cores = .mc_cores,
      .ns = "purrr"
    )
  } else {
    # Single cycle for auto-correlation (IR) call
    diag_call <- rlang::call2(
      "diag_R_exec",
      .IC,
      .ion1 = .ion1,
      .ion2 = .ion2,
      !!! gr_by,
      .args = args,
      .diag_args = diag_args,
      .mc_cores = .mc_cores
    )
  }

  # Execute
  IC <- eval(diag_call)

  # Collapse repeated analysis
  if (.method != "IR") IC <- reduce_diag(IC, .output)

  # Plot data
  if (.plot) {
    if (.method != "IR") {
      plot_call <- rlang::call2(
        "gg_IC",
        IC,
        .ion1 = .ion1,
        .ion2 = .ion2,
        !!! gr_by,
        !!! args[!names(args) %in% c(".species", ".t")],
        .plot_args = plot_args
      )
    } else {
      plot_call <- rlang::call2(
        "gg_IR",
        IC,
        .lag = rlang::parse_expr("lag"),
        .acf = rlang::parse_expr("acf"),
        .flag = rlang::parse_expr("flag"),
        !!! gr_by,
        .sd = rlang::parse_expr("e_acf")
      )
    }
    print(eval(plot_call))
  }

  # Inferences
  if (.output == "inference" | .output == "complete") {

    if (.method %in%
        dplyr::filter(point::names_plot, .data$inference == "eval_diag")$name) {

      eval_call <- rlang::call2(
        "eval_diag",
        IC,
        .ion1 = .ion1,
        .ion2 = .ion2,
        !!! gr_by,
        .nest = enquo(.nest),
        !!! args,
        .flag = rlang::parse_expr("flag"),
        .output = .output,
        .label = .label,
        .mc_cores = .mc_cores
      )

      IC <- eval(eval_call)
    }

    if (.method %in%
        dplyr::filter(point::names_plot, .data$inference == "external")$name) {
      IC <- dplyr::distinct(IC, !!! gr_by, .data$hyp)
    }
  }

  # Return metadata
  if (.meta) IC <- fold(IC, type = ".mt",  meta = meta)
  IC
}

# rerunning the diagnostics for several iterations
rerun_diag_R <- function(out, input, .ion1, .ion2, ..., .args, .diag_args,
                         .output, .mc_cores){

  # Grouping
  gr_by <- enquos(...)

  # Remove white space in ion names
  .ion1 <- ion_trim(.ion1)
  .ion2 <- ion_trim(.ion2)

  # Rare isotope
  X1 <- quo_updt(.args[[".X"]], post = .ion1) # count rate
  N1 <- quo_updt(.args[[".N"]], post = .ion1) # counts

  # Common isotope
  X2 <- quo_updt(.args[[".X"]], post = .ion2) # count rate
  N2 <- quo_updt(.args[[".N"]], post = .ion2) # counts

  # Execute
  out <- diag_R_exec(
    out$IC,
    .ion1,
    .ion2,
    !!! gr_by,
    .args = .args,
    .diag_args = .diag_args,
    .mc_cores = .mc_cores
  )

  # Save augmented data frame for next cycle
  aug <- dplyr::select(
    dplyr::filter(out, .data$flag == "confluent"),
    !!! gr_by, !! .args[[".t"]], !! X1, !! X2, !! N1, !! N2
  ) |>
    tidyr::pivot_longer(
      -c(!!!gr_by, !! .args[[".t"]]),
      names_to = c(".value", as_name(.args[[".species"]])),
      names_sep = "\\.(?=[[:digit:]]{1,2})"
    )

  # Filter all stat variables out, except Poisson prediction for rare isotope
  if (.output == "outlier") {
    except <- quo_updt(.args[[".N"]], pre = "hat_S", post = .ion1)  |>
      as_name()
    out <- dplyr::select(
      out,
      -dplyr::any_of(all_args(.args, .ion1, .ion2, except))
    )
  }
  # Output
  rlang::list2(IC = aug, results = out)
}

# execute
diag_R_exec <- function(.IC, .ion1, .ion2, ..., .args, .diag_args, .mc_cores){

  # Grouping
  gr_by <- enquos(...)

  # Diagnostic call check
  if(!.diag_args[".method"] %in% point::names_plot$name) {
    stop("Method does not exist.", call. = FALSE)
  }

  # Diagnostics (single call)
  diag_call <- rlang::call2(
    .diag_args[[".method"]],
    rlang::expr(tb_R),
    .ion1 = .ion1,
    .ion2 = .ion2,
    !!! gr_by,
    !!! .args,
    .output = "complete",
    !!! .diag_args[names(.diag_args) != ".method"],
    .mc_cores = .mc_cores
  )

  # Descriptive an predictive calls for single ions statistics
  X_call <- rlang::call2(
    "stat_X",
    rlang::expr(.IC),
    !!! gr_by,
    !!! .args,
    .output = "complete"
  )
  # Descriptive an predictive calls for ion ratios statistics
  R_call <- rlang::call2(
    "stat_R",
    rlang::expr(tb_X),
    .ion1 = .ion1,
    .ion2 = .ion2,
    !!! gr_by,
    !!! .args,
    .output = "complete",
    .zero = TRUE
  )

  # Execute
  tb_X <- rlang::eval_tidy(expr = X_call)
  tb_R <- rlang::eval_tidy(expr = R_call)
  rlang::eval_tidy(expr = diag_call)
}

# Reduce the results to a single data frame
reduce_diag <- function(ls, output){
  purrr::transpose(ls) |>
    purrr::pluck(if (output == "augmented") "IC" else "results") |>
    dplyr::bind_rows(.id = "execution") |>
    dplyr::mutate(
      execution =
        as.numeric(.data$execution) - if (output == "augmented") 0 else 1
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
  args_R <- rlang::list2(
    !!! arg_builder(args, "R"),
    R = quo_updt(args[[".X"]], pre = "R"),
    ratio = rlang::parse_quo("ratio.nm", env = rlang::quo_get_env(args[[".X"]]))
  )
  args <- append(args_Xt, args_R)
  if (isTRUE(chr)) {
    vc_args <- sapply(args, as_name)
    vc_args[!vc_args %in% except]
  } else {
    args
  }
}
