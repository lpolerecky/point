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
#' \code{vignette("IC-diagnostics", package = "point")} for examples. The
#' argument \code{.return = "augmented"} can be used to toggle between
#' returning the augmented versions (with the original data labelled as
#' execution 1) of the datasets and the results \code{.return = "results"} for
#' outlier detection as well as associated statistics of the selected procedure.
#'
#' @param .df A tibble containing processed ion count data.
#' @param .ion1 A character string constituting the heavy isotope ("13C").
#' @param .ion2 A character string constituting the light isotope ("12C").
#' @param ... Variables for grouping.
#' @param .method Character string for the method for diagnostics (default =
#' \code{"CooksD"}, see details).
#' @param .reps Numeric setting the number of repeated iterations of outlier
#' detection (default = 1).
#' @param .Xt A variable constituting the ion count rate (defaults to
#' variables generated with \code{read_IC()})
#' @param .N A variable constituting the ion counts (defaults to variables
#' generated with \code{read_IC()}.).
#' @param .species A variable constituting the species analysed (defaults to
#' variables generated with \code{read_IC()}).
#' @param .t A variable constituting the time of the analyses (defaults to
#' variables generated with \code{read_IC()}).
#' @param .output Can be set to \code{"complete"} which returns \code{stat_R()}
#' statistics and outlier detection results following the selected method
#' (see above argument \code{.method}).
#' @param .hyp Hypothesis test. Only usable in combination with a selection of
#' methods (see above argument \code{.method}) and details.
#' @param .return Either the augmented data (\code{"augmented"}) or statistic
#' and outlier detection results (\code{"results"}).
#' @param .meta Logical whether to preserve the metadata as an attribute
#' (defaults to TRUE).
#' @param .plot Logical indicating whether plot is generated.
#' @param .plot_type Character string determining whether the returned a
#' \code{"static"} \code{\link[ggplot2::ggplot2]{ggplot2}} or an
#' \code{"interactive"} plot with \code{\link[plotly::ggplotly]{ggplotly}}.
#' @param .plot_stat Adds a statistic label to the plot (e.g. . \code{"M"}), see
#' \code{point::nm_stat_R} for the full selection of statistics available.
#' @param .plot_iso A character string (e.g. \code{"VPDB"}) for the delta
#' conversion of R \code{?calib_R()} for options.
#'
#' @return A dynamic static plot are returned together with a
#' \code{\link[tibble::tibble]{tibble}} containing \code{stat_R()} statistics
#' and diagnostics associated with the chosen method.
#'
#' @export
#' @examples
#' # Modelled ion count dataset
#'
#' # Cook's D style diagnostic-augmentation of ion count data for
#' # isotope ratios; 3 repeats
#'
#'
diag_R <- function(.df, .ion1, .ion2, ..., .method = "CooksD", .reps = 1,
                   .Xt = Xt.pr, .N = N.pr, .species = species.nm, .t = t.nm,
                   .output = "complete", .hyp = "none", .return = "results",
                   .meta = TRUE, .plot = FALSE, .plot_type = "static",
                   .plot_stat = NULL, .plot_iso = FALSE){

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)

  # stat_R variables
  args <- enquos(.Xt = .Xt, .N = .N, .species = .species, .t = .t)

  # Metadata
  if(.meta) meta <- unfold(xc, merge = FALSE)

  # Repetitions
  if (.method != "IR") {
  max <- .reps + 1

  # plot stats
  if (.output == "flag" & !is.null(.plot_stat)) {
    .plot_stat <- NULL
    warning("If argument `.output = flag`, argument `plot_stat` defaults to
            `NULL`")
  }

  # Empty list for iteration storage
  ls_tb <- rep_named(as.character(1:max), list2())
  ls_tb[[1]] <- list2(df = .df, results = NULL)

  # Execute repeated cycles of augmentation
  ls_tb <- ls_tb  %>%
    purrr::accumulate(
      rerun_diag_R,
      .ion1 = .ion1,
      .ion2 = .ion2,
      !!! gr_by,
      .method = .method,
      .Xt = !!args[[".Xt"]],
      .N = !!args[[".N"]],
      .species = !!args[[".species"]],
      .t = !!args[[".t"]],
      .output = .output
      )

  if (.plot & .return != "results") {
    .return <- "results"
    warning("If `plot` is `TRUE`, `return` defaults to `results`")
  }


  if (.return == "augmented") {
    df <-reduce_diag(ls_tb, type = "df") %>%
      mutate(execution = as.numeric(execution))
    } else {
      df <- reduce_diag(ls_tb, type = "results") %>%
        mutate(execution = as.numeric(execution) - 1)
    }


  # metadata
  if (.meta) df <- fold(df, type = ".mt",  meta = meta)

  if (.plot) {
    df %T>%
      {print(
        gg_IC(., .ion1 = .ion1, .ion2 = .ion2, .method = .method, .plot_type = .plot_type, !!!gr_by,
              .labels = .plot_stat)
        )
      }
    }
  return(df)
  } else {

    # Auto-correlation
    df <- diag_R_exec(.df, .ion1 = .ion1, .ion2 = .ion2, !!! gr_by,
                      .method = "IR", .Xt = !!args[[".Xt"]], .N = !!args[[".N"]],
                      .species = !!args[[".species"]], .t = !!args[[".t"]])
    if (.plot) {

    df %T>%
      {print(gg_IR(., .lag = lag, .acf = acf, .flag = flag, !!! gr_by,
                   .hat = 0, .error = hat_e_acf))
      }
  }
  return(df)
}
}


rerun_diag_R <- function(out, input, .ion1, .ion2, ..., .method,.Xt = Xt.pr,
                         .N = N.pr, .species = species.nm, .t = t.nm,
                         .output){

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)

  # remove white space in ion names
  .ion1 <- str_replace_all(.ion1, "\\s", "")
  .ion2 <- str_replace_all(.ion2, "\\s", "")

  # stat_R variables
  args <- enquos(.Xt = .Xt, .N = .N, .species = .species, .t = .t)

  # Heavy isotope
  Xt1 <- quo_updt(args[[".Xt"]], post = .ion1) # count rate
  N1 <- quo_updt(args[[".N"]], post = .ion1) # counts

  # Ligh isotope
  Xt2 <- quo_updt(args[[".Xt"]], post = .ion2) # count rate
  N2  <- quo_updt(args[[".N"]], post = .ion2) # counts

  # Execute
  out <- diag_R_exec(
    out$df,
    .ion1,
    .ion2,
    !!! gr_by,
    .method = .method,
    .Xt = !!args[[".Xt"]],
    .N = !!args[[".N"]],
    .species = !!args[[".species"]],
    .t = !!args[[".t"]]
    )

  # Save augmented dataframe for next cycle
  df_aug <- select(
    filter(out, flag == "confluent"),
    !!!gr_by, !!args[[".t"]], !!Xt1, !!Xt2, !!N1, !!N2
    ) %>%
    pivot_longer(
      -c(!!!gr_by, !!args[[".t"]]),
      names_to = c(".value", as_name(args[[".species"]])),
      names_sep = "[[:punct:]]{1}(?=[[:digit:]]{1,2})"
      )

  if (.output == "flag") {
    out <- select(out, -ends_with(paste0("R_", as_name(args[[".Xt"]]))))
  }
  # Output
  out <- list2(df = df_aug, results = out)

  return(out)

}


diag_R_exec <- function(.df, .ion1, .ion2, ..., .method, .Xt = Xt.pr,
                        .N = N.pr, .species = species.nm, .t = t.nm){

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)

  # stat_R variables
  args <- enquos(Xt = .Xt, N = .N, species = .species, t = .t)

  # Diagnostic call
  diag_vc <- c("Cameca", "CooksD", "Rm", "CV", "QQ", "norm_E", "IR")
  diag_method <- purrr::map(
    diag_vc,
    call2,
    expr(.),
    .ion1 = expr(.ion1),
    .ion2 = expr(.ion2),
    !!! gr_by,
    .Xt = expr(!! args[["Xt"]]),
    .t = expr(!! args[["t"]]),
    .output = "complete"
    ) %>%
    set_names(nm = diag_vc)

  # Descriptive an predictive statistics for ion ratios
  stat_R(.df, .ion1, .ion2, !!! gr_by, .Xt = !!args[["Xt"]], .N = !!args[["N"]],
         .species = !!args[["species"]], .t = !!args[["t"]],
         .output = "complete", .zero = TRUE) %>%
    eval_tidy(expr = diag_method[[.method]])

  }



reduce_diag <- function(ls, type){

  # Reduce the results to a single dataframe
  ls %>%
    purrr::transpose() %>%
    purrr::pluck(type) %>%
    bind_rows(.id = "execution")

}
