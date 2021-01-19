#' Diagnostics isotope count data
#'
#' \code{diag_R} wrapper function for diagnostics on isotope count data
#'
#' The \code{diag_R} function performs an internal call to stat_R to perform
#' diagnostics on the influence of individual measurements on the blockwise or
#' global (i.e., a complete analysis) statistics. It identifies potentially
#' influential measurements that indicate heterogeneity in the analytical
#' substrate. See \code{vignette("IC-diagnostics", package = "point")} for more
#' information on how to use the function, and possible methods for diagnostics.
#'
#' @param df A tibble containing processed ion count data.
#' @param method Character string for the type of diagnostics. "Cameca" pertains
#' to the default Cameca software setting. "CooksD" pertains to global
#' regression diagnostics based on Cook's D statistics.
#' @param args A list of quosures pertaining to the variables required for a
#' call to stat_R. The function expr_R can be used to simplify setting this
#' argument.
#' @param reps The number of iterations of diagnostic-augmentation cycles.
#' @param ... Variables for grouping.
#' @param plot Logical indicating whether plot is produced, if FALSE a list of
#' tibbles with length `reps` is produced.
#' @param plot_type Character string determining whether the returned a
#' \code{"static"} plot or an \code{"interactive"} plot with
#' \code{\link[plotly:ggplotly]{ggplotly}}.
#' @param iso Logical indicating whether the R value should be converted to the
#' delta scale.
#' @param isoscale If the argument \code{iso} is set to TRUE, a character string
#' (e.g. \code{"VPDB"}) for the delta conversion should be given. See
#' \code{?calib_R()} for options.
#'
#' @return A dynamic plot or static plot. If the argument plot is \code{FALSE},
#' a list containing several \code{\link[tibble:tibble]{tibble}} is returned.
#' The elements named df contain the original dataset and augmented versions
#' from element position number 2 onward. The elements named results contain
#' the flag variable for outlier detection as well as associated statistics of
#' the selected procedure.
#'
#' @export
#' @examples
#' # Modelled ion count dataset
#'
#' # Cook's D style diagnostic-augmentation of ion count data for
#' # isotope ratios; 3 repeats
#' diag_R(sim_IC_extremes, "13C", "12C", simulation, trend, iso_offset, .output = "complete)
#'
diag_R <- function(.df, .ion1, .ion2, ..., .method = "CooksD", .reps = 1,
                   .Xt = Xt.pr, .N = N.pr, .species = species.nm, .t = t.nm,
                   .output = "complete", .hyp = "none", .return = "results",
                   .plot = FALSE, .plot_type = "static", .plot_stat = NULL,
                   .plot_iso = FALSE
                   ){

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)

  # stat_R variables
  args <- enquos(.Xt = .Xt, .N = .N, .species = .species, .t = .t)

  # Repetitions
  if (.method != "IR") {
  max <- .reps + 1

  # plot stats
  if (.output == "flag" & !is.null(.plot_stat)) {
    .plot_stat <- NULL
    warning("If argument `.output = flag`, argument `plot_stat defaults to `NULL`")
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

  if (.plot) {
    df %T>%
      {print(
        gg_IC(., .ion1 = .ion1, .ion2 = .ion2, .method = .method, !!!gr_by,
              .labels = .plot_stat, .plot_type = .plot_type)
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
    filter(out, flag == "good"),
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
         .output = "complete",
         .zero = TRUE) %>%
    eval_tidy(expr = diag_method[[.method]])

  }



reduce_diag <- function(ls, type){

  # Reduce the results to a single dataframe
  ls %>%
    purrr::transpose() %>%
    purrr::pluck(type) %>%
    bind_rows(.id = "execution")

}
