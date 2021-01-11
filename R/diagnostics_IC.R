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
#' # Create expression for mimicking R_stat() call
#' expr_R_stat <- expr_R(Xt = "Xt.sim",
#'                       N = "N.sim",
#'                       species = "species",
#'                       ion1 = "13C",
#'                       ion2 = "12C"
#'                       )
#'
#' # Cook's D style diagnostic-augmentation of ion count data for
#' # isotope ratios; 3 repeats
#' ls.dia <- diag_R(sim_IC_extremes,
#'                  method = "CooksD",
#'                  args = expr_R_stat,
#'                  reps = 2,
#'                  simulation,
#'                  trend,
#'                  iso_offset,
#'                  output = "complete",
#'                  plot = FALSE
#'                  )
diag_R <- function(df,
                   method = "Cameca",
                   args = expr_R(NULL),
                   reps = 1,
                   ...,
                   output = "flag",
                   plot = TRUE,
                   plot_type = "interactive",
                   iso = TRUE,
                   isoscale = NULL
                   ){

  # Check if output is consistent with plot call
  if (output == "flag" & plot == TRUE){
      stop("argument plot = TRUE requires argument output to be set to complete",
           call. = FALSE
           )
    }

  gr_by <- enquos(...)

# repetitions
  max <- reps + 1

# empty list for iteration storage
  ls.tb <- rlang::rep_named(as.character(1:max), lst())

# Remove zeros
  if (method != "Cameca"){

    df <- zeroCt(
      df,
      !! args[["N"]],
      !! args[["species"]],
      as_name(args[["ion1"]]),
      as_name(args[["ion2"]]),
      !!! gr_by,
      warn = FALSE
                 )
  }

# ID for connecting flag to original dataframe
  df <- ID_builder(df, !! args[["species"]], !!! gr_by)
# set initial dataset
  df <- filter(
    df,
    !! args[["species"]] == !! args[["ion1"]] |
    !! args[["species"]] == !! args[["ion2"]]
               )
  ls.tb[[1]] <- lst(df = df, results = NULL)

# Execute repeated cycles of augmentation
  ls.tb <- ls.tb  %>%
    purrr::accumulate(rerun_diag_R,
                      method = method,
                      args = args,
                      !!! gr_by,
                      output = output
                      )
  if (plot) {

    ls.tb %T>%
      {print(plot_diag_R(ls_df = .,
                         args = args,
                         !!! gr_by,
                         plot_title = method,
                         type = plot_type,
                         iso = iso,
                         isoscale = isoscale
      )
      )
      }

    } else {

      return(ls.tb)

      }

}

#' @export
rerun_diag_R <- function(out, input, method, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

  # variables
  ls.vars <- var_fun(out$df, gr_by, args = args)

  # if called this way then output is set fixed (more flexible use with diag_R_exec)
  out <- diag_R_exec(
    out$df,
    method = method,
    args = args,
    !!! gr_by,
    output = output
    )

  # save augmented dataframe for next cycle
  df.aug <- filter(out, flag == "good")

  if (output == "flag") {
    df.aug <- select(out, .data$ID, !!!gr_by, !!! ls.vars[["original"]])
  }
  if (output == "complete") {
    df.aug <- bind_rows(
      select(df.aug ,.data$ID, !!!gr_by, !!! ls.vars[["ion1"]]),
      select(df.aug, .data$ID, !!!gr_by, !!! ls.vars[["ion2"]])
      )
  }

  #if(output == "flag") df.aug <- select(df.aug, .data$ID, !!!gr_by, !!!var_names)

  # save results; flag ad statistics
  if (output == "flag") results <- select(out, -c(!!! ls.vars[["original"]]))
  if (output == "complete") results <- out
  out <- lst(df = df.aug, results = results)

  return(out)

}

#' @export
diag_R_exec <- function(df, method, args = expr_R(NULL), ..., output){

  gr_by <- enquos(...)

  # original variable names
  ls.vars <- var_fun(df, gr_by, args)

  diag_vc <- c("Cameca", "CooksD", "Rm", "CV", "QQ","norm_E")
  diag_method <- purrr::map(
    diag_vc,
    call2,
    expr(.),
    expr(args),
    expr(ls.vars),
    !!!gr_by,
    output = expr(output)
    ) %>%
    set_names(nm = diag_vc)

# Descriptive an predictive statistics for ion ratios
   stat_R(df,
          Xt = !! args[["Xt"]],
          N =!! args[["N"]],
          species = !! args[["species"]],
          ion1 = as_name(args[["ion1"]]),
          ion2 = as_name(args[["ion2"]]),
          !!! gr_by,
          output = "complete"
          ) %>%
  eval_tidy(expr = diag_method[[method]])

# Datafile with flag values associated to diagnostics
  # if (output == "flag") return(left_join(df, tb.aug, by = "ID"))
  # if (output == "complete") return(tb.aug)

  }


#' Reduce diagnostics
#'
#' @export
reduce_diag <- function(ls, type = "df", args = expr_R(NULL)){

  type_reduction <- function(type){
    switch(type,
           results = call2(
             "mutate",
             expr(.),
             execution =
             expr(as.numeric(execution) -1)
             ),
           df = call2(
             "mutate",
             expr(.),
             execution =
             expr(as.numeric(execution))
             )
          )
  }

  # Reduce the results to a single dataframe
  ls %>%
    purrr::transpose() %>%
    purrr::pluck(type) %>%
    bind_rows(.id = "execution") %>%
    eval_tidy(expr = type_reduction(type))

}

#' Create stat_R call quosure
#'
#' \code{expr_R} and \code{expr_cor} functions generate an quosure that mimics a
#' stat_R call or cor_IC call for subsequent usage in dia_R
#'
#' The \code{diag_R} function performs an internal call to stat_R or cor_IC call
#' to perform diagnostics on the influence of individial measurements on the
#' blockwise or global (i.e., a complete analysis) stastics. This function
#' provides a convenient way to enter this call as an quosure into the argument
#' \code{args} of \code{diag_R}.
#'
#' @param Xt A character string constituting the ion count rate.
#' @param N A character string constituting the ion counts.
#' @param species A character string constituting the species analysed.
#' @param ion1 A character string constituting the heavy isotope ("13C").
#' @param ion2 A character string constituting the light isotope ("12C").
#' @param t A character string constituting the time increment of measurement.
#'
#' @return A list containing the input arguments as a
#' \code{\link[rlang:quosure]{quosure}}.
#' @export
#' @examples
#' expr_R(Xt = "Xt.pr", N = "N.pr", species = "species.nm", ion1 = "13C",
#'        ion2 = "12C")
expr_R <- function(Xt, N, species, ion1, ion2){

  as_quosures(lst(Xt = parse_expr(Xt),
                  N = parse_expr(N),
                  species = parse_expr(species),
                  ion1 = ion1,
                  ion2 = ion2
                  ),
              env = caller_env()
              )

}


#-------------------------------------------------------------------------------
# Not exportet helper functions
#-------------------------------------------------------------------------------


# output_lm <- function(df, type){
#
#   if (type == "sum") vars <- quos(data, lm.0, lm.1, sum.lm)
#   if (type == "complete") vars <- quos(lm.0, lm.1, sum.lm)
#
#   switch(type,
#          sum = eval(call2("select", df, expr(-c(!!! vars)))),
#          complete = eval(call2("select", df, expr(-c(!!! vars))))
#          )
#
# }
