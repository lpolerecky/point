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
#' diag_R(sim_IC_extremes,
#'        method = "CooksD",
#'        args = expr_R_stat,
#'        reps = 3,
#'        simulation,
#'        trend,
#'        isoscale = "VPDB"
#'        )
diag_R <- function(df,
                   method = "Cameca",
                   args = expr_R(NULL),
                   reps = 1,
                   ...,
                   plot = TRUE,
                   plot_type = "interactive",
                   iso = TRUE,
                   isoscale = NULL
                   ){

  gr_by <- enquos(...)

# repetitions
  max <- reps + 1

# empty list for iteration storage
  ls.tb <- rlang::rep_named(as.character(1:max), lst())

# Remove zeros
  if (method != "Cameca"){

    df <- zeroCt(df,
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
  df <- filter(df,
               !! args[["species"]] == !! args[["ion1"]] |
               !! args[["species"]] == !! args[["ion2"]]
               )
  ls.tb[[1]] <- lst(df = df, results = NULL)

# Execute repeated cycles of augmentation
  ls.tb <- ls.tb  %>%
    purrr::accumulate(rerun_diag_R,
                      method = method,
                      args = args,
                      !!! gr_by
                      )
  if (plot) {

    ls.tb %T>%
      {print(plot_diag_R(df = .,
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
rerun_diag_R <- function(out,
                         input,
                         method = "Cameca",
                         args = expr_R(NULL),
                         ...
                         ){

  gr_by <- enquos(...)

  # variables to be saved
  results_vars <- set_names(colnames(out$df))
  deselect_vars <- results_vars  %in% append(sapply(gr_by, as_name), "ID")

  variables <- parse_exprs(results_vars[!deselect_vars])

  # if called this way then output is set fixed (more fleible use with diag_R_exec)
  out <- diag_R_exec(out$df,
                     method = method,
                     args = args,
                     !!! gr_by,
                     output = "flag"
                     )

  # save augmented dataframe for next cycle
  df.aug <- filter(out, flag == "good") %>%
    select(.data$ID, !!!gr_by, !!!variables)

  # save results; flag ad statistics
  results <- select(out, -c(!!!variables))

  out <- lst(df = df.aug, results = results)

  return(out)

}

#' @export
diag_R_exec <- function(df,
                        method = "Cameca",
                        args = expr_R(NULL),
                        ...,
                        output
                        ){

  gr_by <- enquos(...)

  diag_vc <- c("Cameca", "CooksD", "Rm", "CV", "QQ","norm_E")
  diag_method <- purrr::map(diag_vc,
                            call2,
                            expr(.),
                            expr(args),
                            !!!gr_by,
                            output = expr(output)
                            ) %>%
    set_names(nm = diag_vc)

# Descriptive an predictive statistics for ion ratios
    tb.aug <- stat_R(df,
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
  if (output == "flag") return(left_join(df, tb.aug, by = "ID"))
  if (output == "complete") return(tb.aug)

  }

#' Evaluate effect size of diagnostics
#'
#' \code{eval_diag} function for the evaluation of effect size of diagnostics
#'
#' @param ls_df A list of tibbles as generated by \code{diag_R}.
#' @param args A list of quosures pertaining to the variables required for a
#' call to stat_R. The function expr_R can be used to simplify setting this
#' argument.
#' @param df2 A dataframe with the ionization tend of the major ion in per
#' mille.
#' @param RSXt A variable constituting the relative standard deviation of the
#' major ion in the dataframe enetered for the argument \code{df2}.
#' @param n A variable indicating the number of observation in the list of
#' dataframes entered for the argument \code{ls_df}
#' @param ... Variables for grouping.
#' @param output A character string for output as summary statistics ("sum") and
#' statistics with the original data ("complete").
#'
#' @export
#' @examples
#'
#' # estimate strength of ionization trend as variation in major ion in per mille
#' tb.var12C <- predict_var(tb.pr, Xt.pr, N.pr, t.rw, species.nm, "12C", file.nm)
#'
#' # evaluation of diagnostics
#' tb.eval <- eval_diag(ls.diag,
#'                      args = expr_R_stat,
#'                      tb.var12C,
#'                      RS_Xt.pr.12C,
#'                      n_R_Xt.pr,
#'                      file.nm,
#'                      bl.mt,
#'                      output = "sum"
#'                      )
eval_diag <- function(ls_df,
                      args,
                      df2,
                      RSXt,
                      n,
                      ...,
                      output = "sum"
                      ){

  RSR <- quo_updt(args[["Xt"]], x = "RS_R")
  RSXt <- enquo(RSXt)
  n <- enquo(n)
  gr_by <- enquos(...)

# connecting the
  gr_vc <- sapply(gr_by, as_name)
# connecting the single ion variable
  uni_gr <- gr_vc[gr_vc %in% colnames(df2)]

  df  <- reduce_diag(ls_df,
                     type = "df",
                     args = args,
                     !!! gr_by
                     )

  ls.trans <- lst(
    quo(!! n),
    quo(c(diff(!! n, 1), NA)),
    quo(c(diff(!! RSR, 1), NA)),
    quo(diff_RS ^ 2 * (diff_n - 1)),
    quo((!! RSR) ^ 2 * (!! n - 1)),
    quo(SSD / SSA),
    quo(n2 * nz)
    )

# Names of parameters for diagnostic evaluation
  ls.names <- c("nz", "diff_n", "diff_RS", "SSD", "SSA", "n2", "chi_n2")

# Set names
  ls.trans <- set_names(ls.trans, nm = ls.names)


# Switch output complete dataset, stats or summary stats
  mod_cal <- function(type) {
    switch(type,
           complete = call2( "mutate", expr(.), !!! ls.trans),
           sum = call2( "transmute",  expr(.), !!! ls.trans)
    )
    }

  ls.preserve <-lst(
    quo(sum(n2, na.rm = TRUE)),
# check this it can be sum, only works for eta (fraction)!!!!
    quo(sum(chi_n2, na.rm = TRUE))
    )

# Names to be preserved
  ls.names2 <- c("n2", "chi_n2")

# Set names
  ls.preserve <- set_names(ls.preserve, nm = ls.names2)

# Switch to determine what extend of the data is preserved
  filter_cal <- function(type) {
    switch(type,
           complete = call2("select",
                            expr(.),
                            expr(-c(ls.names[!ls.names %in% ls.names2]))
                            ),
           sum = call2("summarise",
                       expr(.),
                       !!! ls.preserve
                       )
           )
    }


# Evaluate expressions and calls
    df_eval <- df %>%
      group_by(!!! gr_by) %>%
      arrange(desc(execution)) %>%
      eval_tidy(expr = mod_cal(output)) %>%
      eval_tidy(expr = filter_cal(output))

    if (output == "sum") {
      df <- lst(filter(df, execution == 1), df_eval, df2) %>%
        purrr::reduce2(lst(gr_vc, uni_gr), left_join) %>%
        mutate(crit_val = purrr::map_dbl(!! RSXt, crit_size),
               flag = if_else(.data$chi_n2 > crit_val, "variable", "non-variable")
               )
      return(df)
    }

    if (output == "complete") {
      df <- left_join(df_eval, df2, by = uni_gr) %>%
        mutate(crit_val = purrr::map_dbl(!! RSXt, crit_size),
               flag = if_else(.data$chi_n2 > crit_val, "variable", "non-variable")
               )
      return(df)
    }


  }

#' Reduce diagnostis
#'
#' @export
reduce_diag <- function(ls, type = "df", args = expr_R(NULL), ..., output = "sum"){

  gr_by <- quos(...)

  type_reduction <- function(type){
    switch(type,
           results = call2("bind_rows", expr(.), .id = "execution"),
           df = call2("map_dfr",
                      expr(.),
                      expr(~stat_R(df = .x,
                                   Xt = !! args[["Xt"]],
                                   N = !! args[["N"]],
                                   species = !! args[["species"]],
                                   ion1 = as_name(args[["ion1"]]),
                                   ion2 = as_name(args[["ion2"]]),
                                   !!! gr_by,
                                   output = output
                      )
                      ),
                      .id = "execution",
                      .ns = "purrr"
           )
    )
  }

  # Reduce the results to a single dataframe
  ls %>%
    purrr::transpose() %>%
    purrr::pluck(type) %>%
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

crit_size <- function(RS_Xt.ion2){

  null_dist <- null_dist %>%
    mutate(min.val = abs(RS_Xt.ion2 - parse_number(trend))) %>%
    filter(min.val == min(.data$min.val)) %>%
    pull(y_star)

}
