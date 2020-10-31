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
rerun_diag_R <- function(out,
                         input,
                         method,
                         args = expr_R(NULL),
                         ...,
                         output
                         ){

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
#' # evaluation of diagnostics
#' tb.eval <- eval_diag(ls.dia,
#'                      args = expr_R_stat,
#'                      flag = flag,
#'                      iso_offset,
#'                      simulation,
#'                      trend,
#'                      nest = TRUE,
#'                      group = simulation
#'                      )
eval_diag <- function(ls_df, args, flag, ..., nest = FALSE, group = NULL, output = "sum"){

  gr_by <- enquos(..., .named = TRUE)
  flag <- enquo(flag)
  group <- enquo(group)

  # heavy isotope
  Xt1 <- quo_updt(args[["Xt"]], as_name(args[["ion1"]]))
  # light isotope
  Xt2 <- quo_updt(args[["Xt"]], as_name(args[["ion2"]]))

  # Chi squared light isotope
  Chi <- quo_updt(Xt2, x = "chi2")

  # # analysis quo
  # analysis <- quo(analysis)

# Combine original/augmented datasets and results of diagnostics
  # df <- inner_join(
  #   reduce_diag(ls_df, "df", args, !!!gr_by),
  #   reduce_diag(ls_df, "results",  args, !!!gr_by),
  #   by = c("execution", "ID", sapply(gr_by, as_name))
  #   )

  df <- reduce_diag(ls_df, "results", args)#, !!!gr_by)
#   %>%
# # Create analysis label
#     mutate(analysis = "analysis")

# Check number of levels of bad flag is more than 10
  df <-  df %>%
    group_by(execution, !!!gr_by, !!flag) %>%
    tally() %>%
    filter(!!flag == "bad" & n < 10) %>%
    select(.data$execution, !!!gr_by) %>%
    ungroup() %>%
    anti_join(df, ., by = c("execution", sapply(gr_by, as_name)))

# Check for ionization trend
  if (any(between(pull(df, !! Chi), 0.9, 1.1))) {
    warning("Linear ionization trend absent in some or all analyses; F value might be unreliable.")
  }

# recenter along flag variable
  df <- cstd_var(df, Xt1, quo(hat_Y), flag, !!! gr_by, .data$execution)

# create zero (constrained) model flag and updated model
  df.lm <- df %>%
    tidyr::nest(data = -c(!!! gr_by, .data$execution)) %>%
    mutate(
        lm.1 = purrr::map(data,
                          ~lm_form(.x,
                                   quo(std.var),
                                   Xt2,
                                   flag = flag,
                                   type = "Rm"
                                   )
                           ),
        lm.0 = purrr::map(data,
                          ~lm_form(.x,
                                   quo(std.var),
                                   Xt2,
                                   type = "Rm"
                                   )
                           ),
        # hat_D_R = purrr::map2_dbl(lm.1, lm.0, ~{(unname(coef(.x))[1] / unname(coef(.y))[1] -1) * 1000}  ),
        # new_R = purrr::map_dbl(lm.1, ~{unname(coef(.x))[1] }  ),
        # old_R = purrr::map_dbl(lm.0, ~{unname(coef(.x))[1] }  ),
        sum.lm = purrr::map2(lm.0, lm.1, ~{broom::tidy(anova(.x, .y))}),
        F.val = purrr::map_dbl(sum.lm, ~{(pull(.x, statistic))[2]}),
        p.F = purrr::map_dbl(sum.lm, ~{(pull(.x, p.value))[2]})
        )


# mixed model
  if (nest) {
    df.mlm <- df %>%
      tidyr::nest(data = -c(!!!gr_by[!sapply(gr_by, as_name) %in% as_name(group)])) %>%
      mutate(
        # AM.0 = purrr::map(data, purrr::possibly(lm_form, NA), arg1 = Xt1, arg2 = Xt2, type = "Rm"),
        gls.0 = purrr::map(data, purrr::possibly(lm_form, NA), arg1 = Xt1, arg2 = Xt2, trans = "ppt", type = "GLS"),
        "AM_R.{{group}}" := purrr::map_dbl(gls.0, ~{unname(coef(.x)) / 1000}),
        log_gls.0 = purrr::map(data, purrr::possibly(lm_form, NA), arg1 = Xt1, arg2 = Xt2, trans = "log", type = "GLS"),
        "GM_R.{{group}}" := purrr::map2_dbl(data, log_gls.0, ~{trans_R(.x, Xt2, unname(coef(.y)))}),
        # "GM_R.{{group}}" := purrr::map2_dbl(data, GM.0, ~{trans_R(.x, arg = Xt2, unname(coef(.y)))}),
        R_diff = purrr::map2_dbl(!!quo_updt(quo(GM_R), as_name(group)), !!quo_updt(quo(AM_R), as_name(group)), ~{(.x/.y - 1) * 1000}),
        mlm.inter = purrr::map(data, purrr::possibly(lm_form, NA), arg1 = Xt1, arg2 = Xt2, trans = "ppt", vorce = "inter", nest = group, type = "LME"),
        log_mlm.inter = purrr::map(data, purrr::possibly(lm_form, NA), arg1 = Xt1, arg2 = Xt2, trans = "log", vorce = "inter", nest = group, type = "LME"),
        mlm.intra = purrr::map(data, purrr::possibly(lm_form, NA), arg1 = Xt1, arg2 = Xt2, vorce = "intra", nest = group, type = "LME"),
        RS_R.inter = purrr::map_dbl(mlm.inter, purrr::possibly(mlm_RS, NA_real_), arg = Xt2),
        RS_R_se.inter = purrr::map_dbl(mlm.inter, purrr::possibly(mlm_RS, NA_real_), arg = Xt2, output = "se"),
        log_RS_R.inter = purrr::map_dbl(log_mlm.inter, purrr::possibly(mlm_RS, NA_real_), arg = Xt2),
        log_RS_R_se.inter = purrr::map_dbl(log_mlm.inter, purrr::possibly(mlm_RS, NA_real_), arg = Xt2, output = "se"),
        dR_dt.intra = purrr::map_dbl(mlm.intra, purrr::possibly(mlm_dR, NA_real_), arg = Xt2),
        dR_dt_se.intra = purrr::map_dbl(mlm.intra, purrr::possibly(mlm_dR, NA_real_), arg = Xt2, output = "se"),
        # # # # sum.mlm = if_else(is.na(mlm.0) | is.na(mlm.1), NA_real_, purrr::map2(mlm.0, mlm.1, ~)),
        # # # # L.ratio  = if_else(is.na(mlm.0) | is.na(mlm.1), NA_real_, purrr::map2_dbl(mlm.0, mlm.1, ~{(pull(anova(.x, .y), L.Ratio))[2]})),
        p.inter = purrr::map2_dbl(mlm.inter, gls.0, purrr::possibly(function(x, y){ (pull(anova(x, y), `p-value`))[2] }, NA_real_)),
        log_p.inter = purrr::map2_dbl(log_mlm.inter, log_gls.0, purrr::possibly(function(x, y){ (pull(anova(x, y), `p-value`))[2] }, NA_real_)),
        p.intra = purrr::map_dbl(mlm.intra, purrr::possibly(function(x){ (pull(broom::tidy(car::Anova(x, type = 3)), p.value))[2] }, NA_real_))
            ) %>%
      select(-c(data, mlm.inter, log_mlm.inter, mlm.intra, gls.0, log_gls.0))

# output lm and mlm combined
    if (length(gr_by) == 1) {
      df <- output_lm(df.lm, output)
      df <- bind_cols(df, select(expand_grid(df.mlm, rep = 1:nrow(df)), -rep))
      if (output == "sum")  return(df)
      if (output == "complete")  return(tidyr::unnest(df, data))
    } else {
      df <- left_join(df.lm, df.mlm, by = sapply(gr_by[!sapply(gr_by, rlang::as_name) %in% rlang::as_name(group)], as_name)) %>%
        output_lm(.,output)
      if (output == "sum")  return(df)
      if (output == "complete")  return(tidyr::unnest(df, data))
    }

  }

# output lm only
  df <- output_lm(df.lm, output)
  if (output == "sum")  return(df)
  if (output == "complete")  return(tidyr::unnest(df, data))

  }

#' Reduce diagnostics
#'
#' @export
reduce_diag <- function(ls, type = "df", args = expr_R(NULL)){

  # gr_by <- enquos(...)

  type_reduction <- function(type){
    switch(type,
           results = call2("mutate",
                           expr(.),
                           execution =
                             expr(as.numeric(execution) -1)
                           ),
           df = call2("mutate",
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



# standardizing and re-center independent variable for fit to LM
#' @export
cstd_var <- function(df, Xt1, hat_Y, flag, ...){

  gr_by <- enquos(...)

  df %>%
    group_by(!!! gr_by, !! flag) %>%
    mutate(range = if_else(!!Xt1 >= !!hat_Y, "upper", "lower")) %>%
    group_by(!!! gr_by, !! flag, .data$range) %>%
    mutate(max.range = if_else(.data$range == "upper", max(!!Xt1 - !!hat_Y) ,
                               min(!!Xt1 - !!hat_Y)
                               ),
           min.range = if_else(.data$range == "upper", min(!!Xt1 - !!hat_Y) ,
                               max(!!Xt1 - !!hat_Y)
                               ),
           range.val = abs(.data$max.range - .data$min.range)
           ) %>%
    mutate(std.var = if_else(.data$range == "upper",
                             abs((!! Xt1 - !! hat_Y) -
                                   .data$min.range) / .data$range.val,
                             -abs((!! Xt1 - !! hat_Y) -
                                    .data$min.range) / .data$range.val
                             )
           * .data$range.val + !! hat_Y
           ) %>%
    ungroup() %>%
    select(-c(.data$max.range, .data$min.range, .data$range, .data$range.val))
}


# function to retrieve variable names of original or stats datasets
var_fun <- function(df, grps, args){

  vars <- colnames(df)
  # variables to be discarded
  pat_disc <- purrr::reduce(append(sapply(grps, as_name), "ID"),
                            str_c,
                            sep = "|"
                            )
  #original col names
  var_names <- vars %>%
    set_names() %>%
    purrr::discard(~str_detect(.x,  pat_disc))

  # variables to be saved

    wide_vars <- purrr::map(c(as_name(args[["ion1"]]),
                              as_name(args[["ion2"]])
                              ),
                              ~paste(vars, .x, sep = ".")
                            ) %>%
    purrr::reduce(append) %>%
    purrr::discard(~str_detect(.x, pat_disc))

    wide_vars.ion1 <- wide_vars[str_detect(wide_vars,
                                           as_name(args[["ion1"]])
                                           )
                                ] %>%
      set_names(nm = var_names)
    wide_vars.ion2 <- wide_vars[str_detect(wide_vars,
                                           as_name(args[["ion2"]])
                                           )
                                ] %>%
      set_names(nm = var_names)

 lst(original = var_names, ion1 = wide_vars.ion1, ion2 = wide_vars.ion2)
}

#' temporal trend of the fixed coefficient
#' @export
mlm_dR <- function(sum, arg, output = "value") {

  fix <- nlme::fixed.effects(sum) %>% unname()
  dR <- fix[2] / fix[1]
  if (output == "value") {return(dR * 1000)}
  if (output == "se"){
  fix.sd <- broom.mixed::tidy(sum) %>%
    select(std.error) %>%
    tidyr::drop_na() %>%
    mutate(sd = std.error  * sqrt(nobs(sum)),
           mean = dR)

  dR.se <- ((sqrt(((fix.sd$sd[1]/ fix.sd$mean[1]) ^ 2) + ((fix.sd$sd[2]/ fix.sd$mean[2]) ^ 2)) * dR) /  sqrt(nobs(sum)))
  return(dR.se  * 1000) # per mille
  }

}

#' Conditional coefficient back transformation
#' @export
trans_R <- function(data, arg, coef){

  M_log_pred <- mean(log(pull(data, !!arg)))
  GM_pred <- exp(M_log_pred)
  GM_resp <- exp(M_log_pred * coef)
  GM_resp / GM_pred

}

#' Relative standard deviation of the coefficient
#' @export
mlm_RS <- function(sum, arg, output = "value") {



  #if (trans) { arg <- paste0("I(",as_name(arg),"/1000)")} else { arg <- as_name(arg) }
  ran <- (nlme::VarCorr(sum))[,2] %>%
    tibble::enframe() %>%
    filter(str_detect(name, (as_name(arg))))

  ran <- as.numeric(tibble::deframe(ran[2,2]))
  fix <- nlme::fixed.effects(sum) %>% unname()

  RS <-  ran / fix

  if (output == "value") {return(RS * 1000)} # per mille

  if (output == "se"){
# unequal distribution CI (95 %) converted to sd with deltamethod
    if (!is.null(nrow(sum$apVar))) {
    var_matrix <- sum$apVar
    par <- attr(var_matrix, "Pars")
    ran_sd <- msm::deltamethod(~ exp(x1)^2, par, var_matrix) * sqrt(nobs(sum))
    } else{
      ran_sd <- 0
    }

# fixed effect CI (95 %) converted to sd
    fix_CI <- nlme::intervals(sum, which = "fixed")
    fix_sd <- fix_CI$fixed %>%
      tibble::as_tibble() %>%
      mutate(se = (upper - lower) / 3.92) %>%
      pull(se) * sqrt(nobs(sum))

    attr(fix_sd, "label") <- NULL

# propagation when calculating isotope RS
    RS.se <- ((sqrt(((unique(ran_sd) / unique(ran)) ^ 2) + ((unique(fix_sd) / unique(fix)) ^ 2)) * RS) /  sqrt(nobs(sum))) #* 1000
    return(RS.se * 1000) # per mille
  }

}


output_lm <- function(df, type){

  if (type == "sum") vars <- quos(data, lm.0, lm.1, sum.lm)
  if (type == "complete") vars <- quos(lm.0, lm.1, sum.lm)

  switch(type,
         sum = eval(call2("select", df, expr(-c(!!! vars)))),
         complete = eval(call2("select", df, expr(-c(!!! vars))))
         )

}
