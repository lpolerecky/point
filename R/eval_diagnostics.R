#' Evaluate effect size and significance of outliers on R
#'
#' \code{eval_diag} function for the evaluation of effect size and significance
#' of outliers on R detected with diagnostics, such as Cook's D or sigma
#' rejection (Cameca default method).
#'
#' @param .df A tibble containing ion count data and diagnostics generated with
#' \code{diag_R()}, as a minimum a flag variable is required.
#' @param .ion1 A character string constituting the heavy isotope ("13C").
#' @param .ion2 A character string constituting the light isotope ("12C").
#' @param ... Variables for grouping.
#' @param .nest A variable over which to created the random component of a mixed
#' linear model (default = NULL).
#' @param .Xt A variable constituting the ion count rate (defaults to
#' variables generated with \code{read_IC()})
#' @param .N A variable constituting the ion counts (defaults to variables
#' generated with \code{read_IC()}.).
#' @param .species A variable constituting the species analysed (defaults to
#' variables generated with \code{read_IC()}).
#' @param .t A variable constituting the time of the analyses (defaults to
#' variables generated with \code{read_IC()}).
#' @param .flag A variable constituting the outlier flag (defaults to
#' variables generated with \code{diag_R()}).
#' @param .output A character string for output as summary statistics ("sum")
#' and statistics with the original data ("complete").
#'
#' @return A \code{\link[tibble::tibble]{tibble}} with model output. See
#' \code{point::names_model} for more information on the results.
#'
#' @export
#' @examples
#' # Use point_example() to access the examples bundled with this package
#' # raw data containing 13C and 12C counts on carbonate
#' tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # Processing raw ion count data
#' tb_pr <- cor_IC(tb_rw)
#'
#' # Diagnostics
#' tb_dia <- diag_R(tb_pr, "13C", "12C", file.nm, sample.nm)
#'
#' # evaluation of diagnostics
#' eval_diag(tb_dia, "13C", "12C", file.nm, sample.nm)
eval_diag <- function(.df, .ion1, .ion2, ..., .nest = NULL,
                      .Xt = Xt.pr, .N = N.pr, .species = species.nm,
                      .t = t.nm, .flag = flag, .output = "sum", .tf = "ppt",
                      .label = "none"){

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)
  # Additional arguments
  args <- enquos(Xt = .Xt, N = .N, species = .species, t = .t, flag = .flag,
                 nest = .nest)

  # heavy isotope
  Xt1 <- quo_updt(args[["Xt"]], post = .ion1)
  # light isotope
  Xt2 <- quo_updt(args[["Xt"]], post = .ion2)
  N2 <- quo_updt(args[["N"]], post = .ion2)
  # Mean R analysis
  M_R <- quo_updt(args[["Xt"]], pre = "M_R")
  # Mean R group
  M_R.gr <- quo_updt(expr(M_R), post = as_name(args[["nest"]]))
  # Predicted heavy isotope
  hat_Xt1 <- quo_updt(Xt1, pre = "hat")

  # latex model variable names
  if (length(unique(pull(.df, .data$execution))) > 1) {
    labs <- point::names_model
    } else{
      labs <- point::names_model[1:13,]
    }
  ls_latex <- set_names(labs$name, nm = labs$latex)

  ls_latex[1] <- as_name(M_R)
  ls_latex[7] <- as_name(M_R.gr)

  if (!(as_name(hat_Xt1) %in% colnames(.df))) {
    .df <- mutate(.df, !!hat_Xt1 := !!M_R * !!Xt2)
    }
  # Check number of levels of bad flag is more than 10
  df_n <-  count(.df, .data$execution, !!!gr_by, !!args[["flag"]])

  if (nrow(filter(df_n, !!args[["flag"]] == "divergent" & n > 10)) == 0) {
  stop("Number of flagged outliers in all samples is too small for a reliable
       diagnostic. Execution has stopped.")
    }

  if (nrow(filter(df_n, !!args[["flag"]] == "divergent" & n > 10)) <
      nrow(filter(df_n, !!args[["flag"]] == "bad"))) {
    warning("Number of flagged outliers in some samples is too small for a
            reliable diagnostic. Execution proceeded with remaining samples.")
    # Otherwise filter data-set
    .df <- filter(df_n, !!args[["flag"]] == "bad" & n < 10)  %>%
      select(.data$execution, !!!gr_by) %>%
      anti_join(.df, ., by = c("execution", sapply(gr_by, as_name)))
    }

  # Check for ionization trend
  chi2 <- group_by(filter(.df, .data$execution == 1), !!!gr_by) %>%
    summarise(chi2 = ((sd(!!Xt2) / sqrt(n())) / (sqrt(mean(!!N2) / n()))) ^ 2)

  if (any(between(chi2$chi2, 0.9, 1.1))) {
    warning("Linear ionization trend absent in some or all analyses;
            F statistic might be unreliable.")
    }

  # Re-center along flag variable
  df <- cstd_var(.df, Xt1, hat_Xt1, args[["flag"]], !!! gr_by, execution)

  # Create zero (constrained) model flag and updated model
  df_lm <- tidyr::nest(df, data = -c(!!! gr_by,  .data$execution, !! M_R)) %>%
    mutate(
      lm_out =
        purrr::map(
          data,
          ~lm_fun(.x, .Xt1 = quo(std.var), .Xt2 = Xt2, .flag = args[["flag"]])
          )
      ) %>%
    unnest_wider(lm_out)

  if (is_symbol(get_expr(args[["nest"]]))) {
    # Groups for nested data
    nest_gr <- gr_by[!sapply(gr_by, as_name) %in% as_name(args[["nest"]])]

    df_mlm <- df %>%
    # Nest over nest groups
      tidyr::nest(data = -c(!!! nest_gr))

    df_mlm1 <- mutate(
      df_mlm,
      gls_out =
        purrr::map(
          data,
          purrr::possibly(gls_fun, NA),
          .Xt1 = Xt1,
          .Xt2 = Xt2,
          .tf = .tf,
          .group = args[["nest"]]
          ),
      inter_out =
        purrr::map2(
          data,
          gls_out,
          purrr::possibly(mlm_fun1, NA),
          .Xt1 = Xt1,
          .Xt2 = Xt2,
          .tf = .tf,
          .group = args[["nest"]]
          )
      ) %>%
        hoist(.col = gls_out, "ls_gls") %>%
        select(-c(gls_out, data)) %>%
        unnest_wider(ls_gls) %>%
        unnest_wider(inter_out)

    ls_mlm <- lst(df_lm, df_mlm1)

      # Check if longitudinal analyses can be performed
      if (length(unique(pull(df, .data$execution))) > 1) {

        df_mlm2 <- mutate(
          df_mlm,
          intra_out =
            purrr::map(
              data,
              purrr::possibly(mlm_fun2, NA),
              .Xt1 = Xt1,
              .Xt2 = Xt2,
              .tf = .tf,
              .group = args[["nest"]]
              )
          ) %>%
          select(-data) %>%
          unnest_wider(intra_out)

        ls_mlm$df_mlm2 <- df_mlm2

      } else {
        warning("Longitudinal analyses cannot be peformed.")
      }

    # Prepare output
    df <- purrr::reduce(ls_mlm, left_join, by = sapply(nest_gr, as_name)) %>%
      eval_tidy(expr = output_lm(.output))

    if (.label == "latex") return(rename(df, !!!ls_latex))
    return(df)
  }
  # Prepare output
  df <- df_lm %>%
    eval_tidy(expr = output_lm(.output))
  if (.label == "latex") return(rename(df, !!!ls_latex))
  return(df)
}

#-------------------------------------------------------------------------------
# Not exportet helper functions
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# lm models
#-------------------------------------------------------------------------------

lm_fun <- function(.df, .Xt1, .Xt2, .flag) {

  ls_lm <- lst()

  # full R model
  lm_1 <- lm_form(.df, .Xt1, .Xt2, flag = .flag, type = "Rm")
  # zero R model
  lm_0 <- lm_form(.df, .Xt1, .Xt2, type = "Rm")
  # Effect size
  df_f <- tibble(effectsize::cohens_f(lm_1, model2 =  lm_0))
  ls_lm$f <- pull(df_f , Cohens_f_partial)
  ls_lm$f_cl <- pull(df_f , CI_low)
  ls_lm$f_cu <- pull(df_f , CI_high)
  # Join model hypothesis test
  df_aov <- broom::tidy(anova(lm_0 , lm_1))
  ls_lm$F_vl <- pull(df_aov, statistic)[2]
  ls_lm$p_F <- pull(df_aov, p.value)[2]

  return(ls_lm)
}

#-------------------------------------------------------------------------------
# gls models
#-------------------------------------------------------------------------------

gls_fun <- function(.df, .Xt1, .Xt2, .tf, .group) {

  # group-wise ratio of arithmetic mean
  M_R.gr <- quo_updt(expr(M_R), post = as_name(.group))

  # zero model
  gls_0 <- lm_form(.df, .Xt1, .Xt2, trans = .tf, type = "GLS")
  # log-log zero model
  log_gls_0 <- lm_form(.df, .Xt1, .Xt2, trans = "log", type = "GLS")

  ls_gls <- lst(
    # group-wise ratio of arithmetic mean (R)
    !!M_R.gr := coef_pull(gls_0, .df, .Xt2, .tf),
    # group-wise ratio of geometric mean (logR)
    GM_R = coef_pull(log_gls_0, .df, .Xt2, "log"),
    # comparison of geometric and arithmetic R
    dNorm = (GM_R / !!M_R.gr - 1) * 1000
    )

  return(lst(ls_gls, gls_0 = gls_0))
}

#-------------------------------------------------------------------------------
# mlm model inter R variability
#-------------------------------------------------------------------------------

mlm_fun1 <- function(.df, .gls, .Xt1, .Xt2, .tf, .group) {

  # zero model
  gls_0 <- purrr::pluck(.gls, "gls_0")
  # mlm inter model
  mlm_inter <- lm_form(
    .df,
    .Xt1,
    .Xt2,
    trans = .tf,
    vorce = "inter",
    nest = .group,
    type = "LME"
    )
  # log likelihood test
  df_aov <- anova(mlm_inter, gls_0)

  ls_inter <- lst(
    # model relative standard deviation of group and associated standard error
    RS_R_inter = mlm_RS(mlm_inter, .Xt2),
    RS_R_se_inter = mlm_RS(mlm_inter, .Xt2, output = "se"),
    # test statistic
    dAIC_inter = diff(pull(df_aov, `AIC`)),
    # p value
    p_inter = pull(df_aov, `p-value`)[2]
    )

  return(ls_inter)

}

#-------------------------------------------------------------------------------
# mlm model intra R variability
#-------------------------------------------------------------------------------

mlm_fun2 <- function(.df, .Xt1, .Xt2, .tf, .group) {

  # mlm intra model
  mlm_intra <- lm_form(.df, .Xt1, .Xt2, vorce = "intra", nest = .group,
                      type = "LME"
                      )

  # anova to check whether interaction with execution number is significant
  df_aov <- broom::tidy(car::Anova(mlm_intra, type = 3))

  ls_intra <- lst(
    # model relative standard deviation of group and associated standard error
    RS_R_intra = mlm_dR(mlm_intra, .Xt2),
    RS_R_se_intra = mlm_dR(mlm_intra, .Xt2, output = "se"),
    # p value
    p_intra = pull(df_aov, p.value)[2]
  )

  return(ls_intra)

}

#-------------------------------------------------------------------------------
# output function
#-------------------------------------------------------------------------------

output_lm <- function(.output) {

  switch(
    .output,
    sum = call2( "select", expr(.), expr(-data)),
    complete = call2("unnest", expr(.), cols = expr(data))
    )

}

coef_pull <- function(sum, data, arg, trans){

  cf <- unname(coef(sum))
  if (trans == "ppt") return(cf / 1000)
  if (trans == "log") return(trans_R(data, arg = arg, cf = cf))

}

# standardizing and re-center independent variable for fit to LM
cstd_var <- function(df, Xt1, hat_Y, flag, ...){

  gr_by <- enquos(...)

  group_by(df, !!! gr_by, !! flag) %>%
    mutate(range = if_else(!!Xt1 >= !!hat_Y, "upper", "lower")) %>%
    group_by(!!! gr_by, !! flag, .data$range) %>%
    mutate(
      max.range = if_else(
        .data$range == "upper",
        max(!!Xt1 - !!hat_Y),
        min(!!Xt1 - !!hat_Y)
        ),
      min.range = if_else(
        .data$range == "upper",
        min(!!Xt1 - !!hat_Y) ,
        max(!!Xt1 - !!hat_Y)
        ),
      range.val = abs(.data$max.range - .data$min.range)
      ) %>%
    mutate(
      std.var =
        if_else(
          .data$range == "upper",
          abs((!! Xt1 - !! hat_Y) -.data$min.range) / .data$range.val,
          -abs((!! Xt1 - !! hat_Y) -.data$min.range) / .data$range.val
          ) * .data$range.val + !! hat_Y
      ) %>%
    ungroup() %>%
    select(-c(.data$max.range, .data$min.range, .data$range, .data$range.val))
}

# Temporal trend of the fixed coefficient
mlm_dR <- function(sum, arg, output = "value") {

  fix <- nlme::fixed.effects(sum) %>% unname()
  dR <- fix[2] / fix[1]
  if (output == "value") return(dR * 1000)
  if (output == "se"){
    fix.sd <- broom.mixed::tidy(sum) %>%
      select(std.error) %>%
      tidyr::drop_na() %>%
      mutate(sd = std.error  * sqrt(nobs(sum)),
             mean = dR)

    dR.se <- ((sqrt(((fix.sd$sd[1]/ fix.sd$mean[1]) ^ 2) +
                      ((fix.sd$sd[2]/ fix.sd$mean[2]) ^ 2)) * dR) /
                sqrt(nobs(sum)))
    return(dR.se  * 1000) # per mille
  }
}

# Conditional coefficient back transformation
trans_R <- function(data, arg, cf){

  M_log_pred <- mean(log(pull(data, !!arg)))
  GM_pred <- exp(M_log_pred)
  GM_resp <- exp(M_log_pred * cf)
  GM_resp / GM_pred

}

# Relative standard deviation of the coefficient
mlm_RS <- function(sum, arg, output = "value") {

  ran <- (nlme::VarCorr(sum))[,2] %>%
    tibble::enframe() %>%
    filter(str_detect(name, (as_name(arg))))

  ran <- as.numeric(tibble::deframe(ran[2,2]))
  fix <- nlme::fixed.effects(sum) %>% unname()

  RS <-  ran / fix

  if (output == "value") {return(RS * 1000)} # per mille

  if (output == "se"){
    # unequal distribution CI (95 %) converted to sd with delta-method
    if (!is.null(nrow(sum$apVar))) {
      var_matrix <- sum$apVar
      par <- attr(var_matrix, "Pars")
      ran_sd <- msm::deltamethod(~ exp(x1)^2, par, var_matrix) * sqrt(nobs(sum))
      } else {
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
    RS.se <- ((sqrt(((unique(ran_sd) / unique(ran)) ^ 2) +
                      ((unique(fix_sd) / unique(fix)) ^ 2)) * RS) /
                sqrt(nobs(sum)))
    return(RS.se * 1000) # per mille
  }
}
