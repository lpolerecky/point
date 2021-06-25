#' Evaluate effect size and significance of outliers on R
#'
#' \code{eval_diag} function for the evaluation of effect size and significance
#' of outliers on R detected with diagnostics, such as Cook's D or sigma
#' rejection (Cameca default method).
#'
#' @param .IC A tibble containing ion count data and diagnostics generated with
#' \code{diag_R()}, as a minimum the outlier \code{flag} variable is required.
#' @param .ion1 A character string constituting the rare isotope (e.g. "13C").
#' @param .ion2 A character string constituting the common isotope (e.g. "12C").
#' @param ... Variables for grouping.
#' @param .nest A variable hat identifies a series of analyses to calculate
#' the significance of inter-isotope variability.
#' @param .X A variable constituting the ion count rate (defaults to
#' variables generated with \code{read_IC()})
#' @param .N A variable constituting the ion counts (defaults to variables
#' generated with \code{read_IC()}.).
#' @param .species A variable constituting the species analysed (defaults to
#' variables generated with \code{read_IC()}).
#' @param .t A variable constituting the time of the analyses (defaults to
#' variables generated with \code{read_IC()}).
#' @param .flag A variable constituting the outlier flag (defaults to
#' variables generated with \code{diag_R()}).
#' @param .execution A variable constituting the iterative cycles of diagnostics
#' (defaults to variables generated with \code{diag_R()}).
#' @param .output A character string for output as summary statistics
#' ("inference") and statistics with the original data ("complete").
#' @param .tf Variable transformation as parts per thousand (\code{"ppt"}) or
#' log (\code{"log"}) before mixed linear model application.
#' @param .label A character string indicating whether variable names are latex
#' (\code{"latex"}) or webtex (\code{"webtex"}) compatible. Will be extended in
#' the future \code{default = NULL}.
#' @param .meta Logical whether to preserve the metadata as an attribute
#' (defaults to TRUE).
#'
#' @return A \code{tibble::\link[tibble:tibble]{tibble}()} with model output.
#' See \code{point::names_model} for more information on the model results.
#'
#' @export
#' @examples
#' # Simulated IC data
#' tb_dia <- diag_R(simu_IC, "13C", "12C", type.nm, spot.nm,
#'                  .output = "diagnostic")
#'
#' # Evaluate significance and effect of outliers based on Cook's D
#' eval_diag(tb_dia, "13C", "12C", type.nm, spot.nm, .nest = type.nm)
#'
eval_diag <- function(.IC, .ion1, .ion2, ..., .nest = NULL, .X = NULL,
                      .N = NULL, .species = NULL, .t = NULL, .flag = NULL,
                      .execution = NULL, .output = "inference",
                      .tf = "ppt", .label = NULL, .meta = FALSE){

  # Quoting the call (user-supplied expressions)
  # Additional arguments
  args <- inject_args(
    .IC,
    enquos(.X = .X, .N = .N, .species = .species, .t = .t,
           .execution = .execution, .flag = .flag),
    type = c("processed", "group", "diagnostics"),
    check = FALSE
    )
  # Grouping and nesting
  gr_by <- enquos(...) %>% append(args[[".execution"]])
  nest <- enquo(.nest)

  # Metadata
  if(.meta) meta <- unfold(.IC, merge = FALSE)

  # Update quosures (rare and common isotope)
  args <- list2(
    !!! args,
    !!! all_args(args, .ion1, .ion2, chr = FALSE),
    # Rare isotope
    X1  = quo_updt(args[[".X"]], post = .ion1),
    # Common isotope
    X2  = quo_updt(args[[".X"]], post = .ion2)
    )

  # New quosures (model output)
  model_args <- arg_builder(args, "model")
  args <- list2(
    !!! args,
    # Model arguments
    !!! model_args,
    # Predicted rare isotope count rate
    hat_X1 = quo_updt(args[["X1"]], pre = "hat"),
    # Re-centred X
    Xe = quo_updt(args[[".X"]], post = "e", update_post = TRUE)
    )

  # Predicted rare isotope count rate
  if (!(as_name(args[["hat_X1"]]) %in% colnames(.IC))) {
    .IC <- mutate(.IC, !! args[["hat_X1"]] := !! args[["M_R"]] * !! args[["X2"]])
    }
  # Check number of levels of bad flag is more than 10
  IC_n <-  count(.IC, !! args[[".execution"]], !!! gr_by, !! args[[".flag"]])

  if (nrow(filter(IC_n , !! args[[".flag"]] == "divergent" & n >= 10)) == 0) {
    stop("Number of flagged outliers in all samples is too small for a reliable diagnostic. Execution has stopped.",
         call. = FALSE)
    }

  if (nrow(filter(IC_n, !! args[[".flag"]] == "divergent" & n >= 10)) <
      nrow(filter(IC_n, !! args[[".flag"]] == "divergent"))) {
    warning("Number of flagged outliers in some samples is too small for a reliable diagnostic. Execution proceeded with remaining samples.",
            call. = FALSE)
    # Otherwise filter data-set
    .IC <- filter(IC_n, !! args[[".flag"]] == "divergent" & n < 10)  %>%
      select(!!!gr_by) %>%
      anti_join(.IC, ., by = c(sapply(gr_by, as_name)))
    }

  # Check for ionization efficiency trend
  if (any(between(pull(.IC , !! args[["chi2_N2"]]) , 0.9, 1.1))) {
    warning("Linear ionization trend absent in some or all analyses; F statistic might be unreliable.",
            call. = FALSE)
    }

  # Re-centre residuals along flag variable
  IC <- cstd_var(.IC, gr_by, args)

  # Create zero (constrained) model flag and updated model
  IC_lm <- tidyr::nest(IC, data = -c(!!! gr_by)) %>%
    mutate(
      !! args[["ratio"]] := purrr::map(data, ~distinct(.x, !! args[["ratio"]])),
      !! args[["M_R"]] := purrr::map(data, ~distinct(.x, !! args[["M_R"]])),
      lm_out = purrr::map(data, lm_fun, args)
      ) %>%
    tidyr::unnest_wider(.data$lm_out)

  #
  if (.output == "inference") {
    IC_lm <- tidyr::unnest(IC_lm, cols = c(!! args[["ratio"]], !! args[["M_R"]]))
  } else {
    IC_lm <- select(IC_lm, -c(!! args[["ratio"]] ,!! args[["M_R"]]))
  }

  if (is_symbol(get_expr(nest))) {
    # Groups for nested data
    nest_args <- c(as_name(nest), as_name(args[[".execution"]]))
    nest_gr <- gr_by[!sapply(gr_by, as_name) %in% nest_args]

    # Nest over nest groups
    IC_mlm <- tidyr::nest(IC, data = -c(!!! nest_gr)) %>%
      mutate(
        gls_out =
          purrr::map(.data$data, purrr::possibly(gls_fun, NA), args, .tf = .tf),
        inter_out =
          purrr::map2(
            .data$data,
            .data$gls_out,
            purrr::possibly(mlm_fun, NA),
            args,
            nest,
            .tf = .tf
          )
      ) %>%
      select(-c(.data$gls_out, .data$data)) %>%
      tidyr::unnest_wider(.data$inter_out)

    # Collect
    IC_mlm <- list(IC_lm, IC_mlm)

    # Prepare output
    IC <- purrr::reduce(IC_mlm, left_join, by = sapply(nest_gr, as_name)) %>%
      output_lm(args, model_args, .meta, .label, .output)
    # Return metadata
    if (.meta) return(fold(IC, type = ".mt",  meta = meta)) else return(IC)
  }
  # Prepare output
  IC <- output_lm(IC_lm, args, model_args, .meta, .label, .output)
  # Return metadata
  if (.meta) fold(IC, type = ".mt",  meta = meta) else IC
}

#-------------------------------------------------------------------------------
# Not exportet helper functions
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# lm models
#-------------------------------------------------------------------------------

lm_fun <- function(.IC, args) {

  # full R model
  lm_1 <- lm_form(.IC, args[["Xe"]], args[["X2"]], flag = args[[".flag"]],
                  type = "Rm")
  # zero R model
  lm_0 <- lm_form(.IC, args[["Xe"]], args[["X2"]], type = "Rm")

  # Join model hypothesis test
  IC_aov <- broom::tidy(anova(lm_0 , lm_1))

  lst(
    !! args[["F_R"]] := pull(IC_aov, .data$statistic)[2],
    !! args[["p_R"]] := pull(IC_aov, .data$p.value)[2]
    )
}

#-------------------------------------------------------------------------------
# gls models
#-------------------------------------------------------------------------------

gls_fun <- function(.IC, args, .tf) {

  # zero model
  gls_0 <- lm_form(.IC, args[["X1"]], args[["X2"]], trans = .tf,
                   type = "GLS")

  lst(
    !! args[["hat_M_M_R"]] := coef_pull(gls_0, .IC, !! args[["X2"]], .tf),
    gls_0 = gls_0
    )
}

#-------------------------------------------------------------------------------
# mlm model inter R variability
#-------------------------------------------------------------------------------

mlm_fun <- function(.IC, .gls, args, nest, .tf) {

  # zero model
  gls_0 <- purrr::pluck(.gls, "gls_0")
  # mlm inter model
  mlm_inter <- lm_form(
    .IC,
    args[["X1"]],
    args[["X2"]],
    trans = .tf,
    vorce = "inter",
    nest =  nest,
    type = "LME"
    )
  # log likelihood test
  IC_aov <- anova(mlm_inter, gls_0)

  lst(
    # GLS results
    !! args[["hat_M_M_R"]] := purrr::pluck(.gls, as_name(args[["hat_M_M_R"]])),
    # model relative standard deviation of group and associated standard error
    !! args[["hat_RS_M_R"]] := mlm_RS(mlm_inter, args[["X2"]]),
    # test statistic
    !! args[["dAIC_M_R"]] := diff(pull(IC_aov, .data$`AIC`)),
    # p value
    !! args[["p_M_R"]] :=  zuur_cor(pull(IC_aov, .data$`L.Ratio`)[2])
    )
}

#-------------------------------------------------------------------------------
# output function
#-------------------------------------------------------------------------------

output_lm <- function(IC, args, model_args, meta = NULL, label = NULL, output) {

  # Output transform
  trans_out <- function(IC, output) {
    switch(
      output,
      inference = call2( "select", IC, expr(-.data$data), .ns = "dplyr"),
      complete = call2("unnest", IC, cols = expr(.data$data), .ns = "tidyr")
      )
  }
  IC <- eval(trans_out(IC, output))

  #Latex labels
  if (! is.null(label)) {
    tb_model <- point::names_model
    if (!is_symbol(get_expr(args[[".nest"]]))) {
      tb_model <- filter(point::names_model, .data$type == "Ratio method")
      # Model args augment
      model_args <- model_args[paste(tb_model$name, tb_model$derived, sep = "_")]
    }
    ls_latex <- set_names(
      sapply(model_args, as_name),
      tex_labeller(tb_model, tb_model$name, label)
    )
    IC <- rename(IC, !!! ls_latex)
  }

  # Return metadata
  if (!is.null(meta)) IC <- fold(IC, type = ".mt", meta = meta)

  return(IC)
  }

# get coeffecients from mlm model
coef_pull <- function(sum, data, arg, trans){

  cf <- unname(coef(sum))
  if (trans == "ppt") return(cf / 1000)
  if (trans == "log") return(trans_R(data, arg = arg, cf = cf))

}

# re-center outliers towards the fitted value by subtraction of residual
# extremes
cstd_var <- function(.IC, gr_by, args){
  group_by(.IC, !!! gr_by) %>%
    mutate(
      # Residuals
      E = !! args[["X1"]] - !! args[["hat_X1"]],
      # Divide values in upper and lower sectors
      sector = if_else(!! args[["X1"]] >=  !! args[["hat_X1"]], "upper", "lower")
      ) %>%
    group_by(!!! gr_by, !! args[[".flag"]], .data$sector) %>%
    mutate(
      min_bound = if_else(.data$sector == "upper", min(.data$E), max(.data$E)),
      !! args[["Xe"]] := .data$E - .data$min_bound + !! args[["hat_X1"]],
      .keep = "unused"
      ) %>%
    ungroup() %>%
    select(-c(.data$min_bound, .data$sector))
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
    filter(stringr::str_detect(.data$name, (as_name(arg))))

  # Coefs
  ran <- as.numeric(tibble::deframe(ran[2,2]))
  fix <- nlme::fixed.effects(sum) %>% unname()
  # Relative variance of the slope
  RS <-  ran / fix

  if (output == "value") {return(RS * 1000)} # per mille

}

# correction for testing on the boundary Zuur et al 2008 (CH 5, p 123)
zuur_cor <- function(L) {0.5 * (1 - pchisq(L, 1))}
