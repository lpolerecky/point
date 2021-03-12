#' Analyse raw ion count data
#'
#' \code{stat_Xt} function for descriptive and predictive statistics on single
#' ion precision.
#' \code{stat_R} function for descriptive and predictive statistics on isotope
#' ratios (R) precision with appropriate error propagation.
#'
#' These functions are a convenient wrapper to calculate the statistics
#' pertaining to the precision of pulsed ion count data (e.g. secondary ion mass
#' spectrometry). The statistics can either be calculated for single ions or
#' isotope ratios, and include observed and predicted (Poisson) statistics.
#' Calculations for isotope ratios include proper error propagation. For more
#' information on the usage as well as the mathematics behind these
#' functions see \code{vignette("IC-precision", package = "point")}.
#'
#' @param .IC A tibble containing processed ion count data.
#' @param .ion1 A character string constituting the heavy isotope ("13C").
#' @param .ion2 A character string constituting the light isotope ("12C").
#' @param ... Variables for grouping.
#' @param .Xt A variable constituting the ion count rate (defaults to
#' variables generated with \code{read_IC()}.)
#' @param .N A variable constituting the ion counts (defaults to variables
#' generated with \code{read_IC()}.).
#' @param .species A variable constituting the species analysed (defaults to
#' variables generated with \code{read_IC()}.).
#' @param .t A variable constituting the time of the analyses (defaults to
#' variables generated with \code{read_IC()}.).
#' @param .stat Select statistics (e.g. \code{c("M", "RS"}), see the tables
#' \code{point::names_stat_Xt} and \code{point::names_stat_R} for the full
#' selection of statistics available (default uses all statistic
#' transformations).
#' @param .label A character string indicating whether variable names are latex
#' compatible. Will be extended in the future \code{default = NULL}.
#' @param .output A character string for output as summary statistics ("sum");
#' statistics only ("stat"); and statistics with the original data ("complete")
#' \code{default = "sum"}.
#' @param .zero A character string that determines whether analyses with zero
#' count measurements will be removed from the calculations.
#'
#' @return A \code{tibble::\link[tibble:tibble]{tibble}} containing descriptive
#' and predictive statistics for ion counts and isotope ratios. The naming
#' convention depends on the argument \code{latex}; if set to \code{FALSE},
#' variable names concerning statistics will consist of an abbreviation pasted
#' together with the input variable names of \code{Xt}.
#' @export
#' @examples
#' # Use point_example() to access the examples bundled with this package
#'
#' # raw data containing 13C and 12C counts on carbonate
#' tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # Processing raw ion count data
#' tb_pr <- cor_IC(tb_rw)
#'
#' # Single ion descriptive an predictive statistics for all measured ions
#' stat_Xt(tb_pr, file.nm)
#'
#' # Descriptive an predictive statistics for 13C/12C ratios
#' stat_R(tb_pr, "13C", "12C", file.nm, .zero = TRUE)
#'
stat_Xt <- function(.IC, ..., .Xt = Xt.pr, .N = N.pr, .species = species.nm,
                    .t = t.nm, .stat = point::names_stat_Xt$name, .label = NULL,
                    .output = "sum"){

  stopifnot(tibble::is_tibble(.IC))
  if (.output != "sum"  & !is.null(.label)) {
    stop("Latex labels is not supported for complete datasets.")
  }

  # Quoting the call (user-supplied expressions)
  args <- enquos(.Xt = .Xt, .N = .N, .species = .species, .t = .t)

  # Grouping
  gr_by <- enquos(...)

  # New quosures
  new_args <- arg_builder(args, "Xt")
  args <- append(args, new_args)

  # The statistics
  calcs <- quos(
    # number of measurements
    n(),
    # sum of ion counts (important for external reproducibility calcs)
    sum(!! args[[".N"]]),
    # mean ion count rate
    mean(!! args[[".Xt"]]),
    # standard deviation (SD) count rate
    sd(!! args[[".Xt"]]),
    # RSD ion count rate
    (!! args[["S_Xt"]] / !! args[["M_Xt"]]) * 100,
    # standard error of the mean (SE) count rate
    sd(!! args[[".Xt"]]) / sqrt(n()),
    # predicted SD count rate
    sqrt(mean(!! args[[".N"]])),
    # predicted RSD count rate
    (1 / sqrt(mean(!! args[[".N"]]))) * 100,
    # predicted SE count rate
    sqrt(mean(!! args[[".N"]]) / n()),
    # reduced chi squared
    (!! args[["SeM_Xt"]] / !! args[["hat_SeM_N"]]) ^ 2
    )

  # The statistic names (depend on user-supplied expression)
  ls_nm <- sapply(new_args, as_name)

  # Set statistic names
  calcs <- set_names(calcs , nm = ls_nm)

  # Stat selection
  str_stat <- stringr::str_c(paste0("^", .stat), collapse = "|")
  pos_vars <- purrr::keep(ls_nm, ~stringr::str_detect(., str_stat))
  neg_vars <- purrr::discard(ls_nm, ~stringr::str_detect(., str_stat))

  # Render latex variable names
  if (!is.null(.label)) {
    ls_latex <- set_names(
    pos_vars,
    nm = purrr::map_chr(.stat, ~stat_labeller(stat = .x, label = .label))
    )
  }

  # Evaluate expressions
  IC <- group_by(.IC, !!! gr_by, !! args[[".species"]]) %>%
    eval_tidy(expr = mod_cal(.output, calcs = calcs)) %>%
    ungroup() %>%
    select(-any_of(neg_vars))

  # Output
  if (!is.null(.label)) {
    if (.label == "latex" | .label == "webtex") {
      IC <- mutate(
        IC,
        !! args[[".species"]] :=
          purrr::map_chr(
            !! args[[".species"]],
            ion_labeller,
            label = .label
            )
        ) %>%
        select(!!! gr_by, !! args[[".species"]], !!! ls_latex)
    return(IC)
    }
    }
  return(IC)
  }

#' @rdname stat_Xt
#'
#' @export
stat_R <- function(.IC, .ion1, .ion2, ..., .Xt = Xt.pr, .N = N.pr,
                   .species = species.nm, .t = t.nm, .nest = NULL,
                   .stat = point::names_stat_R$name, .label = NULL,
                   .output = "sum", .zero = FALSE){

  stopifnot(tibble::is_tibble(.IC))

  if (.output != "sum"  & !is.null(.label)) {
    stop("Latex labels is not supported for complete datasets.")
  }

  # Quoting the call (user-supplied expressions)
  args <- enquos(.Xt = .Xt, .N = .N, .species = .species, .t = .t)
  # Grouping
  gr_by <- enquos(...)
  # Nesting
  nest <- enquos(.nest)
  # Set labels for internal reproducibility
  type <- "internal"

  # External precision
  if (!any(sapply(nest, function(x) is.null(get_expr(x))))) {

    # Calculate single ion stat to obtain mean and total ion counts
    data_env <- env(data = .IC)
    IC <- call2("stat_Xt", .IC, !!! gr_by, !!! args, .ns = "point") %>%
      eval(envir = data_env)

    # Updated quotes
    args[[".Xt"]] <- quo_updt(args[[".Xt"]], pre = "M")
    args[[".N"]] <- quo_updt(args[[".N"]], pre = "tot")
    # Updated grouping
    gr_by <- nest
    # Set labels for external reproducibility
    if (!is.null(.label)) type <- "external"

    # Check whether groups have enough observations
    IC <- add_count(IC, !!! nest, !! args[[".species"]])
    if (any(IC$n <= 1)) {
      warning("Some groups have too few observations for a reliable estimation
              of the external precision. These groups have been omitted.")
      .IC <- filter(IC, n > 1)
    }
  }

  # If t column is empty create a manual time increment
  if (!(as_name(args[[".t"]]) %in% colnames(.IC))) {
    .IC <- group_by(.IC, !!!gr_by, !! args[[".species"]]) %>%
      mutate(!! args[[".t"]] := row_number()) %>%
      ungroup()
  }

  # Remove white space in ion names and add underscore for polyatomic species
  .ion1 <- ion_trim(.ion1)
  .ion2 <- ion_trim(.ion2)

  # Update quosures (rare and common isotope)
  args <- list2(
    !!! args,
    # Update quosures (rare isotope)
    Xt1 = quo_updt(args[[".Xt"]], post = .ion1), # count rate
    N1 = quo_updt(args[[".N"]], post = .ion1), # counts
    # Update quosures (common isotope)
    Xt2 = quo_updt(args[[".Xt"]], post = .ion2), # count rate
    N2 = quo_updt(args[[".N"]], post = .ion2) # counts
    )

  # New quosures
  new_args <- arg_builder(args, "R")
  args <- append(args, new_args)

  # The statistics
  calcs <- quos(
    # number of measurements
    n(),
    # mean isotope ratio
    mean(!! args[["Xt1"]]) / mean(!! args[["Xt2"]]),
    # SD isotope ratio
    stat_SDprop(!! args[["Xt1"]], !! args[["Xt2"]]),
    # RSD isotope ratio
    (!! args[["S_R_Xt"]] / !! args[["M_R_Xt"]]) * 1000,
    # SE isotope ratio
    !! args[["S_R_Xt"]] / sqrt(n()),
    # RSE isotope ratio
    (!! args[["SeM_R_Xt"]] / !! args[["M_R_Xt"]]) * 1000,
    # predictive SD isotope ratio
    stat_SDprop(!! args[["N1"]], !! args[["N2"]], predicted = TRUE),
    # predictive RSD isotope ratio
    !! args[["hat_S_R_N"]] / !! args[["M_R_Xt"]] * 1000,
    # predictive SE isotope ratio
    !! args[["hat_S_R_N"]]  / sqrt(n()),
    # predictive RSE isotope ratio
    !! args[["hat_SeM_R_N"]]  / !! args[["M_R_Xt"]] * 1000,
    # reduced chi squared
    (!! args[["SeM_R_Xt"]] / !! args[["hat_SeM_R_N"]]) ^ 2
    )

  # The statistic names (depend on user-supplied expression)
  ls_nm <- sapply(new_args, as_name)

  # Set statistic names
  calcs <- set_names(calcs, nm = ls_nm)

  # Extra arg in case of mutate and transmute (isotope ratio for each time step)
  if (.output != "sum") {
    R <- paste("R", as_name(args[[".Xt"]]), sep = "_")
    calcs[[R]] <- quo(!! args[["Xt1"]] / !! args[["Xt2"]])
  }

  # Stat selection
  str_stat <- stringr::str_c(paste0("^", paste0(.stat, "_R")), collapse = "|")
  pos_vars <- purrr::keep(ls_nm, ~stringr::str_detect(., str_stat))
  neg_vars <- purrr::discard(ls_nm, ~stringr::str_detect(.,  str_stat))

  # Render nice latex variable names
  if (!is.null(.label)) {
    ls_latex <- set_names(
      pos_vars,
      nm = purrr::map_chr(
        .stat,
        ~stat_labeller(var = "R", stat = .x, type = type, label = .label)
        )
      )
  }

  # Evaluate expressions and calls
  IC <- .IC  %>%
    eval_tidy(
      expr =
        zeroCt_cal(.zero, .ion1, .ion2, gr_by, args[[".N"]], args[[".species"]],
                   args[[".t"]])
              ) %>%
    cov_R(c(.ion1, .ion2), !!! gr_by, .species = !! args[[".species"]],
          .t = !! args[[".t"]]) %>%
    group_by(!!! gr_by) %>%
    eval_tidy(expr = mod_cal(.output, calcs = calcs)) %>%
    ungroup() %>%
    select(-any_of(neg_vars))

# Output
  if (!is.null(.label)) {
    if (.label == "latex" | .label == "webtex") {
      IC <- mutate(IC, R.nm = R_labeller(.ion1, .ion2, label = .label)) %>%
        select(!!! gr_by, R.nm = "R.nm", !!! ls_latex)
      return(IC)
      }
    }
  mutate(IC, R.nm = paste(.ion1, .ion2, sep = "/"))
}

#' Propagation of error for isotope ratios
#'
#' \code{stat_SDprop} function for propagation of descriptive and predictive
#' (Poisson) error statistics for isotope ratios (R).
#'
#' Isotope ratios are based on two measured variables, i.e., two isotopes of a
#' single chemical species. The combined variable compounds the uncertainty
#' associated with the individual measurement. Hence error propagation is
#' required to obtain a reliable estimate of the uncertainty associated
#' with the compounded variable (i.e. isotope ratio).
#'
#' @param ion1 A numeric vector constituting the single ion count rate of the
#' heavy isotope if argument predicted is set to \code{FALSE}, otherwise counts
#' are required.
#' @param ion2 A numeric vector constituting the single ion count rate of the
#' light isotope if argument predicted is set to \code{FALSE}, otherwise counts
#' are required.
#' @param type A character string for the type of uncertainty estimate:
#' \code{"sd"}, standard deviation; \code{"rsd"}, relative standard deviation in
#' per mille; \code{"se"}, standard error of the mean; \code{"rse"}, relative
#' standard error of the mean in per mille.
#' @param predicted Logical indicating whether sd is descriptive \code{FALSE},
#' or based on predicted Poisson value \code{TRUE}.
#'
#' @return A numeric vector containing the propagated uncertainty of the isotope
#' ratio.
#' @export
#' @examples
#
#' # Light isotope count rates
#' `32S` <- c(22318.19, 22151.20, 22429.52, 22143.78, 22574.25, 22455.50)
#'
#' # Heavy isotope count rate
#' `34S` <- c(231.4838, 220.3725, 255.5584, 237.0395, 244.4471, 238.8914)
#'
#' # Propagation
#' stat_SDprop(`34S`, `32S`, "rse")
stat_SDprop <- function(ion1, ion2, type = "sd", predicted = FALSE){

  M_ion1 <- mean(ion1) # heavy
  M_ion2 <- mean(ion2) # light
  M_R <- M_ion1 / M_ion2 # R
  n_ion1 <- sum(is.finite(ion1)) # observations heavy
  n_ion2 <- sum(is.finite(ion2)) # observations light
  if (n_ion1 != n_ion2) stop("Unequal number of measurements between isotopes.")

  if (predicted) {

    hat_sd <- sqrt((1 / sum(ion1)) + (1 / sum(ion2)))

    if (type == "sd") return(hat_sd * M_R * sqrt(n_ion1))
    if (type == "rsd") return(hat_sd * sqrt(n_ion1) * 1000)
    if (type == "se") return(hat_sd * M_R / sqrt(n_ion1))
    if (type == "rse") return((hat_sd / sqrt(n_ion1)) * 1000)

    } else {

      S_ion1 <- sd(ion1)
      S_ion2 <- sd(ion2)
      sd <- sqrt(
        ((S_ion1 / M_ion1) ^ 2) + ((S_ion2 / M_ion2) ^ 2) -
          (2 * (cov(ion1, ion2, method = "pearson", use = "everything") /
            (M_ion1 * M_ion2))
          )
        )
  if (type == "sd") return(sd * M_R)
  if (type == "rsd") return(sd * 1000)
  if (type == "se") return(sd * M_R / sqrt(n_ion1))
  if (type == "rse") return((sd / sqrt(n_ion1)) * 1000)

  }
}


#-------------------------------------------------------------------------------
# Helper functions for parsing, testing and validation
#-------------------------------------------------------------------------------

# Switch output complete dataset, stats or summary stats
mod_cal <- function(type, calcs) {
  switch(
    type,
    complete = call2( "mutate", expr(.), !!! calcs),
    stat = call2( "transmute", expr(.), !!! calcs),
    sum = call2("summarize", expr(.), !!! calcs),
  )
}

# Call specific to removal of zero count measurements with ZeroCt
zeroCt_cal <- function(zero_arg, .ion1, .ion2, gr_by, .N, .species, .t){

  if (zero_arg){
    call2(
      "zeroCt",
      expr(.),
      .ion1 = expr(.ion1),
      .ion2 = expr(.ion2),
      !!! gr_by,
      .N = expr(!! .N),
      .species = expr(!! .species),
      .t = expr(!! .t),
      .warn = TRUE
      )
    } else {
      call2("invisible", expr(.))
      }
}

# build new quosures and names for calcs
arg_builder <- function(args, stat){

  if (stat == "Xt") stat_names <- point::names_stat_Xt
  if (stat == "R") {
    stat_names <- mutate(
      point::names_stat_R,
      name = paste(name, "R", sep = "_")
      )
  }

  arg_names <- paste(stat_names$name, stat_names$origin, sep = "_")
  purrr::map2(
    stat_names$origin,
    stat_names$name,
    ~quo_updt(args[[paste0(".", .x)]], pre = .y)
    ) %>%
  set_names(nm = arg_names)
}
