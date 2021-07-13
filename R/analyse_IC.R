#' Analyse raw ion count data
#'
#' \code{stat_X} function for descriptive and predictive statistics on single
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
#' @param .nest A variable hat identifies a series of analyses to calculate
#' external precision.
#' @param .X A variable constituting the ion count rate (defaults to
#' variables generated with \code{read_IC()}.)
#' @param .N A variable constituting the ion counts (defaults to variables
#' generated with \code{read_IC()}.).
#' @param .species A variable constituting the species analysed (defaults to
#' variables generated with \code{read_IC()}.).
#' @param .t A variable constituting the time of the analyses (defaults to
#' variables generated with \code{read_IC()}.).
#' @param .stat Select statistics (e.g. \code{c("M", "RS"}), see the tables
#' \code{point::names_stat_X} and \code{point::names_stat_R} for the full
#' selection of statistics available (default uses all statistic
#' transformations).
#' @param .label A character string indicating whether variable names are latex
#' (\code{"latex"}) or webtex (\code{"webtex"}) compatible. Will be extended in
#' the future \code{default = NULL}.
#' @param .meta Logical whether to preserve the metadata as an attribute
#' (defaults to TRUE).
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
#' tb_rw <- read_IC(point_example("2018-01-19-GLENDON"), meta = TRUE)
#'
#' # Processing raw ion count data
#' tb_pr <- cor_IC(tb_rw)
#'
#' # Single ion descriptive an predictive statistics for all measured ions
#' stat_X(tb_pr, file.nm)
#'
#' # Descriptive an predictive statistics for 13C/12C ratios
#' stat_R(tb_pr, "13C", "12C", file.nm, .zero = TRUE)
#'
#' # Descriptive an predictive statistics for 13C/12C ratios (external)
#' stat_R(tb_pr, "13C", "12C", sample.nm, file.nm, .nest = file.nm,
#'        .zero = TRUE)
stat_X <- function(.IC, ..., .X = NULL, .N =  NULL, .species = NULL,
                   .t =  NULL, .stat = point::names_stat_X$name,
                   .label = "none", .meta = FALSE, .output = "sum") {

  # Checks
  stat_validator(.IC, .stat, point::names_stat_X$name, .output, .label)

  # Metadata
  if(.meta) meta <- unfold(.IC, merge = FALSE)

  # Quoting the call (user-supplied expressions) and complete if NULL
  args <- inject_args(
    .IC,
    enquos(.X = .X, .N = .N, .species = .species, .t = .t),
    type = c("processed", "group")
    )

  # Argument check
  argument_check(.IC, args, "processed")

  # Grouping
  gr_by <- enquos(...)

  # New quosures
  new_args <- arg_builder(args, "X")
  args <- append(args, new_args)

  # The statistics
  calcs <- quos(
    # number of measurements
    n(),
    # sum of ion counts (important for external reproducibility calcs)
    sum(!! args[[".N"]]),
    # mean ion count rate
    mean(!! args[[".X"]]),
    # standard deviation (SD) count rate
    sd(!! args[[".X"]]),
    # RSD ion count rate
    (!! args[["S_X"]] / !! args[["M_X"]]) * 100,
    # standard error of the mean (SE) count rate
    sd(!! args[[".X"]]) / sqrt(n()),
    # predicted SD count rate
    sqrt(mean(!! args[[".N"]])),
    # predicted RSD count rate
    (1 / sqrt(mean(!! args[[".N"]]))) * 100,
    # predicted SE count rate
    sqrt(mean(!! args[[".N"]]) / n()),
    # reduced chi squared
    (!! args[["SeM_X"]] / !! args[["hat_SeM_N"]]) ^ 2
    )

  # The statistic names (depend on user-supplied expression)
  ls_nm <- sapply(new_args, as_name)

  # Set statistic names
  calcs <- set_names(calcs, nm = ls_nm)

  # Stat selection
  vars <- stat_selector(.stat, ls_nm)

  # Render latex variable names
  if (.label == "latex" | .label == "webtex") {
    ls_latex <- set_names(
      vars$pos,
      tex_labeller(point::names_stat_X, .stat, .label)
      )
    }

  # Evaluate expressions
  IC <- group_by(.IC, !!! gr_by, !! args[[".species"]]) %>%
    eval_tidy(expr = call2(mod_call[[.output]], expr(.), !!! calcs)) %>%
    ungroup() %>%
    select(-any_of(vars$neg))

  # Return metadata
  if (.meta) IC <- fold(IC, type = ".mt",  meta = meta)

  # Output
  if (.label == "latex" | .label == "webtex") {
    IC <- mutate(
      IC,
      !! args[[".species"]] :=
        purrr::map_chr(!! args[[".species"]], ion_labeller, label = .label)
      ) %>%
      select(!!! gr_by, !! args[[".species"]], !!! ls_latex)
    }
  IC
  }

#' @rdname stat_X
#'
#' @export
stat_R <- function(.IC, .ion1, .ion2, ..., .nest = NULL, .X = NULL, .N = NULL,
                   .species = NULL, .t = NULL,
                   .stat = point::names_stat_R$name, .label = "none",
                   .output = "sum", .zero = FALSE){

  # Checks
  stat_validator(.IC, .stat, point::names_stat_R$name, .output, .label)

  # Quoting the call (user-supplied expressions)
  args <- inject_args(
    .IC,
    enquos(.X = .X, .N = .N, .species = .species, .t = .t),
    type = c("processed", "group")
    )

  # Argument check
  argument_check(.IC, args, "processed")

  # Grouping
  gr_by <- enquos(...)
  # Nesting
  nest <- enquo(.nest)

  # External precision
  if (!is.null(quo_get_expr(nest))) {

    # Calculate single ion stat to obtain mean and total ion counts
    .IC <- call2("stat_X", .IC, !!! gr_by, !!! args, .ns = "point") %>%
      eval()

    # Updated quotes
    args[[".X"]] <- quo_updt(args[[".X"]], pre = "M")
    args[[".N"]] <- quo_updt(args[[".N"]], pre = "tot")
    # Updated grouping
    gr_by <- gr_by[!sapply(gr_by, as_name) %in% as_name(nest)]

    # Check whether groups have enough observations
    .IC <- add_count(.IC, !!!  gr_by, !! args[[".species"]])
    if (any(.IC$n <= 1)) {
      warning(
        "Some groups have too few observations for a reliable estimation of the external precision. These groups have been omitted.",
         call. = FALSE
         )
      .IC <- filter(.IC, n > 1)
    }
    nest <- TRUE
    }

  # If t column is empty create a manual time increment (important for external
  # precision)
  if (!(as_name(args[[".t"]]) %in% colnames(.IC))) {
    .IC <- group_by(.IC, !!! gr_by, !! args[[".species"]]) %>%
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
    X1 = quo_updt(args[[".X"]], post = .ion1), # count rate
    N1 = quo_updt(args[[".N"]], post = .ion1), # counts
    # Update quosures (common isotope)
    X2 = quo_updt(args[[".X"]], post = .ion2), # count rate
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
    mean(!! args[["X1"]]) / mean(!! args[["X2"]]),
    # SD isotope ratio
    stat_SDprop(!! args[["X1"]], !! args[["X2"]]),
    # RSD isotope ratio
    (!! args[["S_R"]] / !! args[["M_R"]]) * 1000,
    # SE isotope ratio
    !! args[["S_R"]] / sqrt(n()),
    # RSE isotope ratio
    (!! args[["SeM_R"]] / !! args[["M_R"]]) * 1000,
    # predictive SD isotope ratio
    stat_SDprop(!! args[["N1"]], !! args[["N2"]], predicted = TRUE),
    # predictive RSD isotope ratio
    !! args[["hat_S_R"]] / !! args[["M_R"]] * 1000,
    # predictive SE isotope ratio
    !! args[["hat_S_R"]]  / sqrt(n()),
    # predictive RSE isotope ratio
    !! args[["hat_SeM_R"]]  / !! args[["M_R"]] * 1000,
    # reduced chi squared
    (!! args[["SeM_R"]] / !! args[["hat_SeM_R"]]) ^ 2
    )

  # The statistic names (depend on user-supplied expression)
  ls_nm <- sapply(new_args, as_name)

  # Set statistic names
  calcs <- set_names(calcs, nm = ls_nm)

  # Extra arg in case of mutate and transmute (isotope ratio for each time step)
  if (.output != "sum") {
    R <- paste("R", as_name(args[[".X"]]), sep = "_")
    calcs[[R]] <- quo(!! args[["X1"]] / !! args[["X2"]])
    }

  # Stat selection
  vars <- stat_selector(.stat, ls_nm)

  # Render latex variable names
  if (.label == "latex" | .label == "webtex") {
    tb_tex <- point::names_stat_R
    if (isTRUE(nest)) {
      tb_tex <- mutate(
        tb_tex,
        origin = case_when(origin == "X" ~ "M", origin == "N" ~ "Ntot")
        )
    }
    ls_latex <- set_names(
      vars$pos,
      tex_labeller(tb_tex, .stat, .label)
      )
    }

  # Evaluate expressions and calls
  IC <- eval(zeroCt_call(.zero, .IC, .ion1, .ion2, gr_by, args)) %>%
    cov_R(c(.ion1, .ion2), !!! gr_by, .species = !! args[[".species"]],
          .t = !! args[[".t"]]) %>%
    group_by(!!! gr_by) %>%
    eval_tidy(expr = call2(mod_call[[.output]], expr(.), !!! calcs)) %>%
    ungroup() %>%
    select(-any_of(vars$neg))

  # Output
  if (.label == "latex" | .label == "webtex") {
    IC <- mutate(IC, ratio.nm = R_labeller(.ion1, .ion2, label = .label)) %>%
      select(!!! gr_by, .data$ratio.nm, !!! ls_latex)
    }
  mutate(IC, ratio.nm = paste(.ion1, .ion2, sep = "/"))
}

#' Propagation of errors for isotope ratios
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

  if (isTRUE(predicted)) {
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
mod_call <-  c(complete = "mutate", stat = "transmute", sum = "summarize")

# Call specific to removal of zero count measurements with ZeroCt
zeroCt_call <- function(zero, IC, .ion1, .ion2, gr_by, args) {
  if (isTRUE(zero)) {
    call2(
      "zeroCt",
      IC,
      .ion1 = .ion1,
      .ion2 = .ion2,!!!gr_by,
      .N = args[[".N"]],
      .species = args[[".species"]],
      .warn = TRUE,
      .ns = "point"
    )
  } else {
    call2("invisible", IC)
  }
}

# Build new quosures and names for calcs
arg_builder <- function(args, stat, ion = NULL, append = NULL) {

  if (stat == "X") arg_names <- point::names_stat_X
  if (stat == "R") arg_names <- point::names_stat_R
  if (stat == "model") {
    arg_names <- point::names_model
    pre <-  NULL
  } else {
    pre  <- "."
  }

  # no origin of variable names
  if (!"origin" %in% colnames(arg_names)) arg_names$origin <- NA_character_

  arg_names <- mutate(
    arg_names,
    origin = if_else(is.na(.data$origin), .data$derived, .data$origin),
    label =
      if_else(
        .data$origin == .data$derived,
        paste0(paste(.data$name, .data$origin, sep = "_"), append),
        paste0(paste(.data$name, .data$derived, sep = "_"), append)
        ),
    name =
      if_else(
        .data$origin == .data$derived,
        .data$name,
        paste(.data$name, .data$derived, sep = "_")
        )
    )

  # quosure update
  args <- purrr::map2(
    arg_names$origin,
    arg_names$name,
    ~quo_updt(args[[paste0(pre, .x)]], pre = .y)
    )
  # wide format with ions
  if (!is.null(ion)) args <- purrr::map(args, ~quo_updt(.x, post = ion))
  # set names
  set_names(args, nm = arg_names$label)
}

# latex labeller function
tex_labeller <- function(vars, stat, label) {
  if (!"origin" %in% colnames(vars)) vars$origin <- vars$derived
  names_vars <- filter(vars, .data$name %in% stat) %>%
    # if variable has a stat component
    mutate(
      derived =
        if_else(
          stringr::str_detect(.data$derived, "[[:punct:]]"),
          stringr::str_extract("M_R", "(?<=[[:punct:]])[[:alpha:]]"),
          .data$derived
          )
       )
  purrr::pmap_chr(
    list(
      var = names_vars$derived,
      org = names_vars$origin,
      stat = names_vars$name
      ),
    stat_labeller,
    label = label
    )
}

# consistency checks
stat_validator <- function(IC, stat_in = NULL, stat_def = NULL,
                           output = "sum", label = "none") {

  if (!tibble::is_tibble(IC)) {
    stop("Ion count dataset should be a tibble object.", call. = FALSE)
  }
  if (!all(stat_in %in% stat_def)) {
    stop("Unkown statistic.", call. = FALSE)
  }
  if (output != "sum"  & label != "none") {
    stop("Latex labels is not supported for complete datasets.", call. = FALSE)
  }
}

# Stat selection function
stat_selector <- function(stat, vars) {
  str_stat <- stringr::str_c("^", stat, "_", collapse = "|")
  ls_vars <- list()
  ls_vars$pos <- purrr::keep(vars, ~stringr::str_detect(., str_stat))
  ls_vars$neg <- purrr::discard(vars, ~stringr::str_detect(., str_stat))
  ls_vars
}
