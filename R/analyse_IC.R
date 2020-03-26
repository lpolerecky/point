#' Analyse raw ion count data
#'
#' \code{stat_Xt} function for descriptive and predictive statistics on single
#' ion precision.
#' \code{stat_R} function for descriptive and predictive statistics on isotope
#' ratios (R) precision with appropiate error propagation.
#'
#' These functions are a convenient wrapper to calculate the statistics
#' pertaining to the precision of pulsed ion count data (e.g. secondary ion mass
#' spectrometry). The statistics can either be calculated for single ions or
#' isotope ratios, and include observed and predicted (Poisson) statistics.
#' Calculations for isotope ratios include proper error propagation. For more
#' information on the usage as well as the mathematics behind these
#' functions see \code{vignette("IC-precision", package = "point")}.
#'
#' @param df A tibble containing processed ion count data.
#' @param Xt A variable constituting the ion count rate.
#' @param N A variable constituting the ion counts.
#' @param species A variable constituting the species analysed.
#' @param ... Variables for grouping.
#' @param ion1 A character string constituting the heavy isotope ("13C").
#' @param ion2 A character string constituting the light isotope ("12C").
#' @param latex A logical indicating whether variable names are latex
#' compatible.
#' @param output A character string for output as summary statistics ("sum");
#' statistics only ("stat"); and statistics with the original data ("complete").
#' @param zero A character string that determines whether analyses with zero
#' count measurements present will be removed from the calculations.
#'
#' @return A \code{\link[tibble:tibble]{tibble}} containing descriptive and
#' predictive statistics for ion counts and isotope ratios. The naming
#' convention depends on the argument \code{latex}; if set to \code{FALSE},
#' variable names concerning statistics will consist of an abbreviation pasted
#' together with the input variable names of \code{Xt}.
#' @export
#' @examples
#' # Use point_example() to access the examples bundled with this package
#'
#' # raw data containing 13C and 12C counts on carbonate
#' tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # Processing raw ion count data
#' tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt)
#'
#' # Single ion descriptive an predictive statistics for all measured ions
#' tb.Xt <- stat_Xt(tb.pr, Xt.pr, N.pr, species.nm, file.nm)
#'
#' # Descriptive an predictive statistics for 13C/12C ratios
#' tb.R <- stat_R(tb.pr, Xt.pr, N.pr, species.nm, ion1 = "13C",
#'                ion2 = "12C", file.nm)
stat_Xt <- function(df, Xt, N, species, ... , latex = FALSE, output = "sum"){

  stopifnot(tibble::is_tibble(df))

# Quoting the call  (user-supplied expressions)
  Xt <- enquo(Xt)
  N <- enquo(N)
  species <- enquo(species)
  gr_by <- enquos(...)

# The statistics
  args <- list(
    # number of measurements
    quo(n()),
    # sum of ion counts (important for external reproducibility calcs)
    quo(sum(!! N)),
    # mean ion count rate
    quo(mean(!! Xt)),
    # standard deviation (SD) count rate
    quo(sd(!! Xt)),
    # standard error of the mean (SE) count rate
    quo(sd(!! Xt) / sqrt(n())),
    # predicted SD count rate
    quo(sqrt(sum(!!N))),
    # predicted SE count rate
    quo(sqrt(sum(!!N) / n())))

# The statistic names (depend on user-supplied expression)
  ls.names <-paste0(c("n_", "Ntot_", "M_", "S_", "SeM_", "hat_S_", "hat_SeM_"),
                    as_name(Xt))

# To render nice latex variable names in Rmarkdown/Latex
  ls.latex <- set_names(ls.names, nm = c("$n$",
                                        "$N_{tot}$",
                                        "$\\bar{X}$",
                                        "$s_X$",
                                        "$s_\\bar{X}$",
                                        "$\\hat{s}_N$",
                                        "$\\hat{s}_\\bar{N}$"))

# Set statistic names
  args <- set_names(args, nm = ls.names)

# Switch output complete dataset, stats or summary stats
  mod_cal <- function(type) {
    switch(type,
           complete = call2( "mutate", expr(.), !!! args),
           stat = call2( "transmute", expr(.), !!! args),
           sum = call2("summarize", expr(.), !!! args),
          )
  }

# Grouping names
  tb_names  <- sapply(gr_by, as_name) %>% set_names()

# Evaluate expressions
  df <- df %>%
          group_by(!!! gr_by, !! species) %>%
          eval_tidy(expr = mod_cal(output), data = .) %>%
          ungroup()

# Output
  if (latex) {

      df <- df %>%
              mutate(!! species := purrr::map_chr(!! species,
                                                  latex_parser)) %>%
              select(!!! tb_names,
                     !!! ls.latex)

  return(df)

  }

  return(df)

  }

#' @rdname stat_Xt
#'
#' @export
stat_R <- function(df, Xt, N, species, ion1, ion2, ..., latex = FALSE,
                   output = "sum", zero = FALSE){

  stopifnot(tibble::is_tibble(df))

# Quoting the call (user-supplied expressions)
  Xt <- enquo(Xt)
  N <- enquo(N)
  species <- enquo(species)
  gr_by <- enquos(...)

# Update expressions (heavy isotope)
  Xt1 <- quo_updt(Xt, ion1) # count rate
  M.Xt1 <- quo_updt(Xt, ion1, "M") # mean count rate
  S.Xt1 <- quo_updt(Xt, ion1, "S") # sd count rate
  Yt1 <- quo_updt(N, ion1) # counts

# Update expressions (light isotope)
  Xt2 <- quo_updt(Xt, ion2) # count rate
  M.Xt2 <- quo_updt(Xt, ion2, "M") # mean count rate
  S.Xt2 <- quo_updt(Xt, ion2, "S") # sd count rate
  Yt2 <- quo_updt(N, ion2) # counts

# The statistics
  args <- list(
    # number of measurements
    quo(n()),
    # mean isotope ratio
    quo(mean(!!Xt1) / mean(!!Xt2)),
    # SD isotope ratio
    quo(stat_SDprop(M_R = !!quo_updt(my_q = Xt , x = "M_R"),
                    ion1 = !!Xt1,
                    ion2 = !!Xt2,
                    M_ion1 = !!M.Xt1,
                    M_ion2 = !!M.Xt2,
                    S_ion1 = !!S.Xt1,
                    S_ion2 = !!S.Xt2,
                    n = !!quo_updt(my_q = Xt , x = "n_R"),
                    type = "sd")),
    # RSD isotope ratio
    quo((!!quo_updt(my_q = Xt , x = "S_R") /
         !!quo_updt(my_q = Xt , x = "M_R")) * 1000),
    # SE isotope ratio
    quo(!!quo_updt(my_q = Xt , x = "S_R") / sqrt(n())),
    # RSE isotope ratio
    quo((!!quo_updt(my_q = Xt , x = "SeM_R") /
         !!quo_updt(my_q = Xt , x = "M_R")) * 1000),
    # predictive SD isotope ratio
    quo(stat_hSDprop(M_R = !!quo_updt(my_q = Xt , x = "M_R"),
                    ion1 = !!Xt1,
                    ion2 = !!Xt2,
                    N_ion1 = !!Yt1,
                    N_ion2 = !!Yt2,
                    n = !!quo_updt(my_q = Xt , x = "n_R"),
                    type = "sd")),
    # predictive RSD isotope ratio
    quo((!!quo_updt(my_q = Xt , x = "hat_S_R") /
         !!quo_updt(my_q = Xt , x = "M_R")) * 1000),
    # predictive SE isotope ratio
    quo(!!quo_updt(my_q = Xt , x = "hat_S_R") / sqrt(n())),
    # predictive RSE isotope ratio
    quo((!!quo_updt(my_q = Xt , x = "hat_SeM_R") /
           !!quo_updt(my_q = Xt , x = "M_R")) * 1000),
    # reduced chi squared
    quo((!!quo_updt(my_q = Xt , x = "SeM_R") /
         !!quo_updt(my_q = Xt , x = "hat_SeM_R")) ^ 2)
  )

# The statistic names (depend on user-supplied expression)
  ls.names <-paste(c("n", "M", "S", "RS", "SeM", "RSeM",
                     "hat_S", "hat_RS", "hat_SeM", "hat_RSeM",
                     "chi2"), "R", as_name(Xt), sep = "_")

# To render nice latex variable names in Rmarkdown/Latex
  nm <- c("$n$",
           "$\\bar{R}$",
           "$s_{R}$",
           "$\\epsilon_{R} \\,$ (\u2030)",
           "$s_{\\bar{R}}$",
           "$\\epsilon_{\\bar{R}} \\,$ (\u2030)",
           "$\\hat{s}_{R}$",
           "$\\hat{\\epsilon}_{R} \\,$ (\u2030)",
           "$\\hat{s}_{\\bar{R}}$",
           "$\\hat{\\epsilon}_{\\bar{R}} \\,$ (\u2030)",
           "$\\chi^{2}$")

  ls.latex <- set_names(ls.names, nm = nm)

# Set statistic names
  args <- set_names(args, nm = ls.names)

# Extra arg in case of mutate and transmute (isotope ratio for each time step)
  args2 <- append(args,
                  set_names(list(quo(!!Xt1 / !!Xt2)),
                            nm = paste("R", as_name(Xt), sep = "_")))

# Switch output complete dataset, stats or summary stats
  mod_cal <- function(type) {
    switch(type,
           complete = call2( "mutate", expr(.), !!! args2),
           stat = call2( "transmute",  expr(.), !!! args2),
           sum = call2("summarize", expr(.),  !!! args),
           )
  }

# Call specific to removal of zero count measurements with ZeroCt
  zeroCt_cal <- function(zero_arg){

    if (zero_arg){

      call2("zeroCt",
            expr(.),
            N = expr(!! N),
            species = expr(!! species),
            ion1 = expr(ion1),
            ion2 = expr(ion2),
            !!! gr_by)

      } else {

        call2("invisible", expr(.))

      }
  }

# Grouping names
  tb_names  <- sapply(gr_by, as_name) %>% set_names()


# Evaluate expressions and calls
  df <- df  %>%
# Calculate single ion count stats
          eval_tidy(expr = zeroCt_cal(zero)) %>%
          stat_Xt(Xt = !!Xt, N = !!N, species =  !!species,
                  !!! gr_by, output = "complete") %>%
          cov_R(species = !!species,
                ion1 = ion1, ion2 = ion2, !!! gr_by,
                preserve = TRUE) %>%
          group_by(!!! gr_by) %>%
          eval_tidy(expr = mod_cal(output)) %>%
          ungroup()

# Output
  if (latex) {

    df <- df %>%
            mutate(R.nm = latex_parser(ion1, ion2)) %>%
            select(!!! tb_names,
                   R = "R.nm",
                   !!! ls.latex)

    return(df)

  } else {

    df <- df %>%
            mutate(R.nm = paste(ion1, ion2, sep = "/"))

    return(df)

  }
}

#' Propagation of error for isotope ratios
#'
#' \code{stat_SDprop} function for propagation of descriptive error statistics
#' for isotope ratios (R). \code{stat_hSDprop} function for propagation of
#' predictive (Poisson) error statistics for isotope ratios (R).
#'
#' Isotope ratios are based on two measured variables, i.e., two isotopes of a
#' single chemical species. The combined variable compounds the uncertainty
#' associated with the individual measurement. Hence error propagation is
#' required to obtain a reliable estimate of the uncertainty associated
#' with the compounded variable (i.e. isotope ratio).
#'
#' @param M_R A numeric vector constituting the mean isotope ratio.
#' @param ion1 A numeric vector constituting the single ion count rate of the
#' heavy isotope.
#' @param ion2 A numeric vector constituting the single ion count rate of the
#' light isotope
#' @param N_ion1 A numeric vector constituting the single ion counts of the
#' heavy isotope
#' @param N_ion2 A numeric vector constituting the single ion counts of the
#' light isotope.
#' @param M_ion1 A numeric vector constituting the mean ion count of the heavy
#' isotope.
#' @param M_ion2 A numeric vector constituting the mean ion count of the light
#' isotope.
#' @param S_ion1 A numeric vector constituting the standard deviation of the ion
#' counts of the heavy isotope.
#' @param S_ion2 A numeric vector constituting the standard deviation of the ion
#' counts of the light isotope.
#' @param n A numeric vector constituting the length of the count block
#' \emph{n}.
#' @param type A character string for the type of uncertainty estimate:
#' \code{"sd"}, standard deviation; \code{"rsd"}, relative standard deviation in
#' per mille; \code{"se"}, standard error of the mean; \code{"rse"}, relative
#' standard error of the mean in per mille.
#'
#' @return A numeric vector containing the propagated uncertainty of the isotope
#' ratio.
#' @export
#' @examples
#
#' # Light isotope count rates
#' `32S` <- c(22318.19, 22151.20, 22429.52, 22143.78, 22574.25, 22455.50,
#'          22516.73, 22414.68, 22288.50, 22327.47)
#'
#' # Heavy isotope count rate
#' `34S` <- c(231.4838, 220.3725, 255.5584, 237.0395, 244.4471, 238.8914,
#'          238.8914, 237.0395, 262.9660, 264.8179)
#'
#' # Mean 32S
#' M_32S <- mean(`32S`)
#'
#' # Mean 34S
#' M_34S <- mean(`34S`)
#'
#' # Mean R
#' M_R <- M_34S / M_32S
#'
#' # Standard deviation 32S
#' S_32S <- sd(`32S`)
#'
#' # Standard deviation 34S
#' S_34S <- sd(`34S`)
#'
#' # Propagation
#' stat_SDprop(M_R, `34S`, `32S`, M_34S, M_32S, S_34S, S_32S, 10, "rse")
stat_SDprop <- function(M_R, ion1, ion2, M_ion1, M_ion2, S_ion1, S_ion2, n,
                        type = "sd"){

  M_R <- unique(M_R)
  n <- unique(n)

  sd <- sqrt(
              ((unique(S_ion1) / unique(M_ion1)) ^ 2) +
              ((unique(S_ion2) / unique(M_ion2)) ^ 2) -
              (2 *
                  (
                    cov(ion1, ion2, method = "pearson", use = "everything") /
                    (unique(M_ion1) * unique(M_ion2))
                  )
              )
            )

  if (type == "sd") {return(sd * M_R)}
  if (type == "rsd") {return(sd * 1000)}
  if (type == "se") {return(sd * M_R / sqrt(n))}
  if (type == "rse") {return((sd / sqrt(n)) * 1000)}
}

#' @rdname  stat_SDprop
#'
#' @export
stat_hSDprop <- function(M_R, ion1, ion2 ,N_ion1, N_ion2, n, type = "sd"){

  M_R <- unique(M_R)
  n <- unique(n)

  hat_sd <- sqrt(
                  (1 / sum(N_ion1)) +
                  (1 / sum(N_ion2))
                )

  if (type == "sd") {return(hat_sd * M_R * sqrt(n))}
  if (type == "rsd") {return(hat_sd * sqrt(n) * 1000)}
  if (type == "se") {return(hat_sd * M_R / sqrt(n))}
  if (type == "rse") {return((hat_sd / sqrt(n)) * 1000)}
}

#' Remove analytical runs with zero counts
#'
#' \code{zeroCt} removes analytical runs for isotope ratios that contain zero
#' counts.
#'
#' This functions removes analytical runs with zero counts for calculating
#' isotope ratios.
#'
#' @param df A tibble containing processed ion count data.
#' @param N A variable constituting the ion counts.
#' @param species A variable constituting the species analysed.
#' @param ... Variables for grouping.
#' @param ion1 A character string constituting the heavy isotope ("13C").
#' @param ion2 A character string constituting the light isotope ("12C").
#' @param warn A logical indicating whether to produce a warning.
#'
#' @return A \code{\link[tibble:tibble]{tibble}} containing the single ion count
#' dataset for the specified ion ratio. The grouping variable specifies on
#' which level the original dataset is subsetted.
#' @export
#' @examples
#' # Use point_example() to access the examples bundled with this package
#'
#' # raw data containing 13C and 12C counts on carbonate
#' tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # Remove analyses with zero counts
#' # Vectors of isotope ratios
#' ion1 <-  c("13C", "12C 13C", "13C 14N", "12C 14N", "12C")
#' ion2 <-  c("12C", "12C2", "12C 14N", "40Ca 16O", "40Ca 16O")
#' tb.pr <- purrr::map2(ion1, ion2,
#'          ~zeroCt(tb.rw, N.rw, species.nm, .x, .y, file.nm))
zeroCt <- function(df, N, species, ion1, ion2, ..., warn = TRUE){

  N <- enquo(N)
  species <- enquo(species)
  gr_by <- enquos(...)

  df <- df %>%
    filter(!!species == ion1 | !!species == ion2)

  if (warn){
    if (any(df %>% select(!!N) %>% pull(!!N) == 0)) {

      warning("Zero counts present and removed")

    }
  }

  ls.zero <- df %>%
    filter(!! N == 0) %>%
    select(!!! gr_by)

  df.pr <- anti_join(df, ls.zero, by = sapply(gr_by, as_name))

  return(df.pr)

}

#' Latex parsable chemical species names
#'
#' \code{latex_parser} Converts a character string containing chemical species
#' names in a latex parsable string.
#'
#' This functions converts chemical species names of the form, e.g. `"12C"`,
#' `"13C2"`, or `"12C 14N"`, to a character string which can be parsed in Latex
#' to species names with appropiate superscripts on th left for mass and
#' subscripts for stochiometry on the right. If both arguments are given; e.g.,
#' `ion1 = "13C"` and `ion2 = "12C"`, the function will create an isotope ratio
#' seperate.
#'
#' @param ion1 A chemical species character string (for isotope ratios the heavy
#' isotope).
#' @param ion2 A chemical species character string (for isotope ratios the light
#' isotope).
#'
#' @return A character string parsable in Latex
#' @export
latex_parser <- function(ion1, ion2 = NULL){

  lat_pars <- function(ion){

    # Switch str extract method: numeric sub- or superscript or character
    mod_cal <- function(type) {
      switch(type,
             sup = call2("str_extract",
                                quote(ion),
                                "[:digit:]+(?=[:alpha:])"),
             char = call2("str_extract",
                                 quote(ion),
                                 "[:alpha:]+"),
             sub = rlang::call2("str_extract",
                                quote(ion),
                                "(?<=[:alpha:])[:digit:]+")
      )
    }

    # Call substitute to change element in latex expr
    substitute(paste0("$\\phantom{,}^{", sup, "}$", x, "$_{", sub, "}$"),
               list(sup = eval(mod_cal("sup")),
                    x = eval(mod_cal("char")),
                    sub = eval(mod_cal("sub")))) %>%
      eval()
  }

  latex_poly <- function(ion){
    # Seperation for strings for poly atomic species
    ls.str <- as.list(unlist(str_split(ion, "[:blank:]")))

    purrr::map(ls.str, lat_pars) %>% purrr::reduce(paste, sep = "--")

  }

  if(!is.null(ion2)){

    ls.ions <- lst(ion1 , ion2)
    lat <- purrr::map(ls.ions, latex_poly)  %>% purrr::reduce(paste, sep = "/")

  }else{

    if (any(str_detect(ion1, "[:blank:]"))){
      lat <- latex_poly(ion1)

    }else{

      lat <- lat_pars(ion1)

    }
  }

  # Replace NAs from string with empty
  lat <- str_replace_all(lat, "NA", "")
  return(lat)

}

#-------------------------------------------------------------------------------
# Helper functions for parsing, testing and validation
#-------------------------------------------------------------------------------
# Function for covariate convertion of isotope systems
cov_R <- function(df, species, ion1, ion2, ..., preserve = FALSE){

  gr_by <- enquos(...)
  species <- enquo(species)

# CAll that creates an ID Uniqually identifies ion pairs for calculating isotope
# ratios (in case the data file does not contain an ID)
  ID_cal <- function(ID_arg){

    if (ID_arg){

      call2("ID_builder",
            expr(.),
            species = expr(!! species),
            !!! gr_by)

    } else {

      call2("invisible", expr(.))

    }
  }

  df <- df %>%
    eval_tidy(expr = ID_cal(df %>%
                       select(starts_with("ID")) %>%
                       length() == 0))

# Filtering single ion stats minor isotope
  df.13C <- df %>%
    filter(!! species == ion1)

# Filtering single ion stats major isotope
  df.12C <- df %>%
    filter(!! species == ion2) %>%
    select(-c(!!! gr_by))

# When adding the labels for the ion blanks are removed to match quo_updt
  df.R <- full_join(df.13C,  df.12C, by = "ID",
                    suffix = c(paste0(".", str_replace_all(ion1, " ", "")),
                               paste0(".", str_replace_all(ion2, " ", ""))))

  if (preserve){return(df.R)}else{return(select(df.R, -.data$ID))}
  }



# Function which updates quosures for subsequent tidy evaluation
quo_updt <- function(my_q, txt = NULL, x = NULL, sepfun = "_"){

# Get expressions
  old_expr <- get_expr(my_q)
# Modify expression, turn expr into a character string
  if (length(txt) == 0){

    new_chr <-paste(x, expr_text(old_expr), sep = sepfun)

  }else{

    if (length(x) == 0){

      new_chr <- paste(expr_text(old_expr), txt, sep = ".")

    }else{

      new_chr <- paste(paste(x, expr_text(old_expr), sep = sepfun), txt, sep = ".")

      }
  }
# New expression from character (remove whitespace)
  new_expr <- parse_expr(str_replace_all(new_chr, " ", ""))
# Update old quosure
  set_expr(my_q, new_expr)

  }



# Function for summation error propagation (will be updated later on)
sum_sd_prop <- function(x, type = "rse"){

  n <- length(x)
  sum_sd <- sqrt( 1 / sum(x^-2))

}
