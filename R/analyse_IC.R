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
#' functions see the vignette \code{"IC-precision"}.
#'
#' @param df A tibble containing processed ion count data
#' @param Xt A variable constituting the ion count rate
#' @param N A variable constituting the ion counts
#' @param species A variable constituting the species analysed
#' @param ... Variables for grouping
#' @param ion1 A character string constituting the heavy isotope ("13C")
#' @param ion2 A character string constituting the light isotope ("12C")
#' @param latex A logical indicating whether variable names are latex compatible
#' @param output A character string for output as summary statistics ("sum");
#' statistics only ("stat"); and statistics with the original data ("complete")
#'
#' @return A t\code{\link[tibble:tibble]{tibble}} containing descriptive and
#' predictive statistics for ion counts and isotope ratios
#'
#' @examples
#' # Use point_example() to access the examples bundled with this package in the
#' # inst/extdata directory.
#'
#' # raw data containing 13C and 12C counts on carbonate
#' tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # processing raw ion count data
#' tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt)
#'
#' # single ion descriptive an predictive statistics for all measured ions
#' tb.Xt <- stat_Xt(tb.pr, Xt.pr, N.pr, species.nm, file.nm)
#'
#' # descriptive an predictive statistics for 13C/12C ratios
#' tb.R <- stat_R(tb.pr, Xt.pr, N.pr, species.nm, ion1 = "13C",
#'                ion2 = "12C", file.nm)
#'
#' @export
stat_Xt <- function(df, Xt, N, species, ... , latex = FALSE, output = "sum"){

# Quoting the call  (user-supplied expressions)
  Xt <- enquo(Xt)
  N <- enquo(N)
  gr_by <- enquos(..., species)

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
                    quo_name(Xt))

# To render nice latex variable names in Rmarkdown/Latex
  ls.latex <- ls.names %>%
    purrr::set_names(., nm = c("$n$",
                               "$N_{tot}$",
                               "$\\bar{X}$",
                               "$s_X$",
                               "$s_\\bar{X}$",
                               "$\\hat{s}_N$",
                               "$\\hat{s}_\\bar{N}$"))

# Set statistic names
  args <- purrr::set_names(args, nm = ls.names)

# Switch output complete dataset, stats or summary stats
  mod_cal <- function(type) {
    switch(type,
           complete = call2( "mutate", quote(.), quote(!!! args)),
           stat = call2( "transmute", quote(.), quote(!!! args)),
           sum = call2("summarize", quote(.), quote(!!! args)),
          )
  }

# Grouping names
  tb_names  <- sapply(gr_by, as_name) %>% purrr::set_names()

# Evaluate expressions
  df <- df %>%
              group_by(!!! gr_by) %>%
              eval_tidy(expr = mod_cal(output), data = .) %>%
              #summarise(!!! args) %>%
              ungroup()

# Output
  if (latex) {

      df <- df %>%
              mutate(species.nm = purrr::map_chr(species.nm,
                                                 ~latex_parser(.x,
                                                               ion1 = NULL,
                                                               ion2 = NULL))) %>%
              select(!!! tb_names,
                     !!! ls.latex)

  return(df)

  }

  return(df)

  }

#' @rdname stat_Xt
#'
#' @export
stat_R <- function(df, Xt, N, species, ion1, ion2, ..., latex = FALSE, output = "sum", zero = FALSE){

# Quoting the call  (user-supplied expressions)
  Xt <- enquo(Xt)
  N <- enquo(N)
  species <- enquo(species)
  gr_by <- enquos(...)

# Update expressions (heavy isotope)
  Xt1 <- quo_updt(Xt, ion1) # count rate
  M.Xt1 <- quo_updt(Xt, ion1, "M") # mean count rate
  S.Xt1 <- quo_updt(Xt, ion1, "S") # sd count rate
  Yt1 <- quo_updt(N, ion1) #counts

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
    quo(des_SD_prop(M_R = !!quo_updt(my_q = Xt , x = "M_R"),
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
    quo(hat_SD_prop(M_R = !!quo_updt(my_q = Xt , x = "M_R"),
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
                     "chi2"), "R", quo_name(Xt), sep = "_")

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

  ls.latex <- ls.names %>%
                purrr::set_names(nm = nm)

# Set statistic names
  args <- purrr::set_names(args, nm = ls.names)

# Extra arg in case of mutate and transmute (isotope ratio for each time step)
  args2 <- append(args, purrr::set_names(list(quo(!!Xt1 / !!Xt2)) ,
                                         nm = paste("R", quo_name(Xt),
                                                    sep = "_")))

# Switch output complete dataset, stats or summary stats
  mod_cal <- function(type) {
    switch(type,
           complete = call2( "mutate", quote(.), quote(!!! args2)),
           stat = call2( "transmute", quote(.), quote(!!! args2)),
           sum = call2("summarize", quote(.), quote(!!! args)),
           )
  }

# Grouping names
  tb_names  <- sapply(gr_by, as_name) %>% purrr::set_names()

# Evaluate expressions
  df <- df  %>%
# Calculate single ion count stats
          stat_Xt(., Xt = !!Xt, N = !!N, species =  !!species,
                  !!! gr_by, output = "complete")

# Check and remove zero counts
  if (zero == TRUE) {

    df <- df %>%
            zeroCt(., N = !!N, species = !!species,
                   ion1 = ion1, ion2 = ion2, !!! gr_by)
  }

# Isotopes as rowwise observations
  df <- df %>%
          cov_R(. , species = !!species,
                ion1 = ion1, ion2 = ion2, !!! gr_by) %>%
          group_by(!!! gr_by) %>%
          eval_tidy(expr = mod_cal(output), data = .) %>%
          ungroup()

# Output
  if (latex) {

    df <- df %>%
            mutate(R.nm = latex_parser(species = NULL, ion1, ion2)) %>%
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
#' \code{des_SD_prop} function for propagation of descriptive error statistics
#' for isotope ratios (R). \code{hat_SD_prop} function for propagation of
#' predictive (Poisson) error statistics for isotope ratios (R).
#'
#' Isotope ratios are based on two measured variables, i.e., two isotopes of a
#' single chemical species. The combined variable compounds the uncertainty
#' associated with the individual measurement. Hence error propagation is
#' required to obtain a reliable estimate of the uncertainty associated
#' with the compounded variable (i.e. isotope ratio).
#'
#' @param M_R A numeric vector constituting the mean isotope ratio
#' @param ion1 A numeric vector constituting the single ion count rate of the
#' heavy isotope
#' @param ion2 A numeric vector constituting the single ion count rate of the
#' light isotope
#' @param N_ion1 A numeric vector constituting the single ion counts of the
#' heavy isotope
#' @param N_ion2 A numeric vector constituting the single ion counts of the
#' light isotope
#' @param M_ion1 A numeric vector constituting the mean ion count of the heavy
#' isotope
#' @param M_ion2 A numeric vector constituting the mean ion count of the light
#' isotope
#' @param S_ion1 A numeric vector constituting the standard deviation of the ion
#' counts of the heavy isotope
#' @param S_ion2 A numeric vector constituting the standard deviation of the ion
#' counts of the light isotope
#' @param n A numeric vector constituting the length of the count block \emph{n}
#' @param type A character string for the type of uncertainty estimate:
#' \code{"sd"}, standard deviation; \code{"rsd"}, relative standard deviation in
#' per mille; \code{"se"}, standard error of the mean; \code{"rse"}, relative
#' standard error of the mean in per mille
#'
#' @return A numeric vector containing the propagated uncertainty of the isotope
#' ratio
#'
#' @examples
#
#' # light isotope count rates
#' `32S` <- c(22318.19, 22151.20, 22429.52, 22143.78, 22574.25, 22455.50,
#'          22516.73, 22414.68, 22288.50, 22327.47)
#'
#' # heavy isotope count rate
#' `34S` <- c(231.4838, 220.3725, 255.5584, 237.0395, 244.4471, 238.8914,
#'          238.8914, 237.0395, 262.9660, 264.8179)
#'
#' # mean 32S
#' M_32S <- mean(`32S`)
#'
#' # mean 34S
#' M_34S <- mean(`34S`)
#'
#' # mean R
#' M_R <- M_34S / M_32S
#'
#' # standard deviation 32S
#' S_32S <- sd(`32S`)
#'
#' # standard deviation 34S
#' S_34S <- sd(`34S`)
#'
#' # propagation
#' des_SD_prop(M_R, `34S`, `32S`, M_34S, M_32S, S_34S, S_32S, 10, "rse")
#'
#' @export
des_SD_prop <- function(M_R, ion1, ion2, M_ion1, M_ion2, S_ion1, S_ion2, n,
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

#' @rdname  des_SD_prop
#'
#' @export
hat_SD_prop <- function(M_R, ion1, ion2 ,N_ion1, N_ion2, n, type = "sd"){

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



sum_sd_prop <- function(x, type = "rse"){

  n <- length(x)
  sum_sd <- sqrt( 1 / sum(x^-2))

}

#-------------------------------------------------------------------------------
# Function for parsing, testing and validation (NOT EXPORTET)
#-------------------------------------------------------------------------------
# Function to detect and remove zero counts of desired isotope system
zeroCt <- function(df, N, species, ion1, ion2, ...){

  N <- enquo(N)
  species <- enquo(species)
  gr_by <- enquos(...)

  gr.ls <- df %>%
    filter(!!species == ion1 | !!species == ion2)

  if (any(gr.ls %>% select(!!N) %>% pull(!!N) == 0)) {

    warning("zero counts present and removed")

  }

  gr.ls <- gr.ls %>%
    filter(!!N == 0) %>%
    select(!!!gr_by)

  df <- anti_join(df, gr.ls, by = sapply(gr_by, as_name))

}

# Function for covariate convertion of isotope systems
cov_R <- function(df, species, ion1, ion2, ...){

  gr_by <- enquos(...)
  species <- enquo(species)

  df <- df %>%
    group_by(!!! gr_by, !! species) %>%
    mutate(ID = row_number()) %>%
    ungroup() %>%
# Uniqually identifies ion pairs for calculating isotope ratios (in case the
# data file does not contain an ID)
    tidyr::unite(col = ID, !!! gr_by, ID, sep = "/", remove = FALSE)

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

# Function to parse chemical species names in latex
latex_parser <- function(species, ion1, ion2){

  if(length(species) != 0){

    # seperator for poly atomic speices
    if (any(str_detect(species, "[:blank:]"))){
      sep <- "--"
    } else {
      sep <- ""}

    species <- paste(
      paste0(
        paste0("$\\phantom{,}^{",
               str_extract_all(species, "[:digit:]+(?=[:alpha:])")[[1]][1],
               "}$"),
        str_extract_all(species, "[:alpha:]+")[[1]][1],
        paste0("$_{",
               str_extract_all(species, "(?<=[:alpha:])[:digit:]+")[[1]][1],
               "}$")
      ),
      paste0(
        paste0("$\\phantom{,}^{",
               str_extract_all(species, "[:digit:]+(?=[:alpha:])")[[1]][2],
               "}$"),
        str_extract_all(species, "[:alpha:]+")[[1]][2],
        paste0("$_{",
               str_extract_all(species, "(?<=[:alpha:])[:digit:]+")[[1]][2],
               "}$")
      ),
      sep = sep
    )

    species <- str_replace_all(species, "NA", "")

  }else{

    R <- paste(
      paste0(
        paste0("$\\phantom{,}^{",str_extract(ion1, "^\\d+"),"}$"),
        str_extract(ion1, "\\D+"),
        paste0("$_{",str_extract(ion1, "\\d*$"),"}$")
      ),
      paste0(
        paste0("$\\phantom{,}^{",str_extract(ion2, "^\\d+"),"}$"),
        str_extract(ion2, "\\D+"),
        paste0("$_{",str_extract(ion2, "\\d*$"),"}$")
      ), sep = "/")
    }
  }
