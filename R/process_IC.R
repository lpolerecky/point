#' Process raw ion count data
#'
#' \code{cor_IC} function for processing ion count data.
#'
#' The accuracy of pulsed ion counting is influenced by systematic errors which
#' depend on the ion counting system. Deadtime and EM yield are two
#' prominent effects for the electron multiplier systems. The deadtime refer to
#' the timewindow when the system does not register counts; this occurs when
#' incident ions strike the EM in a small enough time window in which the EM
#' channel is electronically paralysed. The EM yield is the ratio between the
#' number of output pulses counted after the EM  discriminator threshold and
#' the number of ions arriving at the EM. The latter can be gauged with the peak
#' height distribution (PHD) which is the probability for an EM output to have a
#' certain voltage amplitude.
#'
#' @param df A tibble containing raw ion count data.
#' @param N A variable constituting the ion counts.
#' @param t A variable constituting the time increments.
#' @param Det A character string or variable for the detection system ("EM" or
#' "FC").
#' @param deadtime A numeric value for the deadtime of the EM system.
#' @param thr_PHD A numeric value for the disrcriminator threshold of the EM.
#' system
#'
#' @return A \code{\link[tibble:tibble]{tibble}} containing the original dataset
#' and adds the variables: \code{Xt.rw}, ion count rates uncorrected for
#' detection device-specific biases; \code{Xt.pr}, ion count rates corrected for
#' detection device-specific biases; and \code{N.pr}, counts corrected for
#' detection device-specific biases.
#' @export
#' @examples
#' # Use point_example() to access the examples bundled with this package
#'
#' # raw data containing 13C and 12C counts on carbonate
#' tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # Processing raw ion count data
#' tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt, deadtime = 44, thr_PHD = 50)
cor_IC <-function(df, N, t, Det, deadtime = 0, thr_PHD = 0){

  stopifnot(tibble::is_tibble(df))
  stopifnot(is.numeric(deadtime))
  stopifnot(is.numeric(thr_PHD))

  N <- enquo(N)
  t <- enquo(t)
  Det <- enquo(Det)

# The corrections
  args <- lst(
              # Time increments (time between measurements minus blanking time)
              quo(min(!! t) - .data$tc.mt),
              # Count rates
              quo(!! N / .data$diff.t),
              # Deadtime correction on count rates
              quo(if_else(!! Det == "EM",
                          cor_DT(.data$Xt.rw, deadtime),
                          .data$Xt.rw)),
              # Deadtime correction on counts
              quo(if_else(!! Det == "EM",
                          cor_DT(.data$Xt.rw, deadtime) * .data$diff.t,
                          !! N)),
              # Yield correction on count rates
              quo(if_else(!! Det == "EM",
                          cor_yield(.data$Xt.pr,
                                    .data$mean_PHD,
                                    .data$SD_PHD,
                                    thr_PHD),
                          .data$Xt.pr)),
              # Yield correction on counts
              quo(if_else(!! Det == "EM",
                          cor_yield(.data$Xt.pr,
                                    .data$mean_PHD,
                                    .data$SD_PHD,
                                    thr_PHD) * .data$diff.t,
                          !! quo_updt2(N, "pr")))
              )

# The correction names (depend on user-supplied expression)
  ls.names <- c("diff.t", "Xt.rw", "Xt.pr",
                as_name(quo_updt2(N, "pr")),
                "Xt.pr",
                as_name(quo_updt2(N, "pr")))

# Set correction names
  args <- set_names(args, nm = ls.names)

# Execute corrections
  tb.pr <- df %>%
             mutate(!!! args)

  return(tb.pr %>% select(-.data$diff.t))

  }

#' Correct ion detection bias
#'
#' \code{cor_yield} function to correct counting bias related to EM Yield.
#' \code{cor_DT} function to correct counting bias related to deadtime
#'
#' The accuracy of pulsed ion counting is influenced by systematic errors which
#' depend on the ion counting system. Deadtime and EM yield are two
#' prominent effects for the electron multiplier systems. The deadtime refer to
#' the timewindow when the system does not register counts; this occurs when
#' incident ions strike the EM in a small enough time window in which the EM
#' channel is electronically paralysed. The EM yield is the ratio between the
#' number of output pulses counted after the EM  discriminator threshold and
#' the number of ions arriving at the EM. The latter can be gauged with the peak
#' height distribution (PHD) which is the probability for an EM output to have a
#' certain voltage amplitude.
#'
#' @param Xt A numeric vector containing raw ion count data.
#' @param mean_PHD A numeric vector coontaining the mean PHD value.
#' @param SD_PHD A numeric vector coontaining the standard deviation of the
#' PHD.
#' @param thr_PHD A numeric value for the disrcriminator threshold of the EM
#' system.
#' @param deadtime A numeric value for the deadtime of the EM system.
#'
#' @return A numeric vector with the corrected count rates.
#' @export
#' @examples
#' # Original count rate of a chemical species
#' x <- 30000
#'
#' # Corrected count rate for EM Yield with a threshold of 50 V
#' cor_yield(x, mean_PHD = 210, SD_PHD = 60, thr_PHD = 50)
#'
#' # Corrected count rate for a deadtime of 44 ns
#' cor_DT(x, 44)
cor_yield <- function(Xt, mean_PHD, SD_PHD, thr_PHD){

# Stop execution if threshold = 0 an dterurn Xt
  if (thr_PHD == 0){ return(Xt) } else {
# Lambda parameter
  lambda <- (2 * mean_PHD^2) / (SD_PHD^2 + mean_PHD)

# Probability parameter
  prob <- (SD_PHD^2 - mean_PHD) / (SD_PHD^2 + mean_PHD)
  prob <- if_else(prob < 0 | prob >= 1, NA_real_, prob, NA_real_)

  l.params <- lst(a = lambda,
                  b = prob,
                  c = rep(thr_PHD, length(.data$b)),
                  d = Xt)

  f.polya <- function(a, b, c, d){

    if (!(is.na(a) | is.na(b))){

      Y <- polyaAeppli::pPolyaAeppli(c,
                                     lambda = a,
                                     prob = b,
                                     lower.tail = FALSE)
      Xt.pr <- d / Y

      return(Xt.pr)

    }else{

      Xt.pr <- d

      return(Xt.pr)

    }
  }
  }

  Xt.pr <- purrr::pmap_dbl(l.params, f.polya)

  return(Xt.pr)

  }

#' @rdname  cor_yield
#'
#' @export
cor_DT <- function(Xt, deadtime) {

  # Stop execution if threshold = 0 an dterurn Xt
  if (deadtime == 0){ return(Xt) } else {

  Xt.pr <- Xt / (1 - (Xt * deadtime * 10^-9))

  return(Xt.pr)
  }
}


# Function which updates quosures for subsequent tidy evaluation
quo_updt2 <- function(my_q, txt = NULL){

  # Get expressions
  old_expr <- get_expr(my_q)

  if(str_detect(old_expr, "[:punct:]")){

    new_chr <- str_split(old_expr, "[:punct:]")[[1]][1]

  }

  # Modify expression, turn expr into a character string
  new_chr <- paste(new_chr, txt, sep = ".")

  # New expression from character (remove whitespace)
  new_expr <- parse_expr(new_chr)

  # Update old quosure
  set_expr(my_q, new_expr)

}

#' Predicting trends in ionization efficiency
#' @export
#' @examples
#' # raw ion counts on carbonate
#' tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))
#'
#' # Processing raw ion count data
#' tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt)
#'
#' # remove zero count analysis
#' tb.zr <- zeroCt(
#'   tb.pr,
#'   N.pr,
#'   species.nm,
#'   "12C",
#'   "40Ca 16O",
#'   file.nm,
#'   warn = FALSE
#'   )
#'
#' # predict ionization trends
#' tb.ion<- predict_ionize(
#'   df = tb.zr,
#'   Xt = Xt.pr,
#'   N = N.pr,
#'   t = t.rw,
#'   species = species.nm,
#'   file.nm
#'   )
predict_ionize <- function(df, Xt, N, t, species, ..., nest = FALSE, group = NULL, .plot = TRUE){

  # species <- enquo(species)
  gr_by <- enquos(..., species)
  group <- enquo(group)
  Xt <- enquo(Xt)
  N <- enquo(N)
  t <- enquo(t)

  # new quosures
  Xt.model <- quo_updt(Xt, "model")
  Xt.l0 <- quo_updt(Xt, "l0")
  N.l0 <-  quo_updt(N, "l0")


  df <- df %>%
    # group_by( %>%
    nest(data = -c(!!! gr_by)) %>%
    mutate(!! Xt.model := purrr::map(data, ~gam_fun(.x, Xt, t))) %>%
    unnest(cols = c(data, !! Xt.model)) %>%
# De-trended timeseries
    mutate(
      "{{Xt}}.l0" := mean(!! Xt.model) + (!! Xt - !! Xt.model),
      "{{N}}.l0" := !! Xt.l0 * (min(!! t) - .data$tc.mt),
      hat_Y = !! Xt.model,
      hat_min = !! Xt.model,
      hat_max = !! Xt.model,
      flag = "good"
           ) %>%
    ungroup()

  if (nest) {
    df<- df %>%
      nest(data = -c(!!! gr_by[!sapply(gr_by, as_name) %in% as_name(group)])) %>%
      mutate(mlm = purrr::map(data, ~gam_fun(.x, Xt, t, group = group, output = "model")))
      # unnest()
    return(df)
  }

  if (.plot) {

    facets_gr <- quos(!!! gr_by)#quos(!!! gr_by, !! species)

    plot_args <- list2(
      df = df,
      stat = NULL,
      y = Xt,
      x = t,
      ion1 = NA,
      ion2 = NA,
      plot_title = "timeseries",
      plot_type = "static",
      args = NULL,
      iso = FALSE,
      isoscale = NULL,
      !!! facets_gr
      )

    # control with environment
    data_env <- env(data = df)
    print(eval(expr(gg_default(!!! plot_args)), data_env))

    return(df %>% select(-c(hat_Y, hat_min, hat_max, flag)))

  }
}

#' Predicting variation in ion based on ionization trend
#' @export
predict_var <- function(df, Xt, N, t, species, ion, ...){

  Xt <- enquo(Xt)
  N <- enquo(N)
  t <- enquo(t)
  species <- enquo(species)
  gr_by <- enquos(...)

  df.ion <- filter(df, !!species == ion)
  df.l0 <- predict_ionize(df.ion, !!Xt, !!N, !!t, !!species, !!!gr_by)

  Xt.l0 <- quo_updt(Xt, "l0")
  N.l0 <- quo_updt(N, "l0")

  tb.Xt <- stat_Xt(df.ion, !!Xt, !!N, !!species, !!!gr_by)
  tb.Xt.l0 <- stat_Xt(df.l0, !!Xt.l0, !!N.l0, !!species, !!!gr_by)

  RS_Xt <- transmute(
    tb.Xt,
    !!!gr_by,
    !!quo_updt(Xt, ion, x = "RS") :=
    !!quo_updt(Xt, x = "RS") -
    pull(tb.Xt.l0, !!quo_updt(Xt.l0, x = "RS"))
    )




  }

gam_fun <- function(data, Xt, t, group = NULL, output = "predict"){


  mthd<- "norm"
  s_call <- call2("s", get_expr(t))
  form_gam <- new_formula(get_expr(Xt), s_call, env = caller_env())
  if (!is.null(group)) {
    form_ran <- list2(!!group := new_formula(NULL, get_expr(t)))
    mthd <- "mix"
  }


  gam_switch <- function(type) {

    switch(
      type,
      norm = eval(call2("gam", form_gam, method = "REML", data = expr(data))),
      mix = eval(call2("gamm", form_gam, random = form_ran, method = "REML", data = expr(data)))
           )
  }


  model <- gam_switch(type = mthd)

  if(output == "predict") return(unname(fitted(model)))
  if(output == "model") return(model)
}

#' QSA test
#'
#' \code{QSA_test} function for processing ion count data.
#'
#' The accuracy of pulsed ion counting is influenced by systematic errors which
#' depend on the ion counting system. Deadtime and EM yield are two
#' prominent effects for the electron multiplier systems. The deadtime refer to
#' the timewindow when the system does not register counts; this occurs when
#' incident ions strike the EM in a small enough time window in which the EM
#' channel is electronically paralysed. The EM yield is the ratio between the
#' number of output pulses counted after the EM  discriminator threshold and
#' the number of ions arriving at the EM. The latter can be gauged with the peak
#' height distribution (PHD) which is the probability for an EM output to have a
#' certain voltage amplitude.
#'
#' @param df A tibble containing raw ion count data.
#' @param N A variable constituting the ion counts.
#' @param t A variable constituting the time increments.
#' @param Det A character string or variable for the detection system ("EM" or
#' "FC").
#' @param deadtime A numeric value for the deadtime of the EM system.
#' @param thr_PHD A numeric value for the disrcriminator threshold of the EM.
#' system
#'
#' @return A \code{\link[tibble:tibble]{tibble}} containing the original dataset
#' and adds the variables: \code{Xt.rw}, ion count rates uncorrected for
#' detection device-specific biases; \code{Xt.pr}, ion count rates corrected for
#' detection device-specific biases; and \code{N.pr}, counts corrected for
#' detection device-specific biases.
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
#' # QSA test
#' tb.QSA <- QSA_test(tb.pr, Xt.pr, N.pr, species.nm, "13C", "12C", file.nm)
QSA_test <- function(.df, .Xt, .N, .species, .ion1, .ion2, ..., .plot = TRUE){

  # # Quoting the call (user-supplied expressions)
  Xt <- enquo(.Xt)
  gr_by <- enquos(...)
  Xt1 <- quo_updt(Xt, .ion1) # count rate
  Xt2 <- quo_updt(Xt, .ion2) # count rate
  R.Xt <- quo_updt(Xt, x = "R")

  df <- .df %>%
    zeroCt(N = {{ .N }},
           species = {{ .species }},
           ion1 = {{ .ion1 }},
           ion2 = {{ .ion2 }},
           !!! gr_by,
           warn = FALSE
           ) %>%
    cov_R(species = {{ .species }},
          ion1 = {{ .ion1 }},
          ion2 = {{ .ion2 }},
          !!! gr_by,
          preserve = TRUE
          ) %>%
    mutate("R_{{.Xt}}" := {{ Xt1 }} / {{ Xt2 }}) %>%
    group_by(!!! gr_by) %>%
    tidyr::nest()  %>%
    mutate(
      QSA_lm = purrr::map(data, ~lm_form(.x, R.Xt, Xt2)),
      B1 = purrr::map(QSA_lm, ~{pull(broom::tidy(.x), estimate)[2]}),
      t.val = purrr::map(QSA_lm, ~{pull(broom::tidy(.x), statistic)[2]}),
      p.val = purrr::map(QSA_lm, ~{pull(broom::tidy(.x), p.value)[2]})
          ) %>%
    tidyr::unnest(cols = c(data, B1, t.val, p.val)) %>%
    mutate(
      # Chi_R2 = R2 * n(),
      # Chi_crit = qchisq(.95, df = 1),
      # flag_QSA = if_else(Chi_R2 > Chi_crit & B1 < 0,
      #                    "Ha QSA",
      #                    "H0 no QSA"
      #                    ),
      flag = if_else(p.val < 0.01,
                     "bad",
                     "good"
                     )
           )

  if (.plot) {

  stat_lab <- stat_select2(df, gr_by)

  plot_args <- list2(df = df,
                     stat = stat_lab,
                     y = R.Xt,
                     x = Xt2,
                     ion1 = .ion1,
                     ion2 = .ion2,
                     plot_title = "QSA",
                     plot_type = "static",
                     args = NULL,
                     iso = FALSE,
                     isoscale = NULL,
                     !!! gr_by
                     )


  # control with environment
  data_env <- env(data = df)
  eval(expr(gg_default(!!! plot_args)), data_env) +
    geom_smooth(method = "lm", formula = 'y ~ x')
  }

}

# add this lateron
#gg_default(df.new, stat_lab = NULL, y = R_Xt.pr, x = Xt.pr.12C, "13C", "12C", plot_title = "QSA", plot_type = "static", iso =FALSE, isoscale = NULL) + facet_wrap(vars(file.nm), scale = "free")

# linear model call
#' @export
lm_form <- function(data, arg1, arg2, flag = NULL, trans = NULL, vorce = NULL, nest = NULL, type = "OLS") {

  # if (type == "Rm" | type == "LME"){
  #   if(is.null(flag)) {
  #     if(trans){
  #       call_lm <- new_formula(quo_get_expr(arg1), parse_expr(paste0("-1 + I(", as_name(arg2), "/ 1000)")))
  #     } else {
  #   call_lm <- new_formula(quo_get_expr(arg1), parse_expr(paste0("-1 +", as_name(arg2))))
  #     }
  #   }
  #   call_lm <- new_formula(quo_get_expr(arg1), parse_expr(paste0("-1 +", as_name(arg2), "*", as_name(flag))))
  # }

  # if (type == "Rm" | type == "LME" | type == "GLS"){
  #   if (is.null(flag)) {
  #     if(trans){
  #       if (vorce == "intra") {
  #         call_lm <- new_formula(quo_get_expr(arg1), parse_expr(paste0("-1 + I(", as_name(arg2), "/ 1000) + I(", as_name(arg2), "/ 1000) : execution")))
  #         } else {
  #       call_lm <- new_formula(quo_get_expr(arg1), parse_expr(paste0("-1 + I(", as_name(arg2), "/ 1000)")))
  #         }
  #     } else {
  #       call_lm <- new_formula(quo_get_expr(arg1), parse_expr(paste0("-1 +", as_name(arg2))))
  #     }
  #   } else {
  #     call_lm <- new_formula(quo_get_expr(arg1), parse_expr(paste0("-1 +", as_name(arg2), "*", as_name(flag))))
  #   }
  # }

  lm_switch <- function(type, itrc = NULL, trns = NULL){

    if (type %in% c("Rm", "LME", "GLS")) type <- "Rm" else type = "OLS"
    if (is_quosure(itrc)) {
      type <- paste(type, deparse(substitute(itrc)), sep = "_")
      } else {
        if (is.null(trns)) {
          type <- paste(type, itrc, sep = "_")
          } else {
            type <- paste(type, trns, sep = "_")
          }
        }

    switch(type,
           Rm_flag = new_formula(quo_get_expr(arg1), parse_expr(paste0("-1 +", as_name(arg2), "*", as_name(flag)))),
           Rm_ppt = new_formula(quo_get_expr(arg1), parse_expr(paste0("-1 + I(", as_name(arg2), "/ 1000)"))),
           Rm_ = new_formula(quo_get_expr(arg1), parse_expr(paste0("-1 +", as_name(arg2)))),
           Rm_log = new_formula(parse_expr(paste0("log(", as_name(arg1), ")")), parse_expr(paste0("-1 + log(", as_name(arg2), ")"))),
           Rm_intra =  new_formula(quo_get_expr(arg1), parse_expr(paste0("-1 + I(", as_name(arg2), "/ 1000) + I(", as_name(arg2), "/ 1000) : execution"))),
           OLS_ = new_formula(quo_get_expr(arg1), quo_get_expr(arg2))
           )

  }

  ran_switch <- function(clst, trns  = NULL){

    if (is_character(trns)) clst <- paste(clst, trns, sep = "_")
    switch(clst,
           intra = new_formula(NULL, parse_expr(paste("-1 + execution", as_name(nest), sep = "|"))),
           inter_ppt = new_formula(NULL, parse_expr(paste(paste0("-1 + I(", as_name(arg2), "/ 1000)"), paste("execution", as_name(nest), sep = "/"), sep = "|"))),
           inter_log = new_formula(NULL, parse_expr(paste(paste0("-1 + log(", as_name(arg2), ")"), paste("execution", as_name(nest), sep = "/"), sep = "|")))
           )

  }

  # if (type == "OLS") call_lm <- new_formula(quo_get_expr(arg1), quo_get_expr(arg2))

#   # if (type == "LME" & is.null(flag)) call_ran <- new_formula(NULL, parse_expr(paste(paste0(as_name(arg2), "-1"), as_name(nest), sep = "|")))
#   if (type == "LME" & vorce == "intra") {
#       #call_ran <- new_formula(NULL, parse_expr(paste(paste0(as_name(arg2), "-1"), "execution", sep = "|")))
# # autoregressive situation of lonitudinal executions
#       call_ran <- new_formula(NULL, parse_expr(paste("-1 + execution", as_name(nest), sep = "|")))
#       # cor_struc <- call2("corAR1", form = call_ran,
#       #                    .ns = "nlme")
#       cor_struc <- NULL
#       mthd <-  "ML"
#       # call_ran <- new_formula(NULL, parse_expr(paste(paste0("-1 + I(", as_name(arg2), "/ 1000)"), paste(as_name(nest), "execution", sep = "/"), sep = "|")))
#       #f_rhs(call_lm) <- parse_expr(paste(paste0(f_name(call_lm), "+ (-1 + log(", as_name(arg2), ")"), paste(as_name(nest), "execution)", sep = ":"), sep = "|"))
#       }
#   if (type == "LME" & vorce == "inter") {
#       call_ran <- new_formula(NULL, parse_expr(paste(paste0("-1 + I(", as_name(arg2), "/ 1000)"), paste("execution", as_name(nest), sep = "/"), sep = "|")))
#       cor_struc <- NULL
#       mthd <-  "REML"
#       #f_rhs(call_lm) <- parse_expr(paste(paste0(f_name(call_lm),"+ (-1 + log(", as_name(arg2), ")"), paste("execution", paste0(as_name(nest), ":"), sep = "/"), sep = "|"))
#       #call_ran <- list2(execution = new_formula(NULL, 1), !! nest := new_formula(NULL, parse_expr(as_name(arg2))))
#       }

  wght_switch <- function(type, trns = NULL){
    if (type %in% c("LME", "GLS")) type <- "LME" else type = "Rm"
    if (!is.null(trns) && trns == "log") type <- paste(type, trns, sep = "_")

    switch(type,
           Rm = 1 / (pull(data, !! arg2)) ,
           LME_log = new_formula(NULL, parse_expr(paste0("I(", paste(1, paste0("log(",as_name(arg2),")"), sep = "/"), ")"))),
           LME = new_formula(NULL, parse_expr(paste0("I(", paste(1, paste0("I(",as_name(arg2),"/ 1000)"), sep = "/"), ")")))
           )

  }

  # switch between types of linear model
    lm_method <- function(type) {
      switch(type,
             OLS = eval(call2("lm", lm_switch("OLS") , data = expr(data))),
             GLS = eval(call2("gls",lm_switch("GLS", trans), data =expr(data), weights = wght_switch("GLS", trans), .ns = "nlme")),
             Rm = eval(call2("lm", lm_switch("Rm", flag), data = expr(data), weights = expr(wght_switch("Rm")))),
             LME = eval(
              call2(
               "lme",
               lm_switch("LME", vorce, trans),
               random = ran_switch(vorce, trans),
               data = expr(data),
               method = if_else(vorce == "inter", "REML", "ML"),
               weights = wght_switch("LME", trans),
               .ns = "nlme"
                   )
              )
             )
    }
    # if(trans){
    #   if (type == "Rm"){
    #   wght <- (pull(data, !! arg2)) /1000
    #   } else {
    #     call_wght <- new_formula(NULL, parse_expr(paste0("I(", paste(1, paste0("I(",as_name(arg2),"/ 1000)"), sep = "/"), ")")))
    #   }
    # } else {
    #   wght <- pull(data, !! arg2)
    # }

  #call_wght <- new_formula(NULL, parse_expr(paste0("I(", paste(1, as_name(arg2), sep = "/"), ")")))

# # switch between OLS and ratio method type linear model
#   lm_method <- function(type) {
#     switch(type,
#            OLS = eval(call2("lm", call_lm, data = expr(data))),
#            GLS = eval(call2("gls", call_lm, data =expr(data), weights = call_wght, .ns = "nlme")),
#            Rm = eval(call2("lm", call_lm, data = expr(data), weights = expr(1 / wght))),
#            LME = eval(call2("lme",
#                             call_lm,
#                             random = call_ran,
#                             data = expr(data),
#                             method = mthd,
#                             correlation = cor_struc,
#                             weights = call_wght,
#                             .ns = "nlme"
#                             )
#                       )
#
#            # LME = eval(call2("lmer", call_lm, data = expr(data), weights = expr(1 / wght), .ns = "lme4"))
#            )
#   }

lm_method(type)


}

stat_lab2 <- function(a, b, c){

    expr_lm <- lst(substitute(beta[1] == a * " \n " ~
                              t[beta[1]] == b * " \n " ~
                              "(" * p == c * ";" ~
                              "H0:"~ beta[1] == 0 * ")"
                              ,
                              lst(a = sprintf("%+.1e", a),
                                  b = sprintf("%+.1e", b),
                                  c = sprintf("%.3f", c)
                                  )
                              )
                   )

    expr_lm %>%
      do.call("expression", .) %>%
      factor(x = as.character(a),
             labels = .
             )
}

# Stat labels selection
stat_select2 <- function(df, facets_gr) {

  df %>%
    distinct(!!! facets_gr, .keep_all = TRUE) %>%
    transmute(trans = "original",
              lb = purrr::pmap(lst(
                a = B1,
                b = t.val,
                c = p.val
                ),
                stat_lab2),
              vjust = 3.5,
              !!! facets_gr
              ) %>%
    tidyr::unnest(cols = lb)

}
