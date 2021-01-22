
# linear model call
lm_form <- function(data, arg1, arg2, flag = NULL, trans = NULL, vorce = NULL, nest = NULL, type = "OLS") {



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



  wght_switch <- function(type, trns = NULL){
    if (type %in% c("LME", "GLS")) type <- "LME" else type = "Rm"
    if (!is.null(trns) && trns == "log") type <- paste(type, trns, sep = "_")

    switch(type,
           Rm = 1 / (pull(data, !! arg2)) ,
           LME_log = new_formula(NULL, parse_expr(paste0("I(", paste(1, paste0("log(",as_name(arg2),")"), sep = "/"), ")"))),
           LME = new_formula(NULL, parse_expr(paste0("I(", paste(1, paste0("I(",as_name(arg2),"/ 1000)"), sep = "/"), ")")))
    )

  }

  # wght <- pull(data, !! arg2)
  # switch between types of linear model
  lm_method <- function(type) {
    switch(type,
           OLS = eval(call2("lm", lm_switch("OLS") , data = expr(data))),
           GLS = eval(call2("gls",lm_switch("GLS", trans), data =expr(data), weights = wght_switch("GLS", trans), .ns = "nlme")),
           Rm = eval(call2("lm", lm_switch("Rm", flag), data = expr(data), weights = wght_switch("Rm"))),#expr(1 / wght))),
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

  lm_method(type)


}
