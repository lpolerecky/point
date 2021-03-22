
gam_fun <- function(IC, args, X.mdl){

  # null model
  mlm_0 <- gam_form(IC, args, output = "model")
  # full model
  mlm_1 <- gam_form(IC, args, args[[".nest"]], output = "model")

  # anova and extract components
  mlm_0_lme <- purrr::pluck(mlm_0, "lme")
  mlm_1_lme <- purrr::pluck(mlm_1, "lme")
  mlm_1_gam <-purrr::pluck(mlm_1, "gam")
  mlm_aov <- anova(mlm_1_lme, mlm_0_lme)
  ls_ran <- nlme::ranef(mlm_1_lme) %>%
    purrr::pluck(2)

  fit <- as.vector(unname(fitted(mlm_1_gam)))

  tb_ran <- tibble::as_tibble(ls_ran) %>%
    rename(ran_in.ml = contains("Intercept"), ran_sh.ml = !! args[[".t"]])

  sq_in <- rep(tb_ran$ran_in.ml, each = length(fit) / nrow(tb_ran))
  sq_sh <- rep(tb_ran$ran_sh.ml, each = length(fit) / nrow(tb_ran))

  tibble(
    ran_in.ml = sq_in,
    ran_sh.m = sq_sh,
    fix_in.ml = unname(nlme::fixef(mlm_1_lme))[1],
    fix_sh.ml = unname(nlme::fixef(mlm_1_lme))[2],
    p_intra.ml = purrr::discard(pull(mlm_aov, `p-value`), is.na),
    dAIC.ml = purrr::discard(diff(pull(mlm_aov, AIC)), is.na),
    !! X.mdl := fit
    )

}

gam_form <- function(data, args, nest = NULL, output = "predict"){

    mthd <- "norm"
    s_call <- call2("s", get_expr(args[[".t"]]))
    form_gam <- new_formula(
      get_expr(args[[".X"]]),
      s_call,
      env = env(data = data)
      )

    if (!is.null(nest)) {
      form_ran <- list2(
        !! nest :=
          new_formula(
            NULL,
            get_expr(args[[".t"]]),
            env = env(data = data)
            )
          )
      mthd <- "mix"
      }

    # Perform model in correct env
    data_env <- env(data = data)
    # Switch between fixed and mixed model
    gam_switch <- function(type) {
      switch(
        type,
        norm = eval(
          call2(
            "gamm",
            form_gam,
            method = "REML",
            data = expr(data)
            ),
          data_env
          ),
        mix = eval(
          call2(
            "gamm",
            form_gam,
            random = form_ran,
            method = "REML",
            data = expr(data)
            ),
          data_env
          )
        )
    }

    model <- gam_switch(type = mthd)

    if(output == "predict") return(unname(fitted(model$gam)))
    if(output == "model") return(model)
}



