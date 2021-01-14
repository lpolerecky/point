
gam_fun <- function(df, Xt, t, group, Xt.mdl){

  # null model
  mlm_0 <- gam_form(df, Xt, t, output = "model")
  # full model
  mlm_1 <- gam_form(df, Xt, t, group, output = "model")

  # anova and extract components
  mlm_0_lme <- purrr::pluck(mlm_0, "lme")
  mlm_1_lme <- purrr::pluck(mlm_1, "lme")
  mlm_1_gam <-purrr::pluck(mlm_1, "gam")
  mlm_aov <- anova(mlm_1_lme, mlm_0_lme)
  ls_ran <- nlme::ranef(mlm_1_lme) %>%
    purrr::pluck(2)

  fit <- as.vector(unname(fitted(mlm_1_gam)))

  tb_ran <- tibble::as_tibble(ls_ran) %>%
    rename(ran_in.ml = contains("Intercept"), ran_sh.ml = !!t)

  sq_in <- rep(tb_ran$ran_in.ml, each = length(fit) / nrow(tb_ran))
  sq_sh <- rep(tb_ran$ran_sh.ml, each = length(fit) / nrow(tb_ran))

  tibble(
    ran_in.ml = sq_in,
    ran_sh.m = sq_sh,
    fix_in.ml = unname(nlme::fixef(mlm_1_lme))[1],
    fix_sh.ml = unname(nlme::fixef(mlm_1_lme))[2],
    p_intra.ml = purrr::discard(pull(mlm_aov, `p-value`), is.na),
    dAIC.ml = purrr::discard(diff(pull(mlm_aov, AIC)), is.na),
    !! Xt.mdl := fit
    )

}

gam_form <- function(data, Xt, t, group = NULL, output = "predict"){

    mthd <- "norm"
    s_call <- call2("s", get_expr(t))
    form_gam <- new_formula(get_expr(Xt), s_call, env = caller_env())

    if (!is.null(group)) {
      form_ran <- list2(!!group := new_formula(NULL, get_expr(t)))
      mthd <- "mix"
    }

    gam_switch <- function(type) {
      switch(
        type,
        norm = eval(
          call2(
            "gamm",
            form_gam,
            method = "REML",
            data = expr(data)
            )
          ),
        mix = eval(
          call2(
            "gamm",
            form_gam,
            random = form_ran,
            method = "REML",
            data = expr(data)
            )
          )
        )
    }

    model <- gam_switch(type = mthd)

    if(output == "predict") return(unname(fitted(model$gam)))
    if(output == "model") return(model)
}



