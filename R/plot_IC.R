# Plotting diagnostics
gg_IC <- function(.df, .ion1, .ion2, .method, .plot_type, ..., .Xt = Xt.pr,
                  .N = N.pr, .species = species.nm, .t = t.nm, .flag = flag,
                  .labels = NULL, .rep = 1, .alpha_level, .plot_iso = FALSE){

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)

  # stat_R variables
  args <- enquos(.Xt = .Xt, .N = .N, .species = .species, .t = .t)

  flag <-  enquo(.flag)
  # Heavy isotope
  Xt1 <- quo_updt(args[[".Xt"]], post = .ion1)
  # Light isotope
  Xt2 <- quo_updt(args[[".Xt"]], post = .ion2)
  # Fitted heavy isotope and variance (sigma)
  hat_Xt1 <- quo_updt(args[[".Xt"]], pre = "hat", post = .ion1)
  hat_s_Xt1 <- quo_updt(args[[".Xt"]], pre = "hat_s", post = .ion1)

  # Filter execution or downsample for animation
  if (.plot_type == "static") .df <- filter(.df, execution == .rep)
  # if (.plot_type == "anim") {
  #   df_sl <- filter(.df, execution == 1)
  #   df_sl <- slice_sample(group_by(df_sl, !!!gr_by, !!flag), prop = 0.05) %>%
  #     ungroup() %>%
  #     select(!!!gr_by, !!args[[".t"]])
  #   .df <- semi_join(.df, df_sl, by = sapply(c(gr_by, args[[".t"]]), as_name))
  #
  #   }
  if (.plot_type == "anim") .df <- slice_sample(group_by(.df, !!!gr_by, !!flag),
                                                prop = 0.05)

  plot_args <- list2(.df = .df, .x = Xt2, .y = Xt1, .flag = flag,
                     .diag_type = .method, .plot_type = .plot_type, !!! gr_by,
                     .labels = .labels, .hat = hat_Xt1, .sd = hat_s_Xt1,
                     .alpha_level = .alpha_level)

  # Residual leverage plot
  if(.method == "norm_E") {
      plot_args[[".y"]] <- quo(studE)
      plot_args[[".x"]] <- quo(hat_Xi)
      plot_args[[".hat"]] <- NULL
      plot_args[[".sd"]] <- NULL
    }

  # Normal QQ plot
  if (.method == "QQ"){
      plot_args[[".y"]] <- quo(RQ)
      plot_args[[".x"]] <- quo(TQ)
      plot_args[[".hat"]] <- quo_updt(plot_args[[".y"]], pre = "hat")
      plot_args[[".sd"]] <- NULL
      plot_args[[".se"]] <- quo_updt(plot_args[[".y"]], pre = "hat_e")
    }

  # Scale location plot
  if (.method  == "CV"){
      plot_args[[".y"]] <- quo(studE)
      plot_args[[".x"]] <- hat_Xt1
      plot_args[[".hat"]] <- 0
      plot_args[[".sd"]] <- NULL
      plot_args[[".cv"]] <- 3.5
  }

  # Execute
  data_env <- env(data = .df)
  p <- eval(expr(gg_base(!!!plot_args)), data_env)

  if (.plot_type == "anim")  {
    p <- plotly::ggplotly(p, tooltip = "text") %>%
      plotly::animation_opts(frame = 150, transition = 0)
    return(p)
    }
  if (.plot_type == "static") return(p)

  }


gg_base <- function(.df, .x, .y, .flag, .diag_type, .plot_type, ...,
                    .labels = NULL, .geom = "point", .hat = NULL, .sd = NULL,
                    .se = NULL, .cv = NULL, .alpha_level, .rug = FALSE){

  gr_by <- enquos(...)
  if (str_detect(as_name(.x), "[[:digit:]]") | .diag_type == "CV") {
    ion2 <- str_split(as_name(.x), "[[:punct:]]")[[1]] %>% tail(1)
    } else {
      ion2 <- NULL
    }
  if (str_detect(as_name(.y), "[[:digit:]]")) {
    ion1 <- str_split(as_name(.y), "[[:punct:]]")[[1]] %>% tail(1)
    } else {
      ion1 <- NULL
    }

  # Adjust point alpha
  alpha_sc <- median(count(.df, !!!gr_by)$n) /
    ifelse(.plot_type == "static", 1e5, 1e3)
  # Filter correct titles
  ttl <- filter(names_diag, name == .diag_type)

  # Uncertainty bounds
  if (!is.null(.sd) | !is.null(.se) | !is.null(.cv)) {
    if (!is.null(.sd)){
      .df <- mutate(
        .df,
        fct  = qnorm((1 - .alpha_level / 2)),
        lower = !!.hat - .data$fct * !!.sd,
        upper = !!.hat + .data$fct * !!.sd
        )
    }
    if (!is.null(.se)) {
      .df <- mutate(
        .df,
        fct  = qt((1 - .alpha_level / 2), n() - 1),
        lower = !!.hat - .data$fct * !!.se,
        upper = !!.hat + .data$fct * !!.se
        )
    }
    if (!is.null(.cv)) {
      .df <- mutate(
        .df,
        lower = !!.hat - !!.cv,
        upper = !!.hat + !!.cv
      )
    }
  .lower <- quo(lower)
  .upper <- quo(upper)
  }

  if (!is.null(.labels)) {
    tb_labs <- distinct(.df, !!!gr_by, .keep_all = TRUE) %>%
      select(!!!gr_by, starts_with(paste0(point::names_stat_R$name,"_"))) %>%
      pivot_longer(
        -c(!!!gr_by),
        names_to = c("stat", ".value"),
        names_sep = "\\_R\\_"
        )
    lbs <- purrr::map2(
      tb_labs$stat,
      pull(select(tb_labs, last_col()), 1),
      ~stat_labeller("R", stat = .x, value = .y, "expr")
      )
    tb_labs <- tibble::add_column(tb_labs, lbs = lbs)
    }

  if (.plot_type == "anim") {
    p <- ggplot(
      data = .df,
      aes(
        x = !!.x,
        y = !!.y,
        color = {{.flag}},
        frame = execution,
        text = eval(text_labs(.diag_type, x = !!.x, y = !!.y, ion1, ion2))
        )
      )
  }
  if (.plot_type == "static") {
    p <- ggplot(data = .df, aes(x = !!.x, y = !!.y, color = {{.flag}}))
    }
  p <- p + facet_wrap(vars(!!!gr_by), scales = "free")
  if (!any(sapply(list(.hat, .sd, .se, .cv), is.null)) &
      .plot_type != "anim") {
    p <- p + geom_ribbon(
      aes(ymin = !!.lower, ymax = !!.upper),
      color = "black",
      fill = "aliceblue",
      alpha = 0.6,
      linetype = 3,
      size = 0.5
      ) +
      geom_line(aes(y = !!.hat), color = "black", linetype = 2, size = 0.5)
    }
  if (.geom == "point") p <- p + geom_point(alpha = alpha_sc)
  if (.geom == "hexbin") p <- p + geom_hexbin()
  if (.geom == "dens2d") {
    p <- p + stat_density_2d(
      aes(fill = ..ndensity..),
      geom = "raster",
      contour = FALSE,
      contour_var = "count",
      show.legend = FALSE
      ) +
      scale_fill_distiller(
        limits = c(0.01, 1),
        breaks = seq(0.01, 1, length.out = 100),
        palette = "YlOrRd",
        direction = 1,
        na.value = "transparent"
        )
    }
  if (.rug) p <- p + geom_rug(sides = "tr", alpha = 0.01)
  if (!is.null(.labels)) {
    p <- p + ggrepel::geom_text_repel(
      data = filter(tb_labs, stat %in% .labels),
      aes(
        x = -Inf,
        y = Inf,
        label = lbs
        ),
      direction = "y",
      segment.color = "transparent",
      parse = TRUE,
      inherit.aes = FALSE
      )
  }
  p + scale_color_manual(values = c(ggplotColours(2)[2], ggplotColours(2)[1])) +
    scale_y_continuous(
      breaks = scales::extended_breaks(),
      labels = scales::label_scientific()
      ) +
    scale_x_continuous(
      breaks = scales::extended_breaks(),
      labels = scales::label_scientific()
      ) +
    labs(
      x = axis_labs(ttl$xaxis, ion2, .plot_type),
      y = axis_labs(ttl$yaxis, ion1, .plot_type),
      title = ttl$label
      ) +
    theme_classic() +
    theme(legend.position = "none")
  }


gg_IR <- function(.df, .lag, .acf, .flag, ..., .hat = NULL,
                  .sd = NULL){

  ggplot(.df, mapping = aes(x = {{.lag}}, y = {{.acf}}, color = {{.flag}})) +
    geom_hline(aes(yintercept = {{.hat}})) +
    geom_segment(mapping = aes(xend = {{.lag}}, yend = {{.hat}})) +
    geom_hline(aes(yintercept = -{{.sd}}), color = "darkblue") +
    geom_hline(aes(yintercept = {{.sd}}), color = "darkblue") +
    scale_color_manual(values = c(ggplotColours(2)[2], ggplotColours(2)[1])) +
    facet_wrap(vars(...), scales = "free") +
    ggtitle("ACF plot") +
    theme_classic() +
    theme(legend.position = "none")
}



axis_labs <- function(type, ion, plot_type){

  if (!(type == "ionct" | type ==  "hat_Y")) ion <- NULL

  if (plot_type == "static"){
  switch(
    type,
    ionct = substitute(a~"(ct/sec)", list(a = ion_labeller(ion, "expr"))),
    studE = expression("studentized residuals (" * italic(e)^"*" * ")"),
    TQ = "Theoretical quantiles",
    SQ = "Sample quantiles",
    Xct = "X (ct/sec)",
    time = "time (sec)",
    hat_Y = substitute("fitted value (" * hat(a) * ")",
                       list(a = ion_labeller(ion, "expr"))),
    hat_Xi = expression("hat-values" (italic(h))),
    R = "R"
    )
      } else {
        switch(
          type,
          ionct = str_replace("ion (ct/sec)", "ion", ion),
          studE = "studentized residuals",
          TQ = "Theoretical quantiles",
          SQ = "Sample quantiles",
          hat_Y = "fitted value ",
          hat_Xi = "hat-values"
          )
      }
}

text_labs <- function(diag_type, x, y, ion1, ion2, iso = FALSE) {

  # if (iso) {
  #   R_not <- paste0("δ", ion1, " (‰): ", delta)
  #   } else {
  #     R_not <- paste0("mean R: ", R)
  #   }
  x <- enexpr(x)
  y <- enexpr(y)
  ion1 <- enexpr(ion1)
  ion2 <- enexpr(ion2)
  if (diag_type == "CooksD" | diag_type == "Cameca") diag_type <- "Rm"
  switch(
      diag_type,
      Rm = expr(
        paste0(
          !!ion2, " (ct/sec): ", round(!!x, 0), '\n',
          !!ion1, " (ct/sec): ", round(!!x, 0), '\n',
          "R :", round(!!y / !!x, 4), '\n'
          )
        ),
      QQ = expr(
        paste0(
          "TQ: ", round(!!x, 2), '\n',
          "SQ: ", round(!!y, 2), '\n'
          )
        ),
      norm_E = expr(
        paste0(
          "e*: ", round(!!x, 4), '\n',
          "hat: ", round(!!y, 2), '\n'
          )
        ),
      CV = expr(
        paste0(
          "fitted value: ", round(!!x, 4), '\n',
          "e*: ", round(!!y, 2), '\n'
              )
        )
      )
}

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360 / n
  hcl(h = (seq(h[1], h[2], length = n)),
      c = 100,
      l = 65
      )
}

