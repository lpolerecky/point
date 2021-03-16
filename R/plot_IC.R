# Plotting diagnostics
gg_IC <- function(.IC, .ion1, .ion2, .method, .plot_type, ..., .Xt = Xt.pr,
                  .N = N.pr, .species = species.nm, .t = t.nm, .flag = flag,
                  .labels = NULL, .rep = 1, .alpha_level, .plot_iso = FALSE){

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)

  # stat_R variables
  args <- enquos(.Xt = .Xt, .N = .N, .species = .species, .t = .t)

  flag <-  enquo(.flag)
  # Common isotope
  N1 <- quo_updt(args[[".Xt"]], post = .ion1)
  Xt1 <- quo_updt(args[[".Xt"]], post = .ion1)
  # rare isotope
  Xt2 <- quo_updt(args[[".Xt"]], post = .ion2)
  # Fitted rare isotope and variance (sigma)
  hat_Xt1 <- quo_updt(args[[".Xt"]], pre = "hat", post = .ion1)
  hat_S_Xt1 <- quo_updt(args[[".N"]], pre = "hat_S", post = .ion1)

  # Filter execution
  if (.plot_type == "static") .IC <- filter(.IC, execution == .rep)

  # Plotting arguments
  plot_args <- list2(.IC = .IC, .x = Xt2, .y = Xt1, .flag = flag,
                     .diag_type = .method, .plot_type = .plot_type, !!! gr_by,
                     .labels = .labels, .hat = hat_Xt1, .sd = hat_S_Xt1 ,
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
  data_env <- env(data = .IC)
  p <- eval(expr(gg_base(!!!plot_args)), data_env)
  if (.plot_type == "static") return(p)
  }


gg_base <- function(.IC, .x, .y, .flag, .diag_type, .plot_type, ...,
                    .labels = NULL, .geom = "point", .hat = NULL, .sd = NULL,
                    .se = NULL, .cv = NULL, .alpha_level, .rug = FALSE){

  # Grouping
  gr_by <- enquos(...)

    # Create ion labels from variables
  ion2 <- detect_ion(.x, .diag_type)
  ion1 <- detect_ion(.y)

  # 2D density
  .IC <- twodens(.IC, .x, .y, gr_by, .flag)

  # Filter correct titles
  ttl <- filter(point::names_diag, name == .diag_type)

  # R statistics labels for geom_text
  if (!is.null(.labels)) tb_labs <- stat_labs(.IC, gr_by)

  # Base plot
  if (.plot_type == "static") p <- ggplot(data = .IC, aes(x = !! .x, y = !! .y))

  plot_width <- as.numeric(
    ggplot2::ggplotGrob(p)$widths[1]) *
    ceiling(log2(nrow(distinct(.IC, !!!gr_by)))
            )

  # Facets
  p <- p + facet_wrap(vars(!!! gr_by), scales = "free")

  # Geom for "point" data
  if (.geom == "point") {
    p <- p + geom_point(aes(color = .data$dens, alpha = .data$alpha_sc))
  }

  if (.geom == "hexbin") p <- p + geom_hexbin()

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

  # Model
  if (!is.null(.hat)) {
    p <- p + geom_line(aes(y = !!.hat), color = "black", linetype = 2,
                       size = 0.5)
  }

  # Model uncertainties
  bounds <- list(sd = .sd, se = .se, cv = .cv)

  if (!all(sapply(bounds, is.null))) {
    # Uncertainty bounds for geom_ribbon
    .IC <- ribbon_stat(.IC, .hat, bounds[!sapply(bounds, is.null)],
                       .alpha_level)
    p <- p + geom_ribbon(
      data = .IC,
      aes(ymin = .data$lower, ymax = .data$upper),
      color = "black",
      fill = "transparent",
      linetype = 3,
      size = 0.5
    )
  }


  p + scale_color_gradientn(
    "",
    # limits = c(-1, 1),
    breaks = range(.IC$dens),
    labels =  c("divergent", "confluent"),
    colors = colorspace::diverge_hcl(7, rev = TRUE),
    na.value = "transparent",
    guide = guide_colourbar(ticks = FALSE, barwidth = plot_width)
    ) +
    # scale_color_distiller(
    # "",
    # breaks = range(.IC$dens),
    # labels =  c("divergent", "confluent"),
    # palette = "YlOrRd",
    # na.value = "transparent",
    # guide = guide_colourbar(ticks = FALSE, barwidth = plot_width)
    # ) +
    scale_alpha_identity(guide = FALSE, limits = c(1e-5, 1)) +
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
    theme(legend.position = "top")
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
          ionct = stringr::str_replace("ion (ct/sec)", "ion", ion),
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

extract_diff_dens <- function(IC, dens_x, dens_y, dens_z, x, y) {
  x <- pull(IC, !!x)
  y <- pull(IC, !!y)
  int_x <- findInterval(x, dens_x)
  int_y <- findInterval(y, dens_y)
  comb_int <- cbind(int_x, int_y)
  dens_z[comb_int]
}

# creat ion labels from variables
detect_ion <- function(var, diag_type = ""){

  if (stringr::str_detect(as_name(var), "[[:digit:]]") | diag_type == "CV") {
    stringr::str_split(as_name(var), "[[:punct:]]")[[1]] %>% tail(1)
    } else {
      NULL
      }
}

# density calculation function (https://themockup.blog/posts/2020-08-28-heatmaps-in-ggplot2/)
twodens <- function(IC, x, y, grps, flag){
  # Ranges for alpha and density
  ran_z <- range(pull(count(IC, !!! grps), n))

  # Calculate density for x and y frame
  diff_dens <- group_by(IC, !!! grps) %>%
    mutate(
      # Group-wise ranges x and y
      min_x = min(!! x),
      max_x = max(!! x),
      min_y = min(!! y),
      max_y = max(!! y),
      # Bandwidths
      h_x = if_else(MASS::bandwidth.nrd(!!x) == 0, 0.1, MASS::bandwidth.nrd(!!x)),
      h_y = if_else(MASS::bandwidth.nrd(!!y) == 0, 0.1, MASS::bandwidth.nrd(!!y)),
      # Calculate bins for density
      bin_n = log(ran_z[2]) / log(n()) * 100
      ) %>%
    ungroup() %>%
    tidyr::nest(data = -c(!!! grps, !! flag)) %>%
    mutate(dens = purrr::map(data, nest_dens, x, y), .keep ="unused") %>%
    tidyr::pivot_wider(names_from = flag, values_from = dens) %>%
    mutate(
      dens_x = purrr::map(confluent, list("x")),
      dens_y = purrr::map(confluent, list("y")),
      # Take care off case when no outliers
      divergent = tidyr::replace_na(divergent, list(list(z = 0))),
      # Cubed density to accentuate the outliers density
      dens_z = purrr::map2(confluent, divergent, ~{(.x$z - .y$z) ^ 3}),
      .keep ="unused"
      )

  left_join(
    tidyr::nest(IC, IC = -c(!!! grps)),
    diff_dens,
    by = sapply(grps, as_name)
    ) %>%
    mutate(
      dens =
        purrr::pmap(
          list(IC, dens_x, dens_y, dens_z),
          extract_diff_dens,
          x = x,
          y = y
          )
      ) %>%
    select(!!!grps, IC, dens) %>%
    tidyr::unnest(cols = c(IC, dens)) %>%
    group_by(!!! grps) %>%
    # Adjust alpha
    mutate(alpha_sc = (ran_z[1] / n()) / 4)

}


nest_dens <- function(IC, x, y) {
    x <- pull(IC, !!x)
    y <- pull(IC, !!y)
    ran_x <- c(unique(IC$min_x), unique(IC$max_x))
    ran_y <- c(unique(IC$min_y), unique(IC$max_y))
    # Calculate density
    MASS::kde2d(x, y, h = c(IC$h_x, IC$h_y), n = IC$bin_n,
                lims = c(ran_x, ran_y))
}

stat_labs <- function(.IC, gr_by){

  tb_labs <- distinct(.IC, !!! gr_by, .keep_all = TRUE) %>%
    select(!!!gr_by, starts_with(paste0(point::names_stat_R$name, "_R"))) %>%
    tidyr::pivot_longer(
      -c(!!! gr_by),
      names_to = c("stat", ".value"),
      names_sep = "\\_R\\_"
    ) %>%
    tidyr::unite(col = "value", - c(!!!gr_by, stat), na.rm = TRUE)
  lbs <- purrr::map2(
    tb_labs$stat,
    pull(tb_labs, value),
    ~stat_labeller("R", stat = .x, value = as.numeric(.y), "expr")
  )
  tibble::add_column(tb_labs, lbs = lbs)
}

# stat for uncertainty intervals
ribbon_stat <- function(IC, hat, bound, alpha_level){

    fct_switch <- function(alpha_level, bound){
      switch(
        bound,
        sd = qnorm((1 - alpha_level / 2)),
        se = qt((1 - alpha_level / 2), n() - 1),
        cv = 1
        )
    }
    mutate(
      IC,
      fct = fct_switch(alpha_level, names(bound)),
      lower = !! hat - .data$fct * !!bound[[1]],
      upper = !! hat + .data$fct * !!bound[[1]]
      )

  }


  # if (!is.null(.sd)){
  #   .IC <- mutate(
  #     .IC,
  #     fct  = qnorm((1 - .alpha_level / 2)),
  #     lower = !!.hat - .data$fct * !!.sd,
  #     upper = !!.hat + .data$fct * !!.sd
  #   )
  # }
  # if (!is.null(.se)) {
  #   .IC <- mutate(
  #     .IC,
  #     fct  = qt((1 - .alpha_level / 2), n() - 1),
  #     lower = !!.hat - .data$fct * !!.se,
  #     upper = !!.hat + .data$fct * !!.se
  #   )
  # }
  # if (!is.null(.cv)) {
  #   .IC <- mutate(
  #     .IC,
  #     lower = !!.hat - !!.cv,
  #     upper = !!.hat + !!.cv
  #   )
  # }




#-------------------------------------------------------------------------------
# anim plot parts
# if (.plot_type == "anim") {
#   p <- ggplot(
#     data = .IC,
#     aes(
#       x = !!.x,
#       y = !!.y,
#       color = {{ .flag }},
#       frame = execution,
#       text = eval(text_labs(.diag_type, x = !!.x, y = !!.y, ion1, ion2))
#     )
#   )
# }

# if (.plot_type == "anim")  {
#   p <- plotly::ggplotly(p, tooltip = "text") %>%
#     plotly::animation_opts(frame = 150, transition = 0)
#   return(p)
# }
