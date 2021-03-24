# Plotting diagnostics
gg_IC <- function(.IC, .ion1, .ion2, .method, .plot_type, ..., .X = Xt.pr,
                  .N = N.pr, .species = species.nm, .t = t.nm, .flag = flag,
                  .labels = NULL, .rep = 1, .alpha_level, .plot_iso = FALSE){

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)

  # stat_R variables
  args <- enquos(.X = .X, .N = .N, .species = .species, .t = .t)

  flag <-  enquo(.flag)
  # Common isotope
  N1 <- quo_updt(args[[".X"]], post = .ion1)
  X1 <- quo_updt(args[[".X"]], post = .ion1)
  # rare isotope
  X2 <- quo_updt(args[[".X"]], post = .ion2)
  # Fitted rare isotope and variance (sigma)
  hat_X1 <- quo_updt(args[[".X"]], pre = "hat", post = .ion1)
  hat_S_N1 <- quo_updt(args[[".N"]], pre = "hat_S", post = .ion1)

  # Filter execution
  if (.plot_type == "static") .IC <- filter(.IC, .data$execution == .rep)

  # Plotting arguments
  plot_args <- list2(.IC = .IC, .x = X2, .y = X1, .flag = flag,
                     .diag_type = .method, .plot_type = .plot_type, !!! gr_by,
                     .labels = .labels, .hat = hat_X1, .sd = hat_S_N1 ,
                     .alpha_level = .alpha_level)
  # environment
  data_env <- env(data = .IC)
  # Residual leverage plot
  if(.method == "norm_E") {
      plot_args[[".y"]] <- parse_quo("studE", env = data_env)
      plot_args[[".x"]] <- parse_quo("hat_Xi", env = data_env)
      plot_args[[".hat"]] <- NULL
      plot_args[[".sd"]] <- NULL
    }

  # Normal QQ plot
  if (.method == "QQ"){
      plot_args[[".y"]] <- parse_quo("RQ", env = data_env)
      plot_args[[".x"]] <- parse_quo("TQ", env = data_env)
      plot_args[[".hat"]] <- quo_updt(plot_args[[".y"]], pre = "hat")
      plot_args[[".sd"]] <- NULL
      plot_args[[".se"]] <- quo_updt(plot_args[[".y"]], pre = "hat_e")
    }

  # Scale location plot
  if (.method  == "CV"){
      plot_args[[".y"]] <- parse_quo("studE", env = data_env)
      plot_args[[".x"]] <- hat_X1
      plot_args[[".hat"]] <- 0
      plot_args[[".sd"]] <- NULL
      plot_args[[".cv"]] <- 3.5
  }

  # Execute
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
  .IC <- twodens(.IC, !!.x, !!.y, !!! gr_by, .flag = !! .flag)

  # Filter correct titles
  ttl <- filter(point::names_diag, .data$name == .diag_type)

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
  if (.geom == "point") p <- dens_point(p, .flag, .IC, plot_width)

  if (.geom == "hex") p <- p + geom_hex()

  if (.rug) p <- p + geom_rug(sides = "tr", alpha = 0.01)

  if (!is.null(.labels)) {
    p <- p + geom_text(
      data = filter(tb_labs, .data$stat %in% .labels),
      aes(
        x = -Inf,
        y = Inf,
        label = .data$lbs
        ),
      vjust = "inward",
      hjust = "inward",
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

  p <- p +
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
    theme_classic()+
    theme(legend.position = "top")

  }


gg_IR <- function(.df, .lag, .acf, .flag, ..., .sd = NULL){

  ggplot(.df, mapping = aes(x = {{.lag}}, y = {{.acf}}, color = {{.flag}})) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = {{.lag}}, yend = 0)) +
    geom_hline(aes(yintercept = -{{.sd}}), color = "darkblue") +
    geom_hline(aes(yintercept = {{.sd}}), color = "darkblue") +
    facet_wrap(vars(...), scales = "free") +
    labs(title = "ACF plot") +
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


# creat ion labels from variables
detect_ion <- function(var, diag_type = ""){

  if (stringr::str_detect(as_name(var), "[[:digit:]]") | diag_type == "CV") {
    stringr::str_split(as_name(var), "[[:punct:]]")[[1]] %>% tail(1)
    } else {
      NULL
      }
}

stat_labs <- function(.IC, gr_by){

  tb_labs <- distinct(.IC, !!! gr_by, .keep_all = TRUE) %>%
    select(!!! gr_by, starts_with(paste0(point::names_stat_R$name, "_R"))) %>%
    tidyr::pivot_longer(
      -c(!!! gr_by),
      names_to = c("stat", ".value"),
      names_sep = "\\_R\\_"
      ) %>%
    tidyr::unite(col = "value", - c(!!! gr_by, .data$stat), na.rm = TRUE)
  lbs <- purrr::map2(
    tb_labs$stat,
    pull(tb_labs, .data$value),
    ~stat_labeller("R", stat = .x, value = as.numeric(.y), label = "expr")
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


dens_point <- function (p, flag, IC, width) {
  # colors
  div_col <- c("#8E063B", "#BB7784", "#D6BCC0", "#E2E2E2", "#BEC1D4", "#7D87B9",
               "#023FA5") # colorspace::diverge_hcl(7, rev = TRUE)
  p <- p + geom_point(aes(color = .data$dens, alpha = .data$alpha_sc))
  if (is_symbol(get_expr(flag))) {
    p <- p + scale_color_gradientn(
        "",
        breaks = range(IC$dens),
        labels =  c("divergent", "confluent"),
        colors = div_col,
        na.value = "transparent",
        guide = guide_colourbar(ticks = FALSE, barwidth = width)
        )
      } else {
      p <- p + scale_color_distiller(
        "",
        breaks = seq(0, 1, length.out = 100),
        palette = "YlOrRd",
        direction = 1,
        na.value = "transparent",
        guide = FALSE
        )
        }
  # alpha
  p <- p + scale_alpha_identity(guide = FALSE, limits = c(1e-5, 1))
  return(p)
}
