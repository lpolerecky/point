# Plotting diagnostics
gg_IC <- function(.IC, .ion1, .ion2, ..., .X = NULL, .N = NULL, .flag = NULL,
                  .rep = 1, .plot_args, .return_call = FALSE) {

  # Quoting the call (user-supplied expressions)
  # Grouping
  gr_by <- enquos(...)

  # stat_R variables
  args <- inject_args(
    .IC,
    enquos(.X = .X, .N = .N, .flag = .flag),
    type = c("processed", "group", "diagnostics"),
    check = FALSE
  )

  # Common isotope
  N1 <- quo_updt(args[[".X"]], post = .ion1)
  X1 <- quo_updt(args[[".X"]], post = .ion1)
  # rare isotope
  X2 <- quo_updt(args[[".X"]], post = .ion2)
  # Fitted rare isotope and variance (sigma)
  hat_X1 <- quo_updt(args[[".X"]], pre = "hat", post = .ion1)
  hat_S_N1 <- quo_updt(args[[".N"]], pre = "hat_S", post = .ion1)

  # Filter execution
  if (.plot_args[[".plot_type"]] == "static") {
    .IC <- dplyr::filter(.IC, .data$execution == .rep)
  }

  # Plotting arguments
  plot_args <- rlang::list2(
    .IC = .IC,
    .x = X2,
    .y = X1,
    .flag = args[[".flag"]],
    !!! gr_by,
    .hat = hat_X1,
    .sd = hat_S_N1,
    !!! .plot_args
  )

  # Residual leverage plot
  if(plot_args[[".method"]] == "norm_E") {
    plot_args[[".y"]] <- rlang::parse_quo("studE", env = rlang::caller_env())
    plot_args[[".x"]] <- rlang::parse_quo("hat_Xi", env = rlang::caller_env())
    plot_args[[".hat"]] <- NULL
    plot_args[[".sd"]] <- NULL
  }

  # Normal QQ plot
  if (plot_args[[".method"]] == "QQ"){
    plot_args[[".y"]] <- rlang::parse_quo("RQ", env = rlang::caller_env())
    plot_args[[".x"]] <- rlang::parse_quo("TQ", env = rlang::caller_env())
    plot_args[[".hat"]] <- quo_updt(plot_args[[".y"]], pre = "hat")
    plot_args[[".sd"]] <- NULL
    plot_args[[".se"]] <- quo_updt(plot_args[[".y"]], pre = "hat_e")
  }

  # Scale location plot
  if (plot_args[[".method"]] == "CV"){
    plot_args[[".y"]] <- rlang::parse_quo("studE", env = rlang::caller_env())
    plot_args[[".x"]] <- hat_X1
    plot_args[[".hat"]] <- 0
    plot_args[[".sd"]] <- NULL
    plot_args[[".cv"]] <- 3.5
  }
  # call
  p_call <- rlang::call2(gg_base, !!! plot_args)
  # Execute
  if (isTRUE(.return_call)) p_call else  eval(p_call)
}

# base plot
gg_base <- function(.IC, .x, .y, .flag, .method, .plot_type, ...,
                    .plot_stat = NULL, .geom = "point", .hat = NULL,
                    .sd = NULL, .se = NULL, .cv = NULL, .alpha_level,
                    .rug = FALSE,
                    .plot_outlier_labs =  c("divergent", "confluent")){

  # Grouping
  gr_by <- enquos(...)

  # Create ion labels from variables
  ion2 <- detect_ion(.x, .method)
  ion1 <- detect_ion(.y)

  # 2D density
  .IC <- twodens(.IC, !! .x, !! .y, !!! gr_by, .flag = !! .flag)

  # Filter correct titles
  ttl <- dplyr::filter(point::names_plot, .data$name == .method)

  # R statistics labels for geom_text
  if (!is.null(.plot_stat)) tb_labs <- stat_labs(.IC, gr_by)

  # Base plot
  if (.plot_type == "static") {
    p <- ggplot2::ggplot(
      data = .IC,
      mapping = ggplot2::aes(x = !! .x, y = !! .y)
    )
  }

  # equalize plot legend width
  plot_width <- as.numeric(
    ggplot2::ggplotGrob(p)$widths[1]) *
    ceiling(log2(nrow(dplyr::distinct(.IC, !!!gr_by))))

  # Facets
  p <- p + ggplot2::facet_wrap(ggplot2::vars(!!! gr_by), scales = "free")

  # Geom for "point" data
  if (.geom == "point") {
    if (requireNamespace("MASS", quietly = TRUE)) {
      p <- dens_point(p, .flag, .IC, plot_width, .plot_outlier_labs)
    } else {
      message(paste0("Install package \"MASS\" for a better visualization of",
                     " potential outliers."))
      p <- p + ggplot2::geom_point(alpha = 0.3)
    }

  } else if (.geom == "hex") {
    p <- p + ggplot2::geom_hex()
  }
  if (.rug) p <- p + ggplot2::geom_rug(sides = "tr", alpha = 0.01)
  if (!is.null(.plot_stat)) {
    p <- p +
      ggplot2::geom_text(
        data = dplyr::filter(tb_labs, .data$stat %in% .plot_stat),
        mapping = ggplot2::aes(
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
    p <- p + ggplot2::geom_line(
      mapping = ggplot2::aes(y = !!.hat),
      color = "black",
      linetype = 2,
      size = 0.5
    )
  }

  # Model uncertainties
  bounds <- list(sd = .sd, se = .se, cv = .cv)

  if (!all(sapply(bounds, is.null))) {
    # Uncertainty bounds for geom_ribbon
    .IC <- ribbon_stat(.IC, .hat, bounds[!sapply(bounds, is.null)],
                       .alpha_level)
    p <- p +
      ggplot2::geom_ribbon(
        data = .IC,
        mapping = ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        color = "black",
        fill = "transparent",
        linetype = 3,
        size = 0.5
      )
  }

  p <- p +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(3),
      labels = scales::label_scientific(2)
      ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(3),
      labels = scales::label_scientific(2)
      ) +
    ggplot2::labs(
      x = axis_labs(ttl$xaxis, ion2, .plot_type),
      y = axis_labs(ttl$yaxis, ion1, .plot_type),
      title = ttl$label
      ) +
    ggplot2::theme_classic()+
    ggplot2::theme(legend.position = "top")
  }

# Autocorrelation plot
gg_IR <- function(.IC, .lag, .acf, .flag, ..., .sd = NULL){
  ggplot2::ggplot(
    data = .IC,
    mapping = ggplot2::aes(x = {{.lag}}, y = {{.acf}}, color = {{.flag}})
    ) +
    ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = 0)) +
    ggplot2::geom_segment(mapping = ggplot2::aes(xend = {{.lag}}, yend = 0)) +
    ggplot2::geom_hline(
      mapping = ggplot2::aes(yintercept = -{{.sd}}),
      color = "darkblue"
    ) +
    # ggplot2::geom_hline(
    #   mapping = ggplot2::aes(yintercept = {{.sd}}),
    #   color = "darkblue"
    # ) +
    ggplot2::facet_wrap(ggplot2::vars(...), scales = "free") +
    ggplot2::labs(title = "ACF plot") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none")
}

# axis labels
axis_labs <- function(type, ion, plot_type){
  if (!(type == "ionct" | type ==  "hat_Y")) ion <- NULL
  if (plot_type == "static") {
    switch(
      type,
      ionct = substitute(
        a ~ "(count sec" ^ "-" * ")", list(a = ion_labeller(ion, "expr"))
        ),
      studE = expression("studentized residuals (" * italic(e) ^ "*" * ")"),
      TQ = "Theoretical quantiles",
      SQ = "Sample quantiles",
      Xct = expression(X ~ "(count sec" ^ "-" * ")"),
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

# Create ion labels from variables
detect_ion <- function(var, diag_type = "") {
  if (stringr::str_detect(as_name(var), "[[:digit:]]") |
      diag_type == "CV") {
    stringr::str_split(as_name(var), "[[:punct:]]")[[1]] %>% tail(1)
  } else {
    NULL
  }
}

# calculate statistics for plotting on plot
stat_labs <- function(.IC, gr_by) {

  tb_labs <- dplyr::distinct(.IC, !!! gr_by, .keep_all = TRUE) %>%
    dplyr::select(
      !!! gr_by,
      dplyr::starts_with(paste0(point::names_stat_R$name, "_R"))
    ) %>%
    tidyr::pivot_longer(
      -c(!!! gr_by),
      names_to = c("stat", ".value"),
      names_sep = "\\_R\\_"
    ) %>%
    tidyr::unite(col = "value", - c(!!! gr_by, .data$stat), na.rm = TRUE)

  lbs <- purrr::map2(
    tb_labs$stat,
    dplyr::pull(tb_labs, .data$value),
    ~stat_labeller("R", stat = .x, value = as.numeric(.y), label = "expr")
    )
  # Return with new column
  tibble::add_column(tb_labs, lbs = lbs)
}

# Statistics for uncertainty intervals
ribbon_stat <- function(IC, hat, bound, alpha_level) {
  fct_switch <- function(alpha_level, bound) {
    switch(
      bound,
      sd = qnorm((1 - alpha_level / 2)),
      se = qt((1 - alpha_level / 2), dplyr::n() - 1),
      cv = 1
    )
  }
  dplyr::mutate(
    IC,
    fct = fct_switch(alpha_level, names(bound)),
    lower = !!hat - .data$fct * !!bound[[1]],
    upper = !!hat + .data$fct * !!bound[[1]]
  )
}

# calculate 2D density
dens_point <- function (p, flag, IC, width, plot_outlier_labs) {

  # colors
  div_col <- c("#8E063B", "#BB7784", "#D6BCC0", "#E2E2E2", "#BEC1D4", "#7D87B9",
               "#023FA5") # colorspace::diverge_hcl(7, rev = TRUE)

  # point with color following 2D density
  p <- p + ggplot2::geom_point(
    mapping = ggplot2::aes(color = .data$dens, alpha = .data$alpha_sc)
  )

  if (rlang::is_symbol(rlang::get_expr(flag))) {
    p <- p +
      ggplot2::scale_color_gradientn(
        "",
        breaks = range(IC$dens),
        labels =  plot_outlier_labs,
        colors = div_col,
        na.value = "transparent",
        guide = ggplot2::guide_colourbar(ticks = FALSE, barwidth = width)
      )
  } else {
    p <- p +
      ggplot2::scale_color_distiller(
        "",
        breaks = seq(0, 1, length.out = 100),
        palette = "YlOrRd",
        direction = 1,
        na.value = "transparent",
        guide = "none"
      )
  }
  # alpha
  p + ggplot2::scale_alpha_identity(guide = "none", limits = c(1e-5, 1))
}
