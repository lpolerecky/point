#' Plotting diagnostics
#'
#' \code{plot_RDiag} function for propagation of uncertainty for single ions
#'
#' The \code{plot_RDiag} visualises the diagnostics
#'
#' @param df A tibble containing processed ion count data.
#' @param method Character string for the type of diagnostics.
#' @param args A list of quosures pertaining to the variables required for a
#' call to stat_R. The function expr_R() can be used to simplify setting this
#' argument.
#' @param ... Variables for grouping.
#'
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} or
#' \code{\link[plotly:ggplotly]{ggplotly}}
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
#' # plotting Cook's D diagnostics
#' plot_diag_R(tb.aug, args = expr_R_stat, file.nm)
#'
gg_IC <- function(.df, .ion1, .ion2, .method, ..., .Xt = Xt.pr, .N = N.pr,
                  .species = species.nm, .t = t.nm, .flag = flag,
                  .labels = NULL, .rep = 1, .plot_type = "static",
                  .plot_iso = FALSE){

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

  # Filter execution
  if (.plot_type == "static") .df <- filter(.df, execution == .rep)

  plot_args <- list2(.df = .df, .x = Xt2, .y = Xt1, .flag = flag,
                     .diag_type = .method, !!! gr_by, .labels = .labels,
                     .hat = hat_Xt1, .error = hat_s_Xt1)

  # Residual leverage plot
  if(.method == "norm_E") {
      plot_args[[".y"]] <- quo(studE)
      plot_args[[".x"]] <- quo(hat_Xi)
      plot_args[[".hat"]] <- NULL
      plot_args[[".error"]] <- NULL
    }

  # Normal QQ plot
  if (.method == "QQ"){
      plot_args[[".y"]] <- quo(RQ)
      plot_args[[".x"]] <- quo(TQ)
      plot_args[[".hat"]] <- quo_updt(plot_args[[".y"]], pre = "hat")
      plot_args[[".error"]] <- quo_updt(plot_args[[".y"]], pre = "hat_e")
    }

  # Scale location plot
  if (.method  == "CV"){
      plot_args[[".y"]] <- quo(studE)
      plot_args[[".x"]] <- hat_Xt1
      plot_args[[".hat"]] <- 0
      plot_args[[".error"]] <- NULL
      plot_args[[".lower"]] <- -3.5
      plot_args[[".upper"]] <- 3.5
  }

  # Execute
  data_env <- env(data = .df)
  p <- eval(expr(gg_base(!!!plot_args)), data_env)

  if (.plot_type == "anim")  {
    p <- p + gganimate::transition_states(
      execution,
      transition_length = 2,
      state_length = 1
      )
    return(p)
    }
  if (.plot_type == "static") return(p)

  }



gg_base <- function(.df, .x, .y, .flag, .diag_type, ..., .labels = NULL,
                    .geom = "point", .hat = NULL, .error = NULL, .lower = NULL,
                    .upper = NULL, .rug = FALSE){

  gr_by <- enquos(...)
  if (str_detect(as_name(.x), "[[:digit:]]")| .diag_type == "CV") {
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
  alpha_sc <- median(count(.df, !!!gr_by)$n) / 1e5
  # Filter correct titles
  ttl <- filter(nm_diag_ttl, nm == .diag_type)


  # Uncertainyt bounds
  # Add bounds to dataframe
  if (!is.null(.error)) {
  .df <- mutate(
    .df,
    lower = !!.hat - 2 * !!.error,
    upper = !!.hat + 2 * !!.error
    )
  .lower <- quo(lower)
  .upper <- quo(upper)
  }

  if (!is.null(.labels)) {
    tb_labs <- distinct(.df, !!!gr_by, .keep_all = TRUE) %>%
      select(!!!gr_by, starts_with(paste0(point::nm_stat_R$nm,"_"))) %>%
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

  p <- ggplot(data = .df, aes(x = !!.x, y = !!.y, color = {{.flag}})) +
    facet_wrap(vars(!!!gr_by), scales = "free")
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
  if (!any(sapply(list(.hat, .lower, .upper), is.null))) {
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
      x = axis_labs(ttl$xaxis, ion2, "static"),
      y = axis_labs(ttl$yaxis, ion1, "static"),
      title = ttl$label
      ) +
    theme_classic() +
    theme(legend.position = "none")
  }


gg_IR <- function(.df, .lag, .acf, .flag, ..., .hat = NULL,
                  .error = NULL){

  ggplot(.df, mapping = aes(x = {{.lag}}, y = {{.acf}}, color = {{.flag}})) +
    geom_hline(aes(yintercept = {{.hat}})) +
    geom_segment(mapping = aes(xend = {{.lag}}, yend = {{.hat}})) +
    geom_hline(aes(yintercept = -{{.error}}), color = "darkblue") +
    geom_hline(aes(yintercept = {{.error}}), color = "darkblue") +
    scale_color_manual(values = c(ggplotColours(2)[2], ggplotColours(2)[1])) +
    facet_wrap(vars(...), scales = "free") +
    ggtitle("ACF plot") +
    theme_classic()
}














# # Combine original/augmented datasets and results of diagnostics
#   df.def <- inner_join(
#     distinct(reduce_diag(ls_df, "df", args), execution, ID, .keep_all = TRUE),
#     df.R,
#     by = c("execution", "ID", sapply(gr_by, as_name))
#   )
#
#   ion1 <- as_name(args[["ion1"]])
#   ion2 <- as_name(args[["ion2"]])
#
# sample fraction if interactive
#   if(type == "interactive") {
#   IDs <- sample_frac(filter(df.def, execution == 1), size = 0.01) %>%
#     select(ID)
#   }
#
#   filter_df <- function(df, plot_type){
#     switch(plot_type,
#            interactive  = call2("semi_join", df, IDs, by = "ID"),
#            static = call2("filter", df, expr(execution == rep))
#            )
#   }
#
#   df.def <- filter_df(df = df.def, plot_type = type)
#   stat_lab  <- filter_df(df = stat_lab , plot_type = type)
#
#   plot_args <- list2(
#     df = df.def,
#     stat = stat_lab,
#     y = Xt1,
#     x = Xt2,
#     ion1 = ion1,
#     ion2 = ion2,
#     plot_title = plot_title,
#     plot_type = type,
#     args = args,
#     iso = iso,
#     isoscale = isoscale,
#     !!! gr_by
#     )
#
# # Residual leverage plot
#   if(plot_title == "norm_E"){
#     plot_args[["y"]] <- quo(studE)
#     plot_args[["x"]] <- quo(hat_Xi)
#   }
#
# # normal QQ plot
#   if(plot_title == "QQ"){
#     plot_args[["y"]] <- quo(RQ)
#     plot_args[["x"]] <- quo(TQ)
#   }
#
# # scale location plot
#   if(plot_title == "CV"){
#     plot_args[["y"]] <- quo(studE)
#     plot_args[["x"]] <- quo(fitted)
#   }
#
# # control with environment
#   data_env <- env(data = df.def)
#   gg <- eval(expr(gg_default(!!!plot_args)), data_env)
#
#   if(type == "static"){
#   return(gg)
#   }
#
#   if(type == "interactive"){
#     return(plotly::ggplotly(gg, tooltip="text") %>%
#              plotly::animation_opts(frame = 150, transition = 0)
#     )
#     }
#
#

#-------------------------------------------------------------------------------
# Not exportet
#-------------------------------------------------------------------------------



# #standard plot
# gg_default <- function(df,
#                        stat_lab,
#                        y,
#                        x,
#                        ion1,
#                        ion2,
#                        plot_title,
#                        plot_type,
#                        args,
#                        iso,
#                        isoscale,
#                        ...
#                        ){
#
#   y <- enquo(y)
#   x <- enquo(x)
#   gr_by <- enquos(...)

  # if (plot_type == "interactive") {
  #
  #   text_vc <- c(
  #     "CooksD" = "Rm",
  #     "Rm" = "Rm",
  #     "QQ" = "QQ",
  #     "norm_E" = "norm_E",
  #     "CV" = "CV"
  #     )
  #
  #   text_type <- text_vc[plot_title]
  #
  #   if (iso) {
  #     df <- mutate(
  #       df,
  #       delta = round(
  #         calib_R(
  #           !! quo_updt(args[["Xt"]], x = "M_R"),
  #           standard = isoscale,
  #           type = "composition",
  #           input = "R",
  #           output = "delta"
  #           ),
  #         1
  #         )
  #       )
  #
  #     R_not <- expr(paste0("δ", ion1, " (‰): ", delta))
  #
  #   }else{
  #     df <- mutate(
  #       df,
  #       R = round(!! quo_updt(args[["Xt"]], x = "M_R"), 4)
  #       )
  #
  #     R_not <- expr(paste0("mean R: ", R))
  #
  #     }
  #
  #   text_expr = lst(
  #     Rm = expr(paste0(ion2, " (ct/sec): ", round(!!x, 0), '\n',
  #                      ion1, " (ct/sec): ", round(!!y, 0), '\n',
  #                      "R :", round(!!y / !!x, 4), '\n',
  #                      eval(R_not)
  #                      )
  #     ),
  #     QQ = expr(paste0("TQ: ", round(!!x, 2), '\n',
  #                      "SQ: ", round(!!y, 2), '\n',
  #                      eval(R_not)
  #                      )
  #     ),
  #     norm_E = expr(paste0("e*: ", round(!!x, 4), '\n',
  #                          "hat: ", round(!!y, 2), '\n',
  #                          eval(R_not)
  #                          )
  #     ),
  #     CV = expr(paste0("fitted value: ", round(!!x, 4), '\n',
  #                      "e*: ", round(!!y, 2), '\n',
  #                      eval(R_not)
  #                      )
  #     )
  #
  #   )
  #
  #   p <-  ggplot(data = df,
  #                aes(y = !!y,
  #                    x = !!x,
  #                    color = flag,
  #                    frame = execution,
  #                    text = eval(text_expr[[text_type]])
  #                    )
  #                )
  #



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
    } else{
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

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360 / n
  hcl(h = (seq(h[1], h[2], length = n)),
      c = 100,
      l = 65
      )
}

