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
#' # expression for R_stat()
#' expr_R_stat <- expr_R(Xt = "Xt.pr",
#'                       N = "N.pr",
#'                       species = "species.nm",
#'                       ion1 = "13C",
#'                       ion2 = "12C"
#'                       )
#'
#' # diagnostics ion count data for isotope ratios
#' tb.aug <- diag_R(tb.pr,
#'                  method = "CooksD",
#'                  args = expr_R_stat,
#'                  reps = 2,
#'                  file.nm
#'                  )
#'
#' # plotting Cook's D diagnostics
#' plot_diag_R(tb.aug, args = expr_R_stat, file.nm)
#'
plot_diag_R <- function(ls_df,
                        args = expr_R(NULL),
                        ...,
                        plot_title,
                        type,
                        rep = "1",
                        iso = TRUE,
                        isoscale = NULL
                        ){

  gr_by <- enquos(...)

  reps <- seq_along(ls_df)
  df.reps <- reps[min(reps):max(reps)-1]
  results.reps <- reps[-1]

# heavy isotope
  Xt1 <- quo_updt(args[["Xt"]], as_name(args[["ion1"]]))
# light isotope
  Xt2 <- quo_updt(args[["Xt"]], as_name(args[["ion2"]]))

# iso stats dataset
  df.R <- reduce_diag(ls_df, "results", args, !!!gr_by)

# labels for point statistics
  stat_lab <- df.R %>%
    group_by(execution) %>%
    tidyr::nest() %>%
    mutate(labels = purrr::map2(
      data,
      execution,
      ~stat_select(.x,
                   .y,
                   Xt = !!args[["Xt"]],
                   ion1 = as_name(args[["ion1"]]),
                   facets_gr = quos(!!!gr_by),
                   iso = iso,
                   isoscale = isoscale
                   )
      )
           ) %>%
    tidyr::unnest(cols = labels) %>%
    ungroup() %>%
    select(-c(data))


  # Combine original/augmented datasets and results of diagnostics
  df.def <- inner_join(
    reduce_diag(ls_df, "df", args, !!!gr_by),
    df.R,
    by = c("execution", "ID", sapply(gr_by, as_name))
  )

  ion1 <- as_name(args[["ion1"]])
  ion2 <- as_name(args[["ion2"]])

# sample fraction if interactive
  if(type == "interactive") {
  IDs <- sample_frac(filter(df.def, execution == 1), size = 0.1) %>%  select(ID)
  }

  filter_df <- function(df, plot_type){
    switch(plot_type,
           interactive  = call2("semi_join", df, IDs, by = "ID"),
           static = call2("filter", df, expr(execution == rep))
           )
  }

  df.def <- filter_df(df = df.def, plot_type = type)
  stat_lab  <- filter_df(df = stat_lab , plot_type = type)

  plot_args <- list2(df = df.def,
                     stat = stat_lab,
                     y = Xt1,
                     x = Xt2,
                     ion1 = ion1,
                     ion2 = ion2,
                     plot_title = plot_title,
                     plot_type = type,
                     args = args,
                     iso = iso,
                     isoscale = isoscale,
                     !!! gr_by
                     )

# Residual leverage plot
  if(plot_title == "norm_E"){
    plot_args[["y"]] <- quo(studE)
    plot_args[["x"]] <- quo(hat_Xi)
  }

# normal QQ plot
  if(plot_title == "QQ"){
    plot_args[["y"]] <- quo(RQ)
    plot_args[["x"]] <- quo(TQ)
  }

# scale location plot
  if(plot_title == "CV"){
    plot_args[["y"]] <- quo(studE)
    plot_args[["x"]] <- quo(fitted)
  }

# control with environment
  data_env <- env(data = df.def)
  gg <- eval(expr(gg_default(!!!plot_args)), data_env)

  if(type == "static"){
  return(gg)
  }

  if(type == "interactive"){
    return(plotly::ggplotly(gg, tooltip="text")  %>%
             plotly::animation_opts(frame = 150, transition = 0)
    )
    }

}

#-------------------------------------------------------------------------------
# Not exportet
#-------------------------------------------------------------------------------

# Stat labels selection
stat_select <- function(df,execution, Xt, ion1, facets_gr, iso, isoscale) {

  Xt <- enquo(Xt)

  if (execution == 1) {
    label <- "original"
    pos <- 2.5
    aug <- FALSE
    }else{
    label <- "augmented"
    pos <- 3.5
    aug <- TRUE
    }

  df %>%
    distinct(!!! facets_gr, .keep_all = TRUE) %>%
    transmute(trans = label,
              lb = purrr:: map2(!!quo_updt(Xt, x = "M_R"),
                                !!quo_updt(Xt, x = "RSeM_R"),
                                ~stat_lab(a = .x,
                                          b = .y,
                                          ion1 = ion1,
                                          aug,
                                          iso = iso,
                                          isoscale = isoscale)
              ),
              vjust = pos,
              !!! facets_gr
              ) %>%
    tidyr::unnest(cols = lb)

}

# function for creation of labels for statistics on the original and augmented dataset
stat_lab <- function(a, b, ion1, aug = FALSE, iso = TRUE, isoscale = "VPDB"){

# turn R into delta value
  if (iso) {
    a <- calib_R(a,
                 standard = isoscale,
                 type = "composition",
                 input = "R",
                 output = "delta"
                 )
    num_iso <- as.numeric(gsub("([0-9]+).*$", "\\1", ion1))
    el_iso <- as.symbol(gsub("^[0-9]+","",ion1))

    expr_R <- list(substitute(delta^a*b*""^c ==
                   ~ d ~ "(" *
                   epsilon[bar(R)]^c ==
                   ~ e * "\u2030"~
                   f*
                   ")",
                   list(a = num_iso,
                        b = el_iso,
                        c = if(aug == TRUE){"*"}else{""},
                        d = sprintf("%.2f", a),
                        e = sprintf("%.2f", b),
                        f = isoscale
                        )
                   )
    )
  } else {

    expr_R <- list(substitute(bar(R)*""^a ==
                                ~ b ~ "(" *
                                epsilon[bar(R)]^a ==
                                ~ c * "\u2030)",
                              list(a = if(aug == TRUE){"*"}else{""},
                                   b = sprintf("%.2f", a),
                                   c = sprintf("%.2f", b)
                                   )
                              )
                  )
    }

  return(
    expr_R %>%
      do.call("expression", .) %>%
      factor(x = as.character(a),
             labels = .
             )
  )
}

# ggplot color generator
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100,
      l = 65)
}


#standard plot
gg_default <- function(df,
                       stat_lab,
                       y,
                       x,
                       ion1,
                       ion2,
                       plot_title,
                       plot_type,
                       args,
                       iso,
                       isoscale,
                       ...
                       ){

  y <- enquo(y)
  x <- enquo(x)
  gr_by <- enquos(...)

  if (plot_type == "interactive") {


    text_vc <- c(
      "CooksD" = "Rm",
      "Rm" = "Rm",
      "QQ" = "QQ",
      "norm_E" = "norm_E",
      "CV" = "CV"
      )

    text_type <- text_vc[plot_title]

    if (iso) {
      df <- mutate(df,
                   delta = round(
                     calib_R(!! quo_updt(args[["Xt"]], x = "M_R"),
                             standard = isoscale,
                             type = "composition",
                             input = "R",
                             output = "delta"),
                     1
                   )
      )

      R_not <- expr(paste0("δ", ion1, " (‰): ", delta))
    }else{
      df <- mutate(df,
                   R = round(!! quo_updt(args[["Xt"]], x = "M_R"), 4)
                   )

      R_not <- expr(paste0("R: ", R))
      }

    text_expr = lst(
      Rm = expr(paste0(ion2, " (ct/sec): ", !!x, '\n',
                       ion1, " (ct/sec): ", !!y, '\n',
                       eval(R_not)
                       )
      ),
      QQ = expr(paste0("TQ: ", round(!!x, 2), '\n',
                       "SQ: ", round(!!y, 2), '\n',
                       eval(R_not)
                       )
      ),
      norm_E = expr(paste0("e*: ", round(!!x, 4), '\n',
                           "hat: ", round(!!y, 2), '\n',
                           eval(R_not)
                           )
      ),
      CV = expr(paste0("fitted value: ", round(!!x, 4), '\n',
                       "e*: ", round(!!y, 2), '\n',
                       eval(R_not)
                       )
      )

    )

    p <-  ggplot(data = df,
                 aes(y = !!y,
                     x = !!x,
                     color = flag,
                     frame = execution,
                     text = eval(text_expr[[text_type]])
                     )
                 )
    } else {
      p <-  ggplot(data = df, aes(y = !!y, x = !!x, color = flag))
    }
  if (plot_type == "interactive") {
    p <- p + geom_point(alpha = 0.5)
  } else {
    p <- p + geom_hex(alpha = 0.5)
    # geom_point(alpha = 0.05)
  }


  if (all(c("hat_min", "hat_max", "hat_Y") %in% colnames(df)) & plot_type == "static"){
    p <- p + geom_ribbon(aes(ymin = hat_min,
                             ymax = hat_max
                            ),
                         color = "black",
                         fill = "transparent",
                         #alpha = 0.1
                         linetype = 3,
                         size = 0.5
                         )
    }


  if (all(c("hat_min", "hat_max", "hat_Y") %in% colnames(df)) & plot_type == "static"){

    p <- p + geom_line(aes(y = hat_Y),
                       color = "black",
                       linetype = 2,
                       size = 0.5)
    }

    p <- p +
      scale_color_manual("",
                       values = c(
                         "good" = ggplotColours(2)[2],
                         "bad" = ggplotColours(2)[1]
                                  ),
                      ) +
      scale_y_continuous(breaks = scales::extended_breaks(), labels = scales::label_scientific()) +
      scale_x_continuous(breaks = scales::extended_breaks(), labels = scales::label_scientific()) +
      geom_rug(sides = "tr", alpha = 0.01)

    if (plot_type == "static") {
      if (!is.null(stat_lab)) {
      p <- p +  geom_text(data = stat_lab,
                           aes(label = lb,
                               vjust = vjust,
                               group = trans
                               ),
                           x = Inf,
                           y = Inf,
                           hjust = 1.2,
                           size = 1.5,
                           inherit.aes = FALSE,
                           parse = TRUE
                           )
    }
    }

    p <- p + facet_wrap(vars(!!! gr_by), scales = "free") #, labeller = label_wrap_gen(width = 15))

    p + labs(!!! lab_updater(plot_title, ion1, ion2, plot_type)) +
    theme_classic() +
    theme(legend.position = "none")
}


lab_updater <- function(method, ion1, ion2, plot_type){

  plot_names <- c(
    "Rm" = "Linear R model (Cameca)",
    "norm_E" = "Residual vs. Leverage",
    "CooksD" = "Linear R model (Cook's D)",
    "QQ" = "Normal QQ plot",
    "CV" = "Scale-location plot",
    "QSA" = "QSA test",
    "timeseries" = "timeseries"
    )

  plot_labs <- lst(
        Rm = lst(
           x = axis_labs("ionct", ion2, plot_type),
           y = axis_labs("ionct", ion1, plot_type),
           title = unname(plot_names[method])
         ),
        norm_E = lst(
           x = axis_labs("studE", NULL, plot_type),
           y = axis_labs("hat_Xi", NULL,  plot_type),
           title = unname(plot_names[method])
         ),
        CooksD = lst(
          x = axis_labs("ionct", ion2, plot_type),
          y = axis_labs("ionct", ion1, plot_type),
          title = unname(plot_names[method])
         ),
        QQ = lst(
          x = axis_labs("TQ", NULL, plot_type),
          y = axis_labs("SQ", NULL, plot_type),
          title = unname(plot_names[method])
        ),
        CV = lst(
          x = axis_labs("hat_Y", ion1, plot_type),
          y = axis_labs("studE", NULL, plot_type),
          title = unname(plot_names[method])
        ),
        QSA = lst(
          x = axis_labs("ionct", ion2, plot_type),
          y = axis_labs("R", NULL, plot_type),
          title = unname(plot_names[method])
        ),
        timeseries = lst(
          x = axis_labs("time", NULL, plot_type),
          y = axis_labs("Xct", NULL, plot_type),
          title = unname(plot_names[method])
        )
        )

  plot_labs[[method]]
}

axis_labs <- function(type, ion = NULL, plot_type){

  if (plot_type == "static"){

  switch(type,
         ionct = substitute(""^a * b ~"(ct/sec)",
                            lst(a = as.numeric(gsub("([0-9]+).*$", "\\1", ion)),
                                b = as.symbol(gsub("^[0-9]+","",ion))
                                )
                            ),
         studE = expression("studentized residuals (" * italic(e)^"*" * ")"),
         TQ = "Theoretical quantiles",
         SQ = "Sample quantiles",
         Xct = "X (ct/sec)",
         time = "time (sec)",
         hat_Y = substitute("fitted value (" * hat(""^a * b) * ")",
                            lst(a = as.numeric(gsub("([0-9]+).*$", "\\1", ion)),
                                b = as.symbol(gsub("^[0-9]+","",ion))
                            )
                            ),
         hat_Xi = expression("hat-values" (italic(h))),
         R = "R"
        )

    } else{

    switch(type,
           ionct = str_replace("ion (ct/sec)", "ion", ion),
           studE = "studentized residuals",
           TQ = "Theoretical quantiles",
           SQ = "Sample quantiles",
           hat_Y = "fitted value ",
           hat_Xi = "hat-values"
           )



  }
}
