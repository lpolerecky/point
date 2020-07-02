#' Plotting diagnostics
#'
#' \code{plot_RDiag} function for propagation of uncertainty for single ions
#'
#' The \code{plot_RDiag} visualises the diagnostics
#'
#' @param df A tibble containing processed ion count data.
#' @param method Character string for the type of diagnostics. Currently only
#' "standard" is supported, which pertains to the default Cameca software
#' setting.
#' @param args A list of quosures pertaining to the variables required for a
#' call to stat_R. The function expr_R can be used to simplify setting this
#' argument.
#' @param ... Variables for grouping.
#'
#' @return A t\code{\link[ggplot2:ggplot]{ggplot}}
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

plot_diag_R <- function(df,
                        args = expr_R(NULL),
                        ...,
                        path = NULL,
                        device = "png",
                        aug = TRUE ,
                        width.out = 20,
                        height.out = 20
                        ){

  # Xt <- enquo(Xt)
  # N <- enquo(N)
  # species <- enquo(species)
  gr_by <- enquos(...)

# heavy isotope
  Xt1 <- quo_updt(args[["Xt"]], as_name(args[["ion1"]]))
# light isotope
  Xt2 <- quo_updt(args[["Xt"]], as_name(args[["ion2"]]))

# labels for point statistics
  stat <- reduce_diag(df, type = "df", args = args, !!! gr_by) %>%
    group_by(execution) %>%
    tidyr::nest() %>%
    mutate(labels = purrr::map2(data, execution,
                ~stat_select(.x, .y, Xt = !!args[["Xt"]], facets_gr = quos(!!!gr_by))
                                )
           ) %>%
    tidyr::unnest(cols = labels) %>%
    ungroup() %>%
    select(-c(execution, data))

# strip grouping
  results <- select(df[[2]]$results, -c(!!!gr_by))

# data with diagnostics
  df.def <-  cov_R(df[[1]]$df,
                   species = !!args[["species"]],
                   ion1 = as_name(args[["ion1"]]),
                   ion2 = as_name(args[["ion2"]]),
                   !!! gr_by,
                   preserve = TRUE
                   ) %>%
    left_join(., results, by = "ID")

  ion1 <- as_name(args[["ion1"]])
  ion2 <- as_name(args[["ion2"]])

  crs <- gg_default(df.def,
                    stat,
                    y = !!Xt1,
                    x = !!Xt2,
                    hat_y = hat_Y,
                    hat_min = hat_Y - 2 * sigma,
                    hat_max =  hat_Y + 2 * sigma,
                    !!!gr_by,
                    z = flag,
                    title = "Cross plot" ,
                    model = TRUE
                    ) +
    xlab(substitute(""^a * b ~"(ct/sec)",
                    lst(a = as.numeric(gsub("([0-9]+).*$", "\\1", ion2)),
                        b = as.symbol(gsub("^[0-9]+","", ion2))
                        )
                    )
         ) +
    ylab(substitute(""^a * b ~"(ct/sec)",
                    lst(a = as.numeric(gsub("([0-9]+).*$", "\\1", ion1)),
                        b = as.symbol(gsub("^[0-9]+","",ion1))
                       )
                    )
         )

  QQ.norm <- gg_default(df.def,
                        stat,
                        y = RQ,
                        x =TQ,
                        hat_y = hat_RQ,
                        hat_min = hat_RQ_min,
                        hat_max = hat_RQ_max,
                        !!! gr_by,
                        z = flag,
                        title = "Normal QQ plot",
                        model = TRUE
                        ) +
    geom_point(shape = 21, alpha = 0.2) +
    ylab(expression("studentized residuals (" * italic(e)^"*" * ")")) +
    xlab("Theoretical quantiles")

  sc.loc <- gg_default(df.def,
                       stat,
                       y = studE,
                       x = hat_Y,
                       hat_y = 0,
                       hat_min = -3.5,
                       hat_max = 3.5,
                       !!!gr_by,
                       z = flag,
                       title = "Scale-Location plot",
                       model = TRUE
                       ) +
    ylim(-5, 5) +
    ylab(expression("studentized residuals (" * italic(e)^"*" * ")")) +
    xlab(substitute("fitted value (" * hat(""^a * b) * ")",
                    lst(a = as.numeric(gsub("([0-9]+).*$", "\\1", ion1)),
                        b = as.symbol(gsub("^[0-9]+","",ion1))
                    )
    )
    )


  rs.lev <- gg_default(df.def,
                       stat,
                       y = studE,
                       x = hat_Xi,
                       z = flag,
                       hat_y = hat_Y,
                       hat_min = NULL,
                       hat_max = NULL,
                       !!!gr_by,
                       title = "Residual vs Leverage plot" ,
                       model = FALSE
                       ) +
    ylab(expression("studentized residuals (" * italic(e)^"*" * ")")) +
    xlab(expression("hat-values" (italic(h))))


  df.acf <- df.def %>%
    group_by(!!! gr_by) %>%
    select(!!! gr_by, ACF) %>%
    # distinct(!!! gr_by) %>%
    tidyr::unnest(cols = c(ACF)) %>%
    ungroup() %>%
    tidyr::unite(facet_gr, !!!gr_by)

  acf <- ggplot(df.acf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) +
    geom_hline(aes(yintercept = ci_upper), color = "darkblue") +
    geom_hline(aes(yintercept = ci_lower), color = "darkblue") +
    facet_wrap(~facet_gr, scales = "free") +
    ggtitle("ACF plot") +
    theme_classic()


  gg.pl <-lst(`1` = crs,
              `2` = QQ.norm,
              `3` = sc.loc,
              `4` = rs.lev,
              `5` = acf
              )

  gg.sa <- function(nm, x, width.out, height.out){

    ggsave(filename = paste0(path,"Figure", nm ,".", device),
           plot = x,
           height = height.out,
           width = width.out,
           units = "cm"
           )
  }

  if (is.null(path)){

    return(gg.pl)

  }else{

    invisible(mapply(gg.sa,
                     names(gg.pl),
                     gg.pl,
                     width.out = width.out,
                     height.out = height.out
                     )
              )
    return(gg.pl)

    }

}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Stat labels selection
stat_select <- function(df, execution, Xt, facets_gr) {

  Xt <- enquo(Xt)
  # facets_gr <- enquos(facets_gr)

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
              lb = purrr:: map2(!!quo_updt(Xt, x = "RSeM_R"),
                                !!quo_updt(Xt, x = "hat_RSeM_R"),
                                ~stat_lab(a = .x, b = .y, aug)
              ),
              vjust = pos,
              !!! facets_gr
              ) %>%
    tidyr::unnest(cols = lb)

}

# funtion for creation of labels for statistics on the original and augmented dataset
stat_lab <- function(a, b, aug = FALSE){

  return(
    list(substitute(epsilon[bar(R)]^c ==
                      ~ a ~ "(" *
                      hat(epsilon)[bar(R)]^c ==
                      ~ b * ")" ~ "\u2030",
                    list(a = sprintf("%.2f", a),
                         b = sprintf("%.2f", b),
                         c = if(aug == TRUE){"*"}else{""} ))) %>%
      do.call("expression", .) %>%
      factor(x = as.character(a),
             labels = .)
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
                       stat,
                       y,
                       x,
                       z,
                       hat_y,
                       hat_min,
                       hat_max,
                       ...,
                       title,
                       model = FALSE
                       ){

  x <- enquo(x)
  y <- enquo(y)
  z <- enquo(z)
  hat_y <- enquo(hat_y)
  hat_min <- enquo(hat_min)
  hat_max <- enquo(hat_max)
  gr_by <- enquos(...)

  # facets groups
  df <- tidyr::unite(data = df, col = facet_gr, !!!gr_by)
  stat <- tidyr::unite(data = stat, col = facet_gr, !!!gr_by)
  # stat.aug <- stat.aug %>% tidyr::unite(facet_gr, !!!gr_by)

  # model
  if(model == TRUE) {
    gg.model <- lst(

      geom_ribbon(aes(ymin = !!hat_min,
                      ymax = !!hat_max),
                  fill = "green",
                  color = "transparent",
                  alpha = 0.1
      ),
      geom_line(aes(y = !!hat_y,
                    x = !!x),
                color = "black",
                linetype = 2,
                size = 0.5
      )
    )

  }else{
    gg.model <- lst(geom_blank())
  }

  ggplot(data = df, aes(y = !!y, x = !!x, color = !!z)) +
    gg.model +
    geom_point(alpha = 0.05) +
    scale_color_manual("Cook's D",
                       values = c("non-influential" = ggplotColours(2)[2],
                                  "influential" = ggplotColours(2)[1])) +
    guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                     linetype = NULL
    )
    )
    ) +
    geom_rug(sides = "tr", alpha = 0.01) +
    geom_text(data = stat,
              aes(label = lb,
                  vjust = vjust,
                  group = trans
                  ),
              x = Inf,
              y = Inf,
              hjust = 1,
              size = 2,
              inherit.aes = FALSE,
              parse = TRUE
              ) +
    # geom_text(data = stat.aug,
    #           aes(label = lb,
    #               vjust = vjust),
    #           x = Inf,
    #           y = Inf,
    #           hjust = 1,
    #           size = 2,
    #           inherit.aes = FALSE,
    #           parse = TRUE) +
    facet_wrap(~facet_gr, scales = "free") +
    ggtitle(title) +
    theme_classic()
}
