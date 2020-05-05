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
#' # plotting Cook's D diagnostics
#' plot_RDiag(xcp, Xt.pr, N.pr, species.nm, "13C", "12C", file.nm)
#'
#'
plot_RDiag <- function(df, Xt, N, species, ion1, ion2, ..., path = NULL, device = "png", aug = TRUE , width.out = 20, height.out = 20){

  Xt <- enquo(Xt)
  N <- enquo(N)
  species <- enquo(species)
  gr_by <- enquos(...)

# heavy isotope
  Xt1 <- quo_updt(Xt, ion1)
# light isotope
  Xt2 <- quo_updt(Xt, ion2)

  # ID for connecting flag to original dataframe
  df <- ID_builder(df, !! species, !!! gr_by)

  df.def <- diag_R(df,
                   method = "CooksD",
                   args = expr_R(Xt = as_name(Xt),
                                 N = as_name(N),
                                 species = as_name(species),
                                 ion1 = ion1,
                                  ion2 = ion2),
                   !!! gr_by ,
                   output = "complete")


# Stat labels selection
  stat_select <- function(df, label, pos, facets_gr, aug) {

    df %>%
      distinct(!!! facets_gr, .keep_all = TRUE) %>%
      transmute(data = label,
                lb = purrr:: map2(!!quo_updt(Xt, x = "RSeM_R"),
                                  !!quo_updt(Xt, x = "hat_RSeM_R"),
                                  ~stat_lab(a = .x, b = .y, aug)),
                vjust = pos,
                file.nm = file.nm) %>%
      tidyr::unnest(cols = lb)


  }

# Stat labels original dataset
  lb.def <- stat_select(df.def, "original", 2.5, facets_gr = gr_by, aug = FALSE)

# Stat labels for augmented dataset
  lb.aug <- left_join(df, df.def %>% select(ID, CooksD_lab), by = "ID") %>%
                filter(CooksD_lab == "non-influential") %>%
                stat_R(Xt = !! Xt,
                       N = !! N,
                       species = !! species,
                       ion1 = ion1,
                       ion2 = ion2,
                       !!! gr_by) %>%
                stat_select("augmented", 3.6, facets_gr = gr_by, aug = TRUE)

# standard ggplot
  gg_default <- function(df, stat.def, stat.aug, y, x, z, hat_y, ... ,title, model = FALSE){

    x <- enquo(x)
    y <- enquo(y)
    z <- enquo(z)
    hat_y <- enquo(hat_y)
    gr_by <- enquos(...)

    # facets groups
    df <- df %>% tidyr::unite(facet_gr, !!!gr_by)
    stat.def <- stat.def %>% tidyr::unite(facet_gr, !!!gr_by)
    stat.aug <- stat.aug %>% tidyr::unite(facet_gr, !!!gr_by)

# ggplot color generator
    ggplotColours <- function(n = 6, h = c(0, 360) + 15){
      if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
      hcl(h = (seq(h[1], h[2], length = n)), c = 100,
          l = 65)
    }

# model
  if(model == TRUE) {gg.model <- lst(geom_line(aes(y = !!hat_y, x = !!x),
                                                 color = "green" ,
                                                 alpha = 0.5,
                                                 size = 1))
  }else{
    gg.model <- lst(geom_blank())
    }

    ggplot(data = df, aes(y = !!y, x = !!x, color = !!z)) +
      geom_point(alpha = 0.05) +
      scale_color_manual("Cook's D",
                         values = c("non-influential" = ggplotColours(2)[2],
                                    "influential" = ggplotColours(2)[1])) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, linetype = NULL))) +
      geom_rug(sides = "tr", alpha = 0.01) +
      gg.model +
      geom_text(data = stat.def,
                aes(label = lb,
                    vjust = vjust),
                x = Inf,
                y = Inf,
                hjust = 1,
                size = 2,
                inherit.aes = FALSE,
                parse = TRUE) +
      geom_text(data = stat.aug,
                aes(label = lb,
                    vjust = vjust),
                x = Inf,
                y = Inf,
                hjust = 1,
                size = 2,
                inherit.aes = FALSE,
                parse = TRUE) +
      facet_wrap(~facet_gr, scales = "free") +
      ggtitle(title) +
      theme_classic()


  }


  crs <- gg_default(df.def, lb.def, lb.aug,  y = !!Xt1, x = !!Xt2,
                    hat_y = hat_Y, !!!gr_by,
                    z = CooksD_lab,
                    title = "Cross plot" , model = TRUE) +
           xlab(substitute(""^a * b ~"(ct/sec)",
                           lst(a = as.numeric(gsub("([0-9]+).*$", "\\1", ion2)),
                               b = as.symbol(gsub("^[0-9]+","",ion2))))) +
           ylab(substitute(""^a * b ~"(ct/sec)",
                           lst(a = as.numeric(gsub("([0-9]+).*$", "\\1", ion1)),
                               b = as.symbol(gsub("^[0-9]+","",ion1)))))

  rs.fit <- gg_default(df.def, lb.def, lb.aug,  y = E, x = hat_Y,
                       hat_y = NULL, !!!gr_by,
                       z = CooksD_lab,
                       title = "Residuals vs Fitted plot" , model = FALSE) +
              ylab(expression("residuals (" * italic(e) * ")"))+
              xlab(expression("fitted value (" * hat(y) * ")"))

  sc.loc <- gg_default(df.def, lb.def, lb.aug,  y = studE, x = hat_Y,
                       hat_y = NULL, !!!gr_by,
                       z = CooksD_lab,
                       title = "Scale-Location plot" , model = FALSE) +
              ylab(expression("studentized residuals (" * italic(e)^"*" * ")")) +
              xlab(expression("fitted value (" * hat(y) * ")"))

  rs.lev <- gg_default(df.def, lb.def, lb.aug,  y = studE, x = hat_Xi,
                       z = CooksD_lab,
                       hat_y = NULL, !!!gr_by,
                       title = "Residual vs Leverage plot" , model = FALSE) +
              ylab(expression("studentized residuals (" * italic(e)^"*" * ")")) +
              xlab(expression("hat-values" (italic(h))))


  gg.pl <-lst(`1` = crs, `2` = rs.fit, `3`= sc.loc, `4`=rs.lev)

  gg.sa <- function(nm, x, width.out, height.out){

    ggsave(filename = paste0(path,"Figure", nm ,".", device),
           plot = x,
           height = height.out,
           width = width.out,
           units = "cm")
  }

  if (is.null(path)){

    return(gg.pl)

  }else{

    invisible(mapply(gg.sa, names(gg.pl), gg.pl,
                     width.out = width.out,
                     height.out = height.out))
    return(gg.pl)

    }

}


