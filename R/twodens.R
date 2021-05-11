#' Plotting highly dens data (inspired by
#' \href{https://themockup.blog/posts/2020-08-28-heatmaps-in-ggplot2/}{Thomas Mock})
#'
#'
#' \code{twodens} function for 2D density calculations of one or two components
#'
#' This function is inspired by code in the blog-post by
#' \href{https://themockup.blog/posts/2020-08-28-heatmaps-in-ggplot2/}{Thomas Mock},
#' which deals with over-plotting of very dens data. This is just the problem
#' that occurs when dealing with SIMS ion count data, which usually entails
#' thousands of data-points (measurements) for just one analysis. The function
#' provided here can deal with density of one or two components, where in the
#' latter case the difference of two 2D density matrices is used to show the
#' area on the plot with a maximum difference.
#'
#' @param .IC A tibble containing ion count data.
#' @param .x The variable plotted on the x-axis.
#' @param .y The variable plotted on the y-axis.
#' @param ... Variables for grouping.
#' @param .flag Variable identifying two components in the data.
#'
#' @return A
#' \code{tibble::\link[tibble:tibble]{tibble}()} with the 2D density of x and y
#' named \code{dens_z}.
#'
#' @export
twodens <- function(.IC, .x, .y, ..., .flag){

  gr_by <- enquos(...)
  flag <- enquo(.flag)
  x <- enquo(.x)
  y <- enquo(.y)

  # Ranges for alpha and density
  ran_z <- range(pull(count(.IC, !!! gr_by), .data$n))

  # Calculate density for x and y frame
  IC_dens <- group_by(.IC, !!! gr_by) %>%
    mutate(
      # Group-wise ranges x and y
      min_x = min(!! x),
      max_x = max(!! x),
      min_y = min(!! y),
      max_y = max(!! y),
      # Bandwidths
      h_x = if_else(
        MASS::bandwidth.nrd(!! x) == 0, 0.1, MASS::bandwidth.nrd(!! x)
        ),
      h_y = if_else(
        MASS::bandwidth.nrd(!! y) == 0, 0.1, MASS::bandwidth.nrd(!! y)
        ),
      # Calculate bins for density
      bin_n = log(ran_z[2]) / log(n()) * 100
      ) %>%
    ungroup() %>%
    tidyr::nest(data = -c(!!! gr_by, !! flag)) %>%
    mutate(dens = purrr::map(.data$data, nest_dens, x, y), .keep ="unused")

  if (is_symbol(get_expr(flag))) {
    IC_dens <- diff_dens(IC_dens, flag)
    } else {
      IC_dens <- tidyr::unnest_wider(IC_dens, .data$dens, names_sep = "_")
      }

  left_join(
    tidyr::nest(.IC, IC = -c(!!! gr_by)),
    IC_dens,
    by = sapply(gr_by, as_name)
    ) %>%
    mutate(
      dens =
        purrr::pmap(
          list(.data$IC, .data$dens_x, .data$dens_y, .data$dens_z),
          extract_dens,
          x = x,
          y = y
          )
        ) %>%
    select(!!!gr_by, .data$IC, .data$dens) %>%
    tidyr::unnest(cols = c(.data$IC, .data$dens)) %>%
    group_by(!!! gr_by) %>%
    # Adjust alpha
    mutate(alpha_sc = (ran_z[1] / n()) / 4)

}

#-------------------------------------------------------------------------------
# Not exportet
#-------------------------------------------------------------------------------
diff_dens <- function(IC, flag) {
  tidyr::pivot_wider(IC, names_from = !! flag, values_from = .data$dens) %>%
    mutate(
      dens_x = purrr::map(.data$confluent, list("x")),
      dens_y = purrr::map(.data$confluent, list("y")),
      # Take care off case when no outliers
      divergent = tidyr::replace_na(.data$divergent, list(list(z = 0))),
      # Cubed density to accentuate the outliers density
      dens_z =
        purrr::map2(
          .data$confluent,
          .data$divergent,
          ~{(.x$z - .y$z) ^ 3}
          ),
      .keep ="unused"
      )
}

extract_dens <- function(IC, dens_x, dens_y, dens_z, x, y) {
  x <- pull(IC, !!x)
  y <- pull(IC, !!y)
  int_x <- findInterval(x, dens_x)
  int_y <- findInterval(y, dens_y)
  comb_int <- cbind(int_x, int_y)
  dens_z[comb_int]
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
