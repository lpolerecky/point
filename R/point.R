#' point: A package for reading, processing, and analysing raw ion count data.
#'
#' The point package provides function for reading processing and analysing
#' pulsed ion count data.
#'
#' The point functions provide a convenient way to work with and evaluate pulsed
#' ion count data. The \code{read_*} functions can load and collate raw ion
#' count data in a \code{tibble::\link[tibble:tibble]{tibble}()}. Systematic
#' biases that are introduced by the ion detection devices can be corrected with
#' the \code{cor_*} functions. The precision of ion count data can be checked
#' with the \code{stat_*} functions. The \code{diag_*} functions provide a
#' convenient way to measure the influence of individual measurements
#' (\emph{N}) on an \emph{n}-series of measurements (or analysis).
#'
#' @import dplyr
#' @importFrom stats cov sd AIC acf anova approx coef dnorm fitted pchisq ppoints qchisq qnorm qt quantile rpois t.test
#' @importFrom ggplot2 ggplot aes guide_colourbar facet_wrap geom_rug geom_hex geom_text geom_point geom_line geom_hline geom_ribbon geom_segment scale_color_manual scale_color_distiller scale_alpha_identity scale_y_continuous scale_x_continuous labs theme_classic theme scale_color_gradientn
#' @importFrom mgcv gamm s
#' @importFrom nlme lme gls
#' @importFrom utils data tail
#'
#' @docType package
#' @name point
NULL

# @importFrom readr read_lines read_tsv read_table read_table2 parse_number parse_guess
# @importFrom  stringr str_detect str_which str_trim str_split str_replace str_subset str_replace_all str_c
# @importFrom tidyr separate separate_rows drop_na pivot_wider
# @importFrom tibble tibble is_tibble
