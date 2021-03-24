#' Calibration of isotope ratios
#'
#' \code{calib_R} function to convert isotope values back-and-forth between R
#' and delta
#'
#' A fundamental in publishing isotope data is the conversion of isotope the delta
#' formulation comparing the obtained R to that of a standard. These values are
#' reported on a per mill scale. However in IC calculations ratios (R) or
#' fractional abundances (F) are most commonly used. This functions provide an easy way
#' of making these transformations for both the average composition but also
#' fractionations as enrichment factors (alpha and epsilon).
#'
#' @param x  A numeric value or vector.
#' @param reference A character or numeric value or vector.
#' @param isotope A character for the name of the isotope (e.g. \code{"13C"}).
#' @param type Type of conversion, isotope composition or enrichment.
#' @param input The type of input value, R, F, and delta.
#' @param output Desired output value, R, F, and delta.
#' @param y A numeric value or vector of the product when calculating enrichment
#' factors, following the convention; x = R(product)/ y = R(substrate). Both x
#' and y should have the same notation, either R, F or delta, as supplied to
#' the argument input.
#' @return A numeric value or vector.
#'
#' @export
#' @examples
#' # R value of 0.0111 to delta on VPDB scale
#' calib_R(0.0111, reference = 0.011237, type = "composition", input = "R",
#'         output = "delta")
#' # or
#' calib_R(0.0111, reference = "VPDB", isotope = "13C",
#'         type = "composition", input = "R", output = "delta")
#'
#' # Fractional abundance for a given R
#' calib_R(0.0111, reference = 0.011237, type = "composition",
#'         input = "R", output = "F")
#'
#' # Alpha enrichment factors can also be calculated based on two delta/R/F
#' # values, where the second value can be supplied as argument y (i.e.,
#' # the product of a reaction)
#' calib_R(-25, reference = "VPDB", isotope = "13C",
#'         type = "enrichment", input = "delta", output = "alpha", y = -105)
#'
calib_R <- function(x, reference, isotope, type = "composition",
                    input = "R", output = "F", y = NULL){

  if (input == output) stop("Input and output are equal.", call. = FALSE)
  if (type == "enrichment") {
    if (!c(output == "epsilon" | output == "alpha")) {
      stop("For enrichment conversion at least output has to be an enrichment factor.", call. = FALSE)
      } else {
        if(!is.null(y) & c(input == "epsilon" | input == "alpha")) {
           stop("Input epsilon or alpha is only meaningfull when converting between enrichment factors.", call. = FALSE)
          }
        if (is.null(y) & !c(input == "epsilon" | input == "alpha")) {
           stop("For enrichment conversion of one value input has to be alpha or epsilon.", call. = FALSE)
          }
        }
      }


  # Transformations for compositions
  conv_composition <- function(method, value, reference, isotope) {
    switch(
      method,
      R_F = calib_R_F(value),
      R_delta = calib_R_delta(value, reference, isotope),
      F_R = calib_F_R(value),
      F_delta = calib_F_delta(value, reference, isotope),
      delta_R = calib_delta_R(value, reference, isotope),
      delta_F = calib_delta_F(value, reference, isotope)
      )
    }

  # Calculations of enrichment factors for kinetic reactions
  conv_enrichment <- function(method, value){
    switch(
      method,
      epsilon_alpha = calib_epsilon_alpha(value),
      alpha_epsilon = calib_alpha_epsilon(value)
      )
    }

  conv_output <- function(type, input, output, value, reference, isotope){
    switch(
      type,
      composition = conv_composition(
        paste(input, output, sep = "_"),
        value, reference,
        isotope
        ),
      enrichment = conv_enrichment(paste(input, output, sep = "_"), value)
      )
  }

  if (is.null(y)) {

    return(conv_output(type, input, output, x, reference, isotope))

    } else {

      R1 <- conv_output("composition", input, "R", x, reference, isotope)
      R2 <- conv_output("composition", input, "R", y, reference, isotope)
      alpha <- R2 / R1
      if (output == "alpha") return(alpha)
      if (output == "epsilon") return(conv_enrichment("alpha_epsilon", alpha))
    }
  }


#-------------------------------------------------------------------------------
# Not exportet
#-------------------------------------------------------------------------------
# transform functions
calib_R_F <- function(x) x / (1 + x)
calib_R_delta <- function(x, reference, isotope){
  1000 * (x / find_reference(reference, isotope) - 1)
  }
calib_F_R <- function(x) x / 1 - x
calib_F_delta <- function(x, reference, isotope){
  calib_R_delta(calib_F_R(x), reference, isotope)
  }
calib_delta_R <- function(x, reference, isotope){
  R <- find_reference(reference, isotope)
  (x / 1000) * R + R
  }
calib_delta_F <- function(x, reference, isotope){
  calib_R_F(calib_delta_R(x, reference, isotope))
  }

# enrichment transform functions
calib_epsilon_alpha <- function(x)  exp((x / 1000))
calib_alpha_epsilon <- function(x)  1000 * log(x)


find_reference <- function(.reference, .isotope){

  # Reference standard database
  tb_ref <- point::reference_R

  if (is_character(.reference)) {
    if (.reference %in% tb_ref$reference) {
      tb_ref <- filter(tb_ref, .data$reference == .reference)
      if (.isotope %in% tb_ref$isotope) {
        return(pull(filter(tb_ref, .data$isotope == .isotope), .data$value))
        } else {
          stop("Isotope unkown for reference standard, provide numeric value manually", call. = FALSE)
          }
      } else {
        stop("Reference standard unkown, provide numeric value manually", call. = FALSE)
        }
    } else {
    return(.reference)
  }
}

