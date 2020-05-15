#' Calibration of isotope ratios
#'
#' \code{calib_R} function to convert isotope values back-and-forth between R and delta
#'
#' A fundamental in publishing isotope data is the conversion of isotope the delta
#' formulation comparing the obtained R to that of a standard. These values are
#' reported on a per mill scale. However in IC calculations most often we use
#' ratios (R) or fractional abundances (F). This functions provide an easy way
#' of making these transformaitons for both the average composition but also
#' fractionations as enrichment factors.
#'
#' @param x  A numeric value or vector.
#' @param y  A numeric value or vector when calculating enrichments (`type = enrichment`).
#' @param standard A character or numeric value or vector.
#' @param type Type of conversion, isotope composition or enrichment.
#' @param input The type of input value, R, F or delta.
#' @param output Desired ouput value, R, F or delta.
#' @return A numeric value or vector.
#'
#' @export
#' @examples
#' # R value of 0.0111 to delta on VPDB scale
#' calib_R(0.0111, standard = 0.011237, type = "composition", input = "R", output = "delta")
#' # or
#' calib_R(0.0111, standard = "VPDB", type = "composition", input = "R", output = "delta")
#'
#' # Fractional abundance for a given R
#' calib_R(0.0111, type = "composition", input = "R", output = "F")
#'
#' # Generates warning as input and output are equal
#' calib_R(0.0111, type = "composition", input = "R", output = "R")
#'
#' # Generates warning as standard is unknown
#' calib_R(0.0111, standard = "somestandard", type = "composition", input = "R", output = "delta")
#'
#' # Find R values 60 per mille apart
#' A <- calib_R(5, standard = "VPDB", type = "composition", input = "delta", output = "R")
#' B <- calib_R(-55, standard = "VPDB", type = "composition", input = "delta", output = "R")
calib_R <- function(x,
                    y,
                    standard,
                    type = c("composition", "enrichment"),
                    input = c("R", "F", "delta"),
                    output = c("R", "F", "delta")
                    ){

  if (input == output){stop("Input and output is equal", call. = FALSE)}

  convert_comp <- function(type){
    switch(type,
           R_F = calib_R_F(x),
           R_delta = calib_R_delta(x, standard),
           F_R = calib_F_R(x),
           F_delta = calib_F_delta(x, standard),
           delta_R = calib_delta_R(x, standard),
           delta_F = calib_delta_F(x, standard)
           )
  }

  if (type == "composition") convert_comp(paste(input, output, sep = "_"))
  }

#-------------------------------------------------------------------------------
# Not exportet
#-------------------------------------------------------------------------------
# transform functions
calib_R_F <- function(x){x / (1 + x)}
calib_R_delta <- function(x, standard){1000 * (x / find_standard(standard) - 1)}
calib_F_R <- function(x){x / 1 - x}
calib_F_delta <- function(x, standard){calib_R_delta(calib_F_R(x))}
calib_delta_R <- function(x, standard){(x / 1000) * find_standard(standard) + find_standard(standard)}
calib_delta_F <- function(x, standard){calib_R_F(calib_delta_R(x))}

find_standard <- function(standard){

# standard database
  vc_standard <- c(VPDB = 0.011237)

  if (is_character(standard)){
    if (standard %in% names(vc_standard)){
      standard <- vc_standard %>% purrr::pluck(standard)
      return(standard)
      } else {
        stop("Standard unkown, provide numeric value manually", call. = FALSE)
      }
    } else {
    return(standard)
  }
}
