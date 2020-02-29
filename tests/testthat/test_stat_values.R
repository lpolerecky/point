context("Ion Count stat validation")
library(point)

# test datasets in extdata directory
sv1 <- system.file("extdata", "2018-01-19-GLENDON", package = "point")
sv2 <- system.file("extdata", "2020-01-17-TREASURE", package = "point")

#-------------------------------------------------------------------------------
# value checks of raw data output. Cross validation of statistics results with
# the processed values as given by the Cameca softwarre. For this, the function
# read_test is used, whiche extracts the Cameca processed data.
#-------------------------------------------------------------------------------

# function for block specific comparison of Cameca and here procesed data
# block 1 is uncorrected data
# block 4 is yield and deadtime corrected EM data and background corrected FC data
validate_IC <- function(ion1, ion2, block){

  # block additional parameters for correcting data
  if (block == 1) {deadtime <- 0
  thr <- 0}
  if (block == 4) {deadtime <- 44
  thr <- 180}

  # perform the routine point workflow for IC data
  tb.rw <- read_IC(sv1)

  # processing raw ion count data
  tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt, deadtime = deadtime, thr = thr)

  # single ion descriptive an predictive statistics for all measured ions
  tb.Xt <- stat_Xt(tb.pr, Xt.pr, N.pr, file.nm, species.nm) %>%
    # round number as provided comparison data is not accurate enough
             mutate(Ntot_Xt.pr = if_else(Ntot_Xt.pr > 10^6,
                                         round(Ntot_Xt.pr, -2),
                                         Ntot_Xt.pr))

  # descriptive an predictive statistics for 13C/12C ratios
  tb.R <- stat_R(tb.pr, Xt.pr, N.pr, ID = "ID", ion1 = "13C", ion2 = "12C",
                 file.nm, species.nm, latex = FALSE) %>%
    # round number as provided comparison data is not accurate enough
    mutate(M_R_Xt.pr = round(M_R_Xt.pr, 7))

  # ion counts in comparison dataset
  tb.test.Xt <- read_test(sv1, block = block)[[1]] %>%
                  mutate(species.nm =
                           recode(num.mt,
                               !!! rlang::set_names(unique(tb.rw$species.nm),
                                                    nm = unique(tb.rw$num.mt))
                               )) %>%
    left_join(tb.Xt, . , by = c("file.nm", "species.nm"))

  # ratios provide in comparison dataset
  R_vec <- rlang::set_names(outer(unique(tb.rw$species.nm),
                                  unique(tb.rw$species.nm),
                                  paste, sep = "/") %>% as.vector(),
                            nm = outer(unique(tb.rw$num.mt),
                                       unique(tb.rw$num.mt),
                                       paste, sep = "/") %>% as.vector())

  R_vec <- R_vec[names(R_vec) %in% read_test(sv1, block = block)[[2]]$R.nm]

  # and loading the test file
  tb.test.R <- read_test(sv1, block = block)[[2]] %>%
    mutate(R.nm = recode(R.nm, !!! R_vec)) %>%
    # combine datasets
    left_join(tb.R, . ,by = c("file.nm", "R.nm"))

  return(list(tb.test.Xt, tb.test.R))
}

# a validation test for 13C and 12C counts and 13C/12C ratios
ion1 <- "13C"
ion2 <- "12C"

val_IC <- validate_IC(ion1 = ion1, ion2 = ion2, block = 1)

# total ion counts
Xt1 <- val_IC[[1]] %>%
  filter(species.nm == ion1 | species.nm == ion2)

# isotope ratios
R1 <- val_IC[[2]] %>%
  filter(R.nm == paste(ion1, ion2, sep = "/"))


test_that("block 1: Ntot of Xt and mean R from EM data", {
  expect_equal(Xt1 %>%
                 pull(Ntot_Xt.test) %>% round(),
               Xt1 %>%
                 pull(Ntot_Xt.pr)
  )
  expect_equal(R1 %>%
                 pull(M_R_Xt.test),
               R1 %>%
                 pull(M_R_Xt.pr)
  )
})


test_that("block 1: Descriptive and predictive sd of R from EM data", {
  expect_equal(R1 %>%
                 pull(RSeM_R_Xt.test) %>% round(., 2),
               R1 %>%
                 pull(RSeM_R_Xt.pr) %>% round(., 2)
  )
  expect_equal(R1 %>%
                 pull(hat_RSeM_R_Xt.test) %>% round(., 2),
               R1 %>%
                 pull(hat_RSeM_R_Xt.pr) %>% round(., 2)
  )
})


val_IC <- validate_IC(ion1 = "13C", ion2 = "12C", block = 4)

# total ion counts
Xt4 <- val_IC[[1]] %>%
  filter(species.nm == ion1 | species.nm == ion2)

# isotope ratios
R4 <- val_IC[[2]] %>%
  filter(R.nm == paste(ion1, ion2, sep = "/"))


test_that("block 4: Ntot of Xt and mean R from EM data", {
  expect_equal(Xt4 %>%
                 pull(Ntot_Xt.test) %>% round(),
               Xt4 %>%
                 pull(Ntot_Xt.pr)
  )
  expect_equal(R4 %>%
                 pull(M_R_Xt.test),
               R4 %>%
                 pull(M_R_Xt.pr)
  )
})
