context("Reading Ion Count Data")
library(point)

# test datasets in extdata directory
sv1 <- system.file("extdata", "2018-01-19-GLENDON", package = "point")
sv2 <- system.file("extdata", "2020-01-17-TREASURE", package = "point")

#-------------------------------------------------------------------------------
# variable class checks and compatibility of txt and stat files
#-------------------------------------------------------------------------------
# error when directory argument class is wrong
test_that("read_IC numeric input error", {
  expect_error(read_validator(1))
  expect_error(read_validator(1L))
})

# warning when empty text files exist
test_that("empty txt file", {
  expect_warning(read_validator(
    system.file("extdata",
                "2018-01-19-GLENDON",
                package = "point")
  ))
})

# warning when zero count exist
test_that("zero counts", {
  expect_warning(
    zeroCt(
      read_IC(sv1), N.rw, "13C", "12C", file.nm)
    )
  }
  )

# testing the class of the retrieved data
test_that("read_IC creates a tibble", {
  expect_is(read_IC(
    system.file("extdata", "2020-01-17-TREASURE", package = "point")
              ), "tbl_df")
})

#-------------------------------------------------------------------------------
# value checks of raw data output. Cross validation of statistics results with
# the processed values as given by the Cameca softwarre. For this, the function
# read_test is used, whiche extracts the Cameca processed data.
#-------------------------------------------------------------------------------


#----------------------------------
# and loading the test file block 1 (uncorrected data)
#----------------------------------

# perform the routine point workflow for EM data but without correcting the data
# raw data containing 13C and 12C counts on carbonate
tb.rw <- read_IC(sv1)

# processing raw ion count data (deadtime to zero, no correction for deadtime)
tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt, deadtime = 44, yield = TRUE)

# single ion descriptive an predictive statistics for all measured ions
tb.Xt <- stat_Xt(tb.pr, Xt.pr, N.pr, file.nm, species.nm)

# descriptive an predictive statistics for 13C/12C ratios
tb.R <- stat_R(tb.pr, Xt.pr, N.pr, ID = "ID", ion1 = "13C", ion2 = "12C",
               file.nm, species.nm, latex = FALSE)


tb.test.12C <- read_test(sv1, block = 4)[[1]] %>%
# filter only 12C
                filter(`Mass#` == "1") %>%
                rename(species.nm = `Mass#`) %>%
                mutate(species.nm = "12C")

tb.comb.12C <- left_join(tb.Xt %>%
                         filter(species.nm == "12C") %>%
# round number as provided data is not accurate enough
                          mutate(Ntot_Xt.pr = round(Ntot_Xt.pr, -2)) %>%
# delete data with insufficient counts
                          filter(Ntot_Xt.pr > 10000),
                       tb.test.12C,
                       by = c("file.nm", "species.nm"))

# and loading the test file
tb.test.R <- read_test(sv1, block = 4)[[2]] %>%
# filter only 13C/12C
               filter(`Ratio#` == "2/1")
# combine datasets
tb.comb.R <- left_join(tb.R %>%
# round number as provided data is not accurate enough
                         mutate(M_R_Xt.pr = round(M_R_Xt.pr, 7)),
                        tb.test.R, by = "file.nm")

test_that("cumulate counts and isotope ratios are similar for EM data", {
  expect_equal(tb.comb.R$M_R_Xt.pr, tb.comb.R$Ratios)
  expect_equal(tb.comb.12C$Ntot_Xt.pr, tb.comb.12C$`Cumulated count`)
})
