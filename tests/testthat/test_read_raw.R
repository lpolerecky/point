context("Reading Ion Count Data")
library(point)

# Test datasets in extdata directory
sv1 <- point_example("2018-01-19-GLENDON")
# sv2 <- point_example("2018-01-20-GLENDON")

# Processing raw ion count data
tb.pr <- suppressWarnings(cor_IC(read_IC(sv1), N.rw, t.rw, det_type.mt,
                                 deadtime = 44, thr_PHD = 50))

#-------------------------------------------------------------------------------
# Variable class checks and compatibility of txt and stat files
#-------------------------------------------------------------------------------
# Error when directory argument class is wrong
test_that("read_IC input error", {
  expect_error(read_IC(1))
  expect_error(read_IC(1L))
  })

# # Warning when empty text files exist
# test_that("empty txt file", {
#   expect_warning(read_IC(sv2))
#   })

# Warning when zero count exist
test_that("zero counts", {
  expect_warning(
    stat_R(tb.pr, Xt.pr, N.pr, species.nm, ion1 = "13C",
           ion2 = "12C", file.nm, zero = TRUE))
  })

# Testing the class of the retrieved data
test_that("read_IC creates a tibble", {
  expect_is(read_IC(sv1), "tbl_df")
  })
