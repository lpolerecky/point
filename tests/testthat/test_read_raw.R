context("Reading Ion Count Data")
library(point)

# Test datasets in extdata directory
sv1 <- point_example("2018-01-19-GLENDON")

# Processing raw ion count data
tb_pr <- suppressWarnings(cor_IC(read_IC(sv1)))

#-------------------------------------------------------------------------------
# Variable class checks and compatibility of txt and stat files
#-------------------------------------------------------------------------------
# Error when directory argument class is wrong
testthat::test_that("read_IC input error", {
  expect_error(read_IC(1))
  expect_error(read_IC(1L))
  })


# # Warning when zero count exist
# testthat::test_that("zero counts", {
#   expect_warning(
#     stat_R(tb_pr, "13C", "13C 14N", file.nm, .zero = TRUE))
#   })

# Testing the class of the retrieved data
testthat::test_that("read_IC creates a tibble", {
  expect_is(suppressWarnings(read_IC(sv1)), "tbl_df")
  })
