context("Reading Ion Count Data")
library(point)

# test datasets in extdata directory
sv1 <- system.file("extdata", "2018-01-20-GLENDON", package = "point")
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
  expect_warning(read_validator(sv1
  ))
})

# warning when zero count exist
test_that("zero counts", {
  expect_warning(
    zeroCt(
      read_IC(sv1), N.rw, species.nm, "13C", "12C", file.nm)
    )
  }
  )

# testing the class of the retrieved data
test_that("read_IC creates a tibble", {
  expect_is(read_IC(sv2), "tbl_df")
            }
  )
