context("Reading Ion Count Data")
library(point)

# additional test dataset system.file("extdata", "2020-01-17-TREASURE", package = "point")

sv1 <- system.file("extdata", "2018-01-19-GLENDON", package = "point")
sv2 <- system.file("extdata", "2020-01-17-TREASURE", package = "point")

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
