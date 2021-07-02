#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------
# Test datasets in extdata directory
path_point <- point_example("2018-01-19-GLENDON")
test_that("Check change over time in read_IC", {
  expect_snapshot(read_IC(path_point))
  expect_snapshot(read_meta(path_point))
})

#-------------------------------------------------------------------------------
# Helper functions
#-------------------------------------------------------------------------------

# Error when directory argument class is wrong
test_that("read_IC input error", {
  expect_error(read_validator(tempdir()), "`directory` does not contain required filetypes: .is_txt, .chk_is, and .stat")
  expect_warning(unfold(read_IC(path_point), type = "other"), "Attribute unavailable.")
  })

# Testing the class of the retrieved data
test_that("read_IC creates a tibble", {
  expect_s3_class(read_IC(path_point), "tbl_df")
})

# Directory check
test_that("directory check", {
  expect_snapshot(ICdir_chk(point_example("2018-01-19-GLENDON")))
})

