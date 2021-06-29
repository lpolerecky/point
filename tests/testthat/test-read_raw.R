#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------
# Test datasets in extdata directory
path_point <- point_example("2018-01-19-GLENDON")
test_that("Check change over time in read_IC", {
  expect_snapshot(read_IC(path_point))
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

#-------------------------------------------------------------------------------
# benchmark deprecated
#-------------------------------------------------------------------------------

directory <- point_example("2018-01-19-GLENDON")
read_IC_deprecated <- function(directory, hide = TRUE){

  # List files
  ls_IC <- read_validator(directory, "is_txt")[["ion"]]
  n_max <- Inf

  # Collecting measurement data
  purrr::map2_dfr(
    ls_IC,
    n_max,
    ~readr::read_tsv(
      .x,
      col_names = c("t.nm", "N.rw"),
      col_types = "-cc",
      comment = "B",
      skip = 1,
      # n-max is n times number of species
      n_max = .y
    ),
    .id = "file.nm"
  ) %>%
    # Remove old column headers
    filter(.data$t.nm != "X", .data$N.rw  != "Y") %>%
    # Coercion to numeric values
    mutate(
      t.nm = as.numeric(.data$t.nm),
      N.rw = as.numeric(.data$N.rw)
    )
}

microbenchmark::microbenchmark(
  point::read_IC(point_example("2018-01-19-GLENDON"), meta = FALSE),
  read_IC_deprecated(point_example("2018-01-19-GLENDON")),
  times = 20
)
