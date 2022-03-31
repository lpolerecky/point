#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("simple QSA test is consistent", {
  # QSA test
  expect_snapshot(
    QSA_test(real_IC, "13C", "12C", file.nm)
  )
})

test_that("grouped QSA test is consistent", {

  skip_if_not_installed("broom.mixed")

  # grouped
  expect_snapshot(
    QSA_test(real_IC, "13C", "12C", sample.nm, file.nm, .nest = file.nm)
  )
})
