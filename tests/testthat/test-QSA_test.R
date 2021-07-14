#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("consistency of the QSA test", {
  # QSA test
  expect_snapshot(QSA_test(real_IC, "13C", "12C", file.nm))
  # grouped
  expect_snapshot(QSA_test(real_IC, "13C", "12C", sample.nm, file.nm, .nest = file.nm))
})
