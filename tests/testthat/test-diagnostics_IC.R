#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------
test_that("consistency of diagnostics wrapper on synthetic data", {
  # residual based augmentation of ion count data for isotope ratios
  expect_snapshot(diag_R(simu_IC, "13C", "12C", type.nm, spot.nm))
})
