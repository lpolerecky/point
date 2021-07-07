#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------
test_that("consistency of diagnostics wrapper on synthetic data", {
  # residual based augmentation of ion count data for isotope ratios
  expect_snapshot(diag_R(simu_IC, "13C", "12C", type.nm, spot.nm))
})


#-------------------------------------------------------------------------------
# Is metadata preserved
#-------------------------------------------------------------------------------

testthat("Keep metadata", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
  # Processing raw ion count data
  tb_pr <- cor_IC(tb_rw)
  diag_R(tb_pr, "13C", "12C", sample.nm.nm, file.nm)
})
