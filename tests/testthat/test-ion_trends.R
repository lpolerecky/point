#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("consistency of diagnostics for single ion trends", {
  # raw ion counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
  # Processing raw ion count data
  tb_pr <- cor_IC(tb_rw)
  # remove zero count analysis
  tb_0 <- zeroCt(tb_pr, "12C", "40Ca 16O", sample.nm, file.nm, .warn = FALSE)
  # predict ionization trends
  expect_snapshot(predict_ionize(tb_0, file.nm, .plot = FALSE))
})
