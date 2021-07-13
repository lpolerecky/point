#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("consistency of diagnostics for single ion trends", {
  # remove zero count analysis
  tb_0 <- zeroCt(real_IC, "12C", "40Ca 16O", sample.nm, file.nm, .warn = FALSE)
  # predict ionization trends
  expect_snapshot(predict_ionize(tb_0, file.nm, .plot = FALSE))
})
