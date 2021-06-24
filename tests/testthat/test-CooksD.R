#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("consistency of residual diagnostics on internal dataset", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
  # Processing raw ion count data
  tb_pr <- cor_IC(tb_rw)
  # Descriptive an predictive statistics for 13C/12C ratios (note .output
  # argument and remove zero count analysis)
  tb_R <- stat_R(tb_pr, "13C", "12C", file.nm, sample.nm, .output = "complete",
                 .zero = TRUE)
  # residual based augmentation of ion count data for isotope ratios
  expect_snapshot(CooksD(tb_R, "13C", "12C", file.nm, .output = "flag"))
  expect_snapshot(Rm(tb_R, "13C", "12C", file.nm, .output = "flag"))
  expect_snapshot(CV(tb_R, "13C", "12C", file.nm, .output = "flag"))
  expect_snapshot(QQ(tb_R, "13C", "12C", file.nm, .output = "flag"))
  expect_snapshot(norm_E(tb_R, "13C", "12C", file.nm, .output = "flag"))
  expect_snapshot(IR(tb_R, "13C", "12C", file.nm, .output = "flag"))
})
