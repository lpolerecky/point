#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("consistency of the zeroCt and cov_R", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
  # Processing raw ion count data
  tb_pr <- cor_IC(tb_rw)
  # Remove analyses with zero counts
  expect_snapshot(zeroCt(tb_pr, "13C", "12C", file.nm))
  expect_snapshot(cov_R(tb_pr, c("13C", "12C"), file.nm))
})
