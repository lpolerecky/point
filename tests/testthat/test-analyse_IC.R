#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("consistency of precision estimates on internal dataset", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
  # Processing raw ion count data
  tb_pr <- cor_IC(tb_rw)
  # Single ion internal precision
  expect_snapshot(stat_X(tb_pr, file.nm))
  expect_snapshot(stat_X(tb_pr, file.nm, .label = "latex"))
  # Isotope ratio internal precision
  expect_snapshot(stat_R(tb_pr, "13C", "12C", file.nm, .zero = TRUE))
  expect_snapshot(stat_R(tb_pr, "13C", "12C", file.nm, .zero = TRUE,
                         .label = "latex"))
  # Isotope ratio external precision
  expect_snapshot(stat_R(tb_pr, "13C", "12C", sample.nm, file.nm,
                         .nest = file.nm, .zero = TRUE))
})

#-------------------------------------------------------------------------------
# errors
#-------------------------------------------------------------------------------
# try to trick the function by having a similar variable in the global env (evaluation ambiguity)
test_that("global environment object with similar name and misc. errors", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
  # Processing raw ion count data
  tb_pr <- cor_IC(tb_rw)
  # Global environment variable similar to data env variable
  Xt.pr <- 1
  # Single ion internal precision
  expect_equal(stat_X(tb_pr, file.nm, .X = Xt.pr), stat_X(tb_pr, file.nm))
  # Error for unknown variables
  expect_error(stat_X(tb_pr, file.nm, .X = x), "Try explicitly supplying the variables.")
  expect_error(stat_X(tb_pr, "13C", "12C", file.nm, .X = x), "Try explicitly supplying the variables.")
  # Error for unknown statistic
  expect_error(stat_X(tb_pr, file.nm, .stat = "x"), "Unkown statistic.")
  expect_error(stat_X(tb_pr, "13C", "12C", file.nm, .stat = "x"), "Unkown statistic.")
  expect_error(stat_R(tb_pr, "13C", "12C", file.nm, .zero = TRUE,
                      .label = "latex", .output = "complete"),
               "Latex labels is not supported for complete datasets.")
  })

