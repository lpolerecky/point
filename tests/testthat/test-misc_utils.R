#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("the output of zeroCt and cov_R is consistent", {
  # Remove analyses with zero counts
  expect_snapshot(
    zeroCt(real_IC, "13C", "12C", file.nm)
  )
  # wide format ion data frame
  expect_snapshot(
    cov_R(real_IC, c("13C", "12C"), file.nm)
  )
})


test_that("complex elemental configruations can be parsed", {
  expect_equal(
    ion_labeller("12C2-40Ca","webtex"),
    "${}^{12}\\mathrm{C}_{2}$${}^{40}\\mathrm{Ca}_{}$"
  )
  expect_equal(
    ion_labeller("12C2", "expr"),
    substitute("" ^ a * b[c] , env = list(a = "12", b = "C", c = "2"))
  )
})

test_that("convert to wide format ion count data works", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))

  # Processing raw ion count data
  tb_pr <- cor_IC(tb_rw)

  # wide format
  expect_snapshot(
    cov_R(tb_pr, c("13C", "12C"), file.nm)
  )
})
