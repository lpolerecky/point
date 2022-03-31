test_that("isotope values can be converted", {
  expect_snapshot(
    calib_R(0.0111, reference = 0.011237, type = "composition", input = "R",
            output = "delta")
  )
  expect_snapshot(
    calib_R(0.0111, reference = "VPDB", isotope = "13C",
            type = "composition", input = "R", output = "delta")
  )
  expect_snapshot(
    calib_R(0.0111, reference = 0.011237, type = "composition",
            input = "R", output = "F")
  )
  expect_snapshot(
    calib_R(-25, reference = "VPDB", isotope = "13C",
            type = "enrichment", input = "delta", output = "alpha", y = -105)
  )
})
