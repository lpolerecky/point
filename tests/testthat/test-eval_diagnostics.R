#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("consistency of the evaluation of diagnostics on synthetic data", {
  # Simulated IC data
  tb_dia <- diag_R(simu_IC, "13C", "12C", type.nm, spot.nm,
                   .output = "diagnostic")
  # Evaluate significance and effect of outliers based on Cook's D
  expect_snapshot(eval_diag(tb_dia, "13C", "12C", type.nm, spot.nm, .nest = type.nm))
})
