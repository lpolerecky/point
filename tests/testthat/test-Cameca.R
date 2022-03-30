#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("Cameca diagnostics on internal dataset are consistent", {
  # Descriptive an predictive statistics for 13C/12C ratios (note .output
  # argument and remove zero count analysis)
  tb_R <- stat_R(real_IC, "13C", "12C", file.nm, sample.nm, .output = "complete",
                 .zero = TRUE)
  # CAMECA style augmentation of ion count data for isotope ratios
  expect_snapshot(Cameca(tb_R, "13C", "12C", file.nm, .X = Xt.pr, .N = N.pr,
                         .species = species.nm, .t = t.nm, .output = "flag"))
})
