#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("consistency of residual diagnostics on internal dataset", {
  # Descriptive an predictive statistics for 13C/12C ratios (note .output
  # argument and remove zero count analysis)
  tb_R <- stat_R(real_IC, "13C", "12C", file.nm, .output = "complete",
                 .zero = TRUE)
  # residual based augmentation of ion count data for isotope ratios
  expect_snapshot(CooksD(tb_R, "13C", "12C", file.nm, .X = Xt.pr, .N = N.pr,
                         .species = species.nm, .t = t.nm, .output = "flag"))
  expect_snapshot(Rm(tb_R, "13C", "12C", file.nm, .X = Xt.pr, .N = N.pr,
                     .species = species.nm, .t = t.nm, .output = "flag"))
  expect_snapshot(CV(tb_R, "13C", "12C", file.nm, .X = Xt.pr, .N = N.pr,
                     .species = species.nm, .t = t.nm, .output = "flag"))
  expect_snapshot(QQ(tb_R, "13C", "12C", file.nm, .X = Xt.pr, .N = N.pr,
                     .species = species.nm, .t = t.nm, .output = "flag"))
  expect_snapshot(norm_E(tb_R, "13C", "12C", file.nm, .X = Xt.pr, .N = N.pr,
                         .species = species.nm, .t = t.nm, .output = "flag"))
  expect_snapshot(IR(tb_R, "13C", "12C", file.nm, .X = Xt.pr, .N = N.pr,
                     .species = species.nm, .t = t.nm, .output = "flag"))
})
