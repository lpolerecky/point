#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("consistency of plotting", {
  # only use for development purposes
  skip_on_cran()
  skip_on_ci()
  skip_on_covr()

  # complete diagnostics
  tb_dia <- diag_R(simu_IC, "13C", "12C", type.nm, spot.nm,
                   .output = "diagnostics")
  plot_args <- list(.method = "CooksD", .alpha_level = 0.05,
                    .plot_type = "static", .plot_stat = NULL)
  p <- print(gg_IC(tb_dia, "13C", "12C", type.nm, spot.nm, .X = Xt.pr, .N = N.pr,
                   .flag = NULL, .plot_args = plot_args))
  # ggplot plotting
  expect_true(ggplot2::is.ggplot(p))
  # snapshot
  vdiffr::expect_doppelganger("ggplot2 with two-D density", p)
})
