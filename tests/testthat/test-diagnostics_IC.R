#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------
test_that("consistency of diagnostics wrapper on synthetic data", {
  # residual based augmentation of ion count data for isotope ratios
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm)
  )
  # change method
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .method = "QQ")
  )
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .method = "IR")
  )
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .method = "CV")
  )
  # change output (inference is default)
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .output = "complete")
  )
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .output = "augmented")
  )
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .output = "diagnostic")
  )
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .output = "outlier")
  )
  # latex labs
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .label = "latex")
  )
  # with nesting and latex labs
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .nest = type.nm,
           .label = "latex")
  )
})


#-------------------------------------------------------------------------------
# Is metadata preserved
#-------------------------------------------------------------------------------

test_that("Keep metadata", {
  expect_snapshot(diag_R(real_IC, "13C", "12C", file.nm))
  expect_snapshot(diag_R(real_IC, "13C", "12C", file.nm, .meta = TRUE)
                  %>% unfold())
})

#-------------------------------------------------------------------------------
# errors
#-------------------------------------------------------------------------------

test_that("errors in diag_R call", {
  # Missing variables
  expect_error(
    diag_R(real_IC, "13C", "12C", file.nm, .X = L),
    "Tibble does not contain the supplied variables!"
  )
  expect_error(
    diag_R(select(real_IC, -N.pr), "13C", "12C", file.nm),
    "Tibble does not contain the default variables!"
  )
})
