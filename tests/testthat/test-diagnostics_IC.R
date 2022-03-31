#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------
test_that("diagnostics wrapper on synthetic data is consistent", {
  # residual based augmentation of ion count data for isotope ratios
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm)
  )
  # change method
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

test_that("QQ diagnostic on synthetic data is consistent", {
  skip_if_not_installed("nortest")
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .method = "QQ")
  )
})

test_that("IR diagnostic on synthetic data is consistent", {
  skip_if_not_installed("stats")
  expect_snapshot(
    diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .method = "IR")
  )
})

#-------------------------------------------------------------------------------
# specific error while using raster data in pointapply
#-------------------------------------------------------------------------------

test_that("case-study with Utrecht raster dataset", {

  skip_if_offline()
  skip_on_cran()
  skip_on_covr()
  skip_on_ci()
  skip_if_not_installed("pointapply")
  skip_if_not_installed("readmat")

  library(pointapply)

  # download data from ZENODO repo
  # download_point(type = "raw")
  # MEX_files <- list.files(get_matlab("2020-08-20-GLENDON"), full.names = TRUE)
  # MEX <- purrr::map(MEX_files, ~readmat::read_mat(.x))
  # grid_aggregate(MEX, c("height", "width", "depth"), grid_cell = 64,
                 # species = c("12C", "13C"), title = "MEX",
                 # name = "map_sum_grid", corrected = TRUE, save = TRUE)
  load_point("map_sum_grid", "MEX", 64)

  # diagnostics
  expect_snapshot(
    diag_R(
      map_sum_grid_64_MEX,
      "13C",
      "12C",
      dim_name.nm,
      sample.nm,
      file.nm,
      grid.nm,
      .nest = grid.nm
    )
  )

})

#-------------------------------------------------------------------------------
# Is metadata preserved
#-------------------------------------------------------------------------------

test_that("Keep metadata", {
  expect_snapshot(
    diag_R(real_IC, "13C", "12C", file.nm)
  )
  expect_snapshot(
    diag_R(real_IC, "13C", "12C", file.nm, .meta = TRUE) |> unfold()
  )
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
    diag_R(dplyr::select(real_IC, -N.pr), "13C", "12C", file.nm),
    "Tibble does not contain the default variables!"
  )
})
