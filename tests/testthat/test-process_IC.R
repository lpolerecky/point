#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------
test_that(paste0("systematic corrections on internal dataset ",
                 "with custom settings are consistent"), {

  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"), meta = FALSE)

  # Processing raw ion count data defaults to metadata values
  expect_snapshot(
    cor_IC(tb_rw)
  )
  # User supplied values
  expect_snapshot(
    cor_IC(tb_rw, .bl_t = 200)
  )
  expect_snapshot(
    cor_IC(tb_rw, .bl_t = 200, .deadtime = 44, .det = "EM")
  )
})

test_that(paste0("systematic PHD corrections on internal dataset ",
                 "with custom settings are consistent"), {

  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"), meta = FALSE)

  skip_if_not_installed("polyaAeppli")

  expect_snapshot(
    cor_IC(tb_rw, .bl_t = 200, .thr_PHD = 150, .M_PHD = 210,
           .SD_PHD = 60, .det = "EM")
  )
})

test_that(paste0("systematic corrections on internal dataset ",
                 "with internal settings are consistent"), {

  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"), meta = TRUE)

  # Processing raw ion count data defaults to metadata values
  expect_snapshot(
    cor_IC(tb_rw)
  )
  # User supplied values
  expect_snapshot(
    cor_IC(tb_rw, .bl_t = 200)
  )
  expect_snapshot(
    cor_IC(tb_rw, .bl_t = 200, .deadtime = 44)
  )
})

test_that(paste0("systematic PHD corrections on internal dataset ",
                   "with internal settings are consistent"), {

  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"), meta = FALSE)

  skip_if_not_installed("polyaAeppli")

  expect_snapshot(
    cor_IC(tb_rw, .bl_t = 200, .thr_PHD = 150)
  )
})

#-------------------------------------------------------------------------------
# errors
#-------------------------------------------------------------------------------
# try to trick the function by having a similar variable in the global env
# (evaluation ambiguity)
test_that("global environment object with similar name and misc. errors", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"), meta = FALSE)
  N.rw <- 1
  # Single ion internal precision
  expect_equal(
    cor_IC(tb_rw, .N = N.rw),
    cor_IC(tb_rw)
  )
  # Not a tibble
  expect_error(
    cor_IC(data.frame(x = 1:10)),
    "Ion count dataset should be a tibble object."
  )
  # Missing variables
  expect_error(
    cor_IC(dplyr::select(tb_rw, - N.rw)),
    "Tibble does not contain the default variables!"
  )
  expect_error(
    cor_IC(tb_rw, .N = Xt),
    "Tibble does not contain the supplied variables!"
  )
  # supply detector type
  expect_error(
    cor_IC(tb_rw, .bl_t = 200, .deadtime = 44),
    "Supply a variable or character string for the detector type to `.det`"
  )
})

test_that(paste0("global environment object for PHD with similar name"), {

  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"), meta = FALSE)
  N.rw <- 1

  skip_if_not_installed("polyaAeppli")

  # supply all variables for correction
  expect_error(
    cor_IC(tb_rw, .bl_t = 200, .thr_PHD = 150, .M_PHD = 210),
    "PHD correction requires .M_PHD .SD_PHD"
  )
})

#-------------------------------------------------------------------------------
# default argument injection
#-------------------------------------------------------------------------------

test_that(paste0("arguments are injected for ion count variables if not ",
                 "explictely provided"), {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"), meta = TRUE)
  # NULL
  expect_output(
    str(
      inject_args(
        unfold(tb_rw),
        rlang::quos(.N = N.rw, .t = NULL, .bl_t = NULL, .det = NULL,
                    .M_PHD = NULL, .SD_PHD = NULL),
        type = c("raw", "group", "meta")
      )
    ),
    "List of 6"
    )
  # NULL
  expect_output(
    str(
      inject_args(
        unfold(tb_rw),
        rlang::quos(.N = N.rw, .t = NULL, .bl_t = 1, .det = NULL, .M_PHD = NULL,
                    .SD_PHD = NULL),
        type = c("raw", "group", "meta")
      )
    ),
    "List of 6"
  )
})
