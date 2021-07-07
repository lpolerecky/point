#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("consistency of systematic corrections on internal dataset with custom settings", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"), meta = FALSE)
  # Processing raw ion count data defaults to metadata values
  expect_snapshot(cor_IC(tb_rw))
  # User supplied values
  expect_snapshot(cor_IC(tb_rw, .bl_t = 200))
  expect_snapshot(cor_IC(tb_rw, .bl_t = 200, .deadtime = 44, .det = "EM"))
  expect_snapshot(cor_IC(tb_rw, .bl_t = 200, .thr_PHD = 150, .M_PHD = 210,
                         .SD_PHD = 60, .det = "EM"))
})

test_that("consistency of systematic corrections on internal dataset with internal settings", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"), meta = TRUE)
  # Processing raw ion count data defaults to metadata values
  expect_snapshot(cor_IC(tb_rw))
  # User supplied values
  expect_snapshot(cor_IC(tb_rw, .bl_t = 200))
  expect_snapshot(cor_IC(tb_rw, .bl_t = 200, .deadtime = 44))
  expect_snapshot(cor_IC(tb_rw, .bl_t = 200, .thr_PHD = 150))
  # change name PHD parms see if fhunction change the colnames
  # tb_rw2 <- rename(unfold(tb_rw), mean_PHD.mt = "M_PHD.mt")
  # col_nms <- colnames(tb_rw2)
  # expect_identical(
  #   select(cor_IC(fold(tb_rw2, ".mt"), .bl_t = 200, .thr_PHD = 150, .M_PHD = mean_PHD.mt,
  #                 .hide = FALSE), -c("Y.mt", "dt.rw", "Xt.pr", "N.pr")) %>%
  #     colnames(),
  #   col_nms
  #   )
})

#-------------------------------------------------------------------------------
# errors
#-------------------------------------------------------------------------------
# try to trick the function by having a similar variable in the global env (evaluation ambiguity)
test_that("global environment object with similar name and misc. errors", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"), meta = FALSE)
  N.rw <- 1
  # Single ion internal precision
  expect_equal(cor_IC(tb_rw, .N = N.rw), cor_IC(tb_rw))
  # Not a tibble
  expect_error(cor_IC(data.frame(x = 1:10)), "Ion count dataset should be a tibble object.")
  # supply detector type
  expect_error(cor_IC(tb_rw, .bl_t = 200, .deadtime = 44), "Supply a variable or character string for the detector type to `.det`")
  # supply all variables for correction
  expect_error(cor_IC(tb_rw, .bl_t = 200, .thr_PHD = 150, .M_PHD = 210), "PHD correction requires .M_PHD .SD_PHD")
})

#-------------------------------------------------------------------------------
# default argument injection
#-------------------------------------------------------------------------------

test_that("Inject arguments for ion count variables if not explictely given", {
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
