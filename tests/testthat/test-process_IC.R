#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("consistency of systematic corrections on internal dataset", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
  # Processing raw ion count data defaults to metadata values
  expect_snapshot(cor_IC(tb_rw))
  # User supplied values
  expect_snapshot(cor_IC(tb_rw, .bl_t = 0.2))
  expect_snapshot(cor_IC(tb_rw, .bl_t = 0.2, .deadtime = 44))
  expect_snapshot(cor_IC(tb_rw, .bl_t = 0.2, .thr_PHD = 150))
})

#-------------------------------------------------------------------------------
# errors
#-------------------------------------------------------------------------------
# try to trick the function by having a similar variable in the global env (evaluation ambiguity)
test_that("global environment object with similar name and misc. errors", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
  N.rw <- 1
  # Single ion internal precision
  expect_equal(cor_IC(tb_rw, file.nm, .N = N.rw), cor_IC(tb_rw, file.nm))
  # Error for unknown variables
  expect_error(cor_IC(tb_rw, file.nm, .N = x), "Try explicitly supplying the variables.")
  expect_error(cor_IC(tb_rw, file.nm, .bl_t = x), "Try explicitly supplying the variables.")
  expect_error(cor_IC(tb_rw, file.nm, .bl_t = "f"), "Try explicitly supplying the variables.")
  expect_error(cor_IC(tb_rw, file.nm, .det = 1), "Try explicitly supplying the variables.")
  # Not a tibble
  expect_error(cor_IC(data.frame(x = 1:10)), "Ion count dataset should be a tibble object.")
})

#-------------------------------------------------------------------------------
# default argument injection
#-------------------------------------------------------------------------------

test_that("Inject arguments for ion count variables if not explictely given", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
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
