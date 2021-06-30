#-------------------------------------------------------------------------------
# minimum functionality (over time with snapshots)
#-------------------------------------------------------------------------------

test_that("consistency of precision estimates on internal dataset", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
  # Processing raw ion count data
  tb_pr <- cor_IC(tb_rw)
  # Single ion internal precision
  expect_snapshot(stat_X(tb_pr, file.nm))
  expect_snapshot(stat_X(tb_pr, file.nm, .label = "latex"))
  # Isotope ratio internal precision
  expect_snapshot(stat_R(tb_pr, "13C", "12C", file.nm, .zero = TRUE))
  expect_snapshot(stat_R(tb_pr, "13C", "12C", file.nm, .zero = TRUE,
                         .label = "latex"))
  # Isotope ratio external precision
  expect_snapshot(stat_R(tb_pr, "13C", "12C", sample.nm, file.nm,
                         .nest = file.nm, .zero = TRUE))
  # Isotope ratio external precision and statistic selection
  expect_snapshot(stat_R(tb_pr, "13C", "12C", sample.nm, file.nm,
                         .nest = file.nm, .zero = TRUE, .label = "latex",
                         .stat = c("M", "RS")))
})

#-------------------------------------------------------------------------------
# errors
#-------------------------------------------------------------------------------
# try to trick the function by having a similar variable in the global env (evaluation ambiguity)
test_that("global environment object with similar name and misc. errors", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
  # Processing raw ion count data
  tb_pr <- cor_IC(tb_rw)
  # Global environment variable similar to data env variable
  Xt.pr <- 1
  # Single ion internal precision
  expect_equal(stat_X(tb_pr, file.nm, .X = Xt.pr), stat_X(tb_pr, file.nm))
  # Error for unknown variables
  expect_error(stat_X(tb_pr, file.nm, .X = x), "Try explicitly supplying the variables.")
  expect_error(stat_X(tb_pr, "13C", "12C", file.nm, .X = x), "Try explicitly supplying the variables.")
  # Error for unknown statistic
  expect_error(stat_X(tb_pr, file.nm, .stat = "x"), "Unkown statistic.")
  expect_error(stat_X(tb_pr, "13C", "12C", file.nm, .stat = "x"), "Unkown statistic.")
  expect_error(stat_R(tb_pr, "13C", "12C", file.nm, .zero = TRUE,
                      .label = "latex", .output = "complete"),
               "Latex labels is not supported for complete datasets.")
  })

#-------------------------------------------------------------------------------
# Compare to Cameca output
#-------------------------------------------------------------------------------

# single ion precision
test_that("consistency of single ion precision estimates compared to Cameca output", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
  # Processing raw ion count data
  tb_pr <- cor_IC(tb_rw)
  # stat X
  tb_X <- stat_X(tb_pr, file.nm, .stat = "tot", .meta = TRUE)
  tb_meta <- unfold(tb_X, merge = FALSE) %>%
    distinct(species.nm, det.mt) %>%
    mutate(det.mt = readr::parse_number(det.mt))
  tb_X <- left_join(tb_X, tb_meta, by = "species.nm") %>% arrange(file.nm, det.mt)
  # Cameca reference
  ref_X <- filter(cameca_stat_X, `correction block` == 1)$`Cumulated count`
  # tolerance lowered as cameca stat file only contains six signficant digits
  expect_equal(tb_X$tot_N.pr, ref_X, tolerance = 10)
})

# R precision
test_that("consistency of isotope ratio precision estimates compared to Cameca output", {
  # raw data containing 13C and 12C counts on carbonate
  tb_rw <- read_IC(point_example("2018-01-19-GLENDON"))
  # Processing raw ion count data
  tb_pr <- cor_IC(tb_rw)
  # stat X
  tb_R <- stat_R(tb_pr, "13C", "12C", file.nm,
                 .stat = c("M", "RSeM", "hat_RSeM"), .zero = TRUE)
  # Cameca reference
  ref_R <- filter(cameca_stat_R, `correction block` == 1, `Ratio#` == "2/1")
  # tolerance lowered as Cameca stat file only contains six, three and three
  # significant digits
  expect_equal(tb_R$M_R_Xt.pr, ref_R$Ratios, tolerance = 1e-5)
  expect_equal(tb_R$RSeM_R_Xt.pr, ref_R$`Err_mean(%)` * 10, tolerance = 1e-2)
  expect_equal(tb_R$hat_RSeM_R_N.pr, ref_R$`Err_Poisson(%)` * 10, tolerance = 1e-2)
})

#-------------------------------------------------------------------------------
# helpers
#-------------------------------------------------------------------------------
tb_tex <- mutate(
  point::names_stat_R,
  origin = case_when(origin == "X" ~ "M", origin == "N" ~ "Ntot")
  )



test_that("Latex printing of variables for tables", {
  # Latex printing of variables
  # single stat selection (internal precision)
  expect_equal(tex_labeller(point::names_stat_R, "M", "latex"), "$\\bar{R}$")
  # single stat selection (external precision)
  expect_equal(tex_labeller(tb_tex , "M", "latex"), "$\\bar{\\bar{R}}$")
  # multiple  stat selection (external precision)
  expect_equal(
    tex_labeller(tb_tex , c("M", "RS"), "latex"),
    c("$\\bar{\\bar{R}}$", "$\\epsilon_{\\bar{R}}$ (\\text{\\textperthousand})")
    )
})




test_that("Argument selector based on stat selection", {
  ls_nm <- c(n_R  = "n_R_t.nm", M_R = "M_R_Xt.pr",  RS_R = "RS_R_Xt.pr", RSeM_R ="RSeM_R_Xt.pr")
  stat <- "RS"
  # argument selector based on stat selection
  expect_equal(stat_selector(stat, ls_nm)$pos, c(RS_R = "RS_R_Xt.pr"))
  expect_equal(stat_selector(stat, ls_nm)$neg, c(n_R  = "n_R_t.nm", M_R = "M_R_Xt.pr", RSeM_R ="RSeM_R_Xt.pr"))
  })

