# context("Cross-validation of Ion Count statistics with Cameca output")
# library(point)
#
# # Test datasets in extdata directory
# sv1 <- point_example("2018-01-19-GLENDON")
#
# #-------------------------------------------------------------------------------
# # Value checks of raw data output. Cross validation of statistics results with
# # the processed values as given by the Cameca softwarre. For this, the function
# # read_test is used, whiche extracts the Cameca processed data.
# #-------------------------------------------------------------------------------
#
# # function to read CAMECA output to validate point output
# read_test <- function(directory, block = 1){
#
#   # block for checking
#   if (block == 1) {bl <- 2}
#   if (block == 4) {bl <- max}
#
#   # extract stat files with diagnostics of the machine and statistics
#   l.s <- dir(directory,
#              pattern = ".stat$") %>%
#     # set names for subsequent storage
#     purrr::set_names() %>%
#     # remove transect files
#     purrr::discard(., stringr::str_detect(., "transect.stat"))
#
#   max_n <- lapply(
#     purrr::map(
#       l.s,
#       ~readr::read_lines(paste0(directory, "/", .), n_max = -1)
#       ),
#     stringr::str_which, "Mass#") %>%
#     purrr::map(., bl) %>%
#     purrr::flatten_dbl()
#
#   # Test dataset
#   l.R <- lst(a = l.s, b =  max_n  + 10 , c = 5) # isotope
#   l.N <- lst(a = l.s, b =  max_n, c = 7) # cumulative count
#
#   # Function for R
#   fun_read.R <- function(a, b, c) {
#
#     readr::read_table(
#       paste0(directory, "/", a),
#       skip =  b,
#       n_max = c,
#       col_names = c("R.nm", "M_R_X.test", "hat_RSeM_R_N.test", "RSeM_R_X.test",
#                     "chi2_R_N.test"),
#       col_types = "cdddd-"
#       ) %>%
#       mutate(
#         RSeM_R_X.test = RSeM_R_X.test * 10,
#         hat_RSeM_R_N.test = hat_RSeM_R_N.test * 10
#         )
#   }
#
#   # Function for N
#   fun_read.N <- function(a, b, c) {
#
#     readr::read_table(
#       paste0(directory, "/", a),
#       skip = b,
#       n_max = c ,
#       col_names = c("num.mt", "Ntot_X.test"),
#       col_types = "cd"
#       ) %>%
#       tidyr::drop_na()
#   }
#
#   tb.test.R <- purrr::pmap_dfr(l.R, fun_read.R, .id = "file.nm") %>%
#     mutate(
#       file.nm =
#         stringr::str_sub(
#           file.nm,
#           end = (stringr::str_length(file.nm) - 5)
#           )
#       )
#
#   tb.test.N <- purrr::pmap_dfr(l.N, fun_read.N, .id = "file.nm") %>%
#     mutate(
#       file.nm =
#         stringr::str_sub(
#           file.nm,
#           end = (stringr::str_length(file.nm) - 5)
#           )
#       )
#
#   return(list(tb.test.N = tb.test.N, tb.test.R = tb.test.R))
# }
#
#
#
# # function for block specific comparison of Cameca and here procesed data
# # block 1 is uncorrected data
# # block 4 is yield and deadtime corrected EM data and background corrected FC data
# validate_IC <- function(ion1, ion2, block){
#
# # block additional parameters for correcting data
#   if (block == 1) {
#     deadtime <- 0
#     thr <- 0
#     }
#   if (block == 4) {
#     deadtime <- 44
#     thr <- 50
#     }
#
# # Perform the routine point workflow for IC data
#   tb_rw <- suppressWarnings(read_IC(sv1))
#
# # Processing raw ion count data
#   tb_pr <- cor_IC(tb_rw)
#
# # Remove zero count analysis
#   tb_rw <- zeroCt(tb_pr, "13C", "12C", file.nm, .warn = FALSE)
#
# # single ion count statistics
#   tb_X <- stat_X(tb_pr, file.nm)
#
# # Normal descriptive an predictive statistics for 13C/12C ratios (global)
#   tb_R <- stat_R(tb_pr, .ion1 = "13C", .ion2 = "12C", file.nm)
#
# #--------------------------------------
# # Reading CAMECA stat files
# #--------------------------------------
#
# # Ion counts in comparison dataset
#   tb_test_X <- read_test(sv1, block = block)[[1]] %>%
#     mutate(species.nm =
#              recode(
#                num.mt,
#                !!! set_names(
#                  unique(tb_rw$species.nm),
#                  nm = unique(tb_rw$num.mt))
#                 )
#            )
# # Combine datasets
#   tb_test_X <- left_join(tb_X, tb_test_X, by = c("file.nm", "species.nm"))
#
# # Ratios provide in comparison dataset
#   R_vec <- set_names(
#              outer(unique(tb.rw$species.nm),
#                    unique(tb.rw$species.nm),
#                    paste, sep = "/") %>%
#                as.vector(),
#               nm = outer(unique(tb.rw$num.mt),
#                          unique(tb.rw$num.mt),
#                          paste, sep = "/") %>%
#                as.vector())
#
#   R_vec <- R_vec[names(R_vec) %in% read_test(sv1, block = block)[[2]]$R.nm]
#
# # Loading the test file for isotope ratios
#   tb.test.R <- read_test(sv1, block = block)[[2]] %>%
#                  mutate(R.nm = recode(R.nm, !!! R_vec)) %>%
# # Combine datasets
#                  left_join(tb.R, . ,by = c("file.nm", "R.nm"))
#
#   return(list(tb.test.Xt, tb.test.R))
# }
#
# #-------------------------------------------------------------------------------
# # Validation test for 13C and 12C counts and 13C/12C ratios
# #-------------------------------------------------------------------------------
#
# # Ions
# ion1 <- "13C"
# ion2 <- "12C"
#
#
# #--------------------------------------
# # Unit 1:
# #--------------------------------------
#
# # Call validation function for block 1
# val_IC <- suppressWarnings(validate_IC(ion1 = ion1, ion2 = ion2, block = 1))
#
# # Total ion counts
# Xt1 <- val_IC[[1]] %>%
#   filter(species.nm == ion1 | species.nm == ion2)
#
# # Isotope ratios
# R1 <- val_IC[[2]] %>%
#   filter(R.nm == paste(ion1, ion2, sep = "/"))
#
# test_that("Block 1:", {
#   expect_equal(Xt1 %>%
#                  filter(species.nm == "12C") %>%
#                  pull(Ntot_Xt.pr),
#                Xt1 %>%
#                  filter(species.nm == "12C") %>%
#                  pull(Ntot_Xt.test),
#                tolerance = 1e-6, # relative difference
#                label = "Ntot of 12C",
#                expected.label = "Ntot of 12C in the Cameca stat-file")
#   expect_equal(Xt1 %>%
#                  filter(species.nm == "13C") %>%
#                  pull(Ntot_Xt.pr),
#                Xt1 %>%
#                  filter(species.nm == "13C") %>%
#                  pull(Ntot_Xt.test),
#                tolerance = 1e-6, # relative difference
#                label = "Ntot of 13C",
#                expected.label = "Ntot of 13C in the Cameca stat-file")
#   expect_equal(R1 %>%
#                  pull(M_R_Xt.pr),
#                R1 %>%
#                  pull(M_R_Xt.test),
#                tolerance = 1e-7, # absolute difference
#                label = "Mean of 13C/12C",
#                expected.label = "mean of 13C/12C in the Cameca stat-file")
#   expect_equal(R1 %>%
#                  pull(RSeM_R_Xt.pr),
#                R1 %>%
#                  pull(RSeM_R_Xt.test),
#                tolerance = 1e-1, # relative difference within a percent
#                label = "Relative standard error of 13C/12C",
#                expected.label = c("relative standard error of  13C/12C in the
#                                   Cameca stat-file"))
#   expect_equal(R1 %>%
#                  pull(hat_RSeM_R_Xt.pr),
#                R1 %>%
#                  pull(hat_RSeM_R_Xt.test),
#                tolerance = 1e-2, # relative difference within a percent
#                label = "Predicted relative standard error of 13C/12C",
#                expected.label = c("predicted relative standard error of  13C/12C
#                                   in the Cameca stat-file"))
# })
#
