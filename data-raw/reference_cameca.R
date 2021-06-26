#-------------------------------------------------------------------------------
# Value checks of raw data output. Cross validation of statistics results with
# the processed values as given by the Cameca software. For this, the function
# read_test is used, which extracts the Cameca processed data.
#-------------------------------------------------------------------------------

# function to read CAMECA output to validate point output
read_cameca <- function(directory, pattern_begin, table_depth,
                        pattern_end = NULL, file_type, col_types = NULL,
                        table_duplicates) {

  ls_files <- ICdir_chk(directory,  file_type)
  min_row <- row_scanner(directory, ls_files[[file_type]], pattern_begin,
                         return_line = TRUE)
  if (!is.null(pattern_end)) {
    max_row <- row_scanner(directory, ls_files[[file_type]], pattern_end)
    table_row <- purrr::map2(max_row, min_row[[1]], ~.x - .y - 3)
  } else if (!is.null(table_depth)) {
  table_row <- purrr::map(min_row[[1]], ~rep(table_depth, length(.x)))
  }
  var_names <- stringr::str_split(min_row[[2]], "\\s(?=[[:upper:]])")[[1]] %>%
    purrr::keep(stringr::str_detect, pattern = "\\b") %>%
    stringr::str_trim()

  # Function for R
  read_cameca_table <- function(cameca_file, skip_rows, table_rows, var_names, col_types) {
    path_stat <- fs::path(directory, cameca_file)
    purrr::map2_dfr(
      skip_rows,
      table_rows,
      ~vroom::vroom_fwf(
        path_stat,
        vroom::fwf_empty(
          path_stat,
          skip = .x,
          n = .y,
          col_names = var_names
          ),
        skip = .x,
        n_max = .y,
        col_types = col_types
        ),
      .id = table_duplicates
      )
  }

  purrr::map2_dfr(
    ls_files[[file_type]],
    names(ls_files[[file_type]]),
    ~read_cameca_table(
      .x,
      min_row[[1]][[.y]],
      table_row[[.y]],
      var_names,
      col_types
      ),
    .id = "file.nm"
    )
}

# cumulated count according to Cameca software
cameca_stat_X <- read_cameca(
  directory = point_example("2018-01-19-GLENDON"),
  pattern_begin = "Mass# Cumulated count",
  pattern_end = "Ratio#    Ratios",
  file_type = "stat",
  col_types = "dd",
  table_duplicates = "correction block"
  )

usethis::use_data(cameca_stat_X, overwrite = TRUE)


# isotope ratios Cameca software
cameca_stat_R <- read_cameca(
  directory = point_example("2018-01-19-GLENDON"),
  pattern_begin = "Ratio#    Ratios",
  table_depth = 5,
  file_type = "stat",
  col_types = "cddddd",
  table_duplicates = "correction block"
)

usethis::use_data(cameca_stat_R, overwrite = TRUE)
