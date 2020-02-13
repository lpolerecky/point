#' Read raw ion count data
#'
#' @param directory A character string
#'
#' @return A tibble containing raw ion count data
#' @examples
#' read_IC("2018-01-19-GLENDON")
#'
#' @importFrom dplyr
#' @importFrom purrr
#' @importFrom stringr
#' @importFrom tibble
#' @importFrom tidyr
#'
#' @export


read_IC <- function(directory){

  l.c <- read_validator(directory)

# collecting measurement data and metadata
  x.tb <- left_join(
    bind_cols(
# raw count data of measrument
      map_df(l.c, ~read_tsv(paste0(directory, "/", .),
                            col_names = c("X1", "t", "N"),
                            col_types = cols(X1 = col_skip(),
                                             t = col_character(),
                                             N = col_character()
                            ),
                            comment = "B",
                            skip = 1
      ), .id = "file_nm") %>%
# remove old column headers
        filter(t != "X", N  != "Y") %>%
# coercion to numeric values
        mutate(t = as.numeric(t), N = as.numeric(N)),

# saving the metadata of measurement
      map_df(l.c, ~(read_lines(paste0(directory, "/", .),
                               skip_empty_rows = TRUE)  %>%
                      enframe(name = "n") %>%
                      filter(str_detect(value,"(^\\nB) ([^ ]+) ([^ ]+)")) %>%
                      mutate(n = unique(
                        diff(unique(n) -
                               cumsum(
                                 c(2,
                                   rep(2,
                                       length(
                                         unique(n)) - 1)
                                 )
                               )
                        )
                      )
                      ) %>%
                      rowwise() %>%
                      mutate(value = list(rep(value, n))) %>%
                      unnest(cols = c(value))))



    ),
# saving the date of measurement
    map_df(l.c, ~read_lines(paste0(directory, "/", .),
                            skip_empty_rows = TRUE)  %>%
                   enframe(name = "n"),
           .id = "file_nm") %>%
      filter(n == 1) %>%
      select(file_nm, sample_date = value),

  by = "file_nm"

  ) %>%

    # creating readible output
    separate(sample_date,
             into = c("sample", "date"),
             sep = "\\t\\t\\t\\t") %>%
    mutate(file_nm = str_sub(file_nm, end = (str_length(file_nm) - 7)),
           date = as.POSIXct(date, format = c("%d.%m.%y / %H:%M")),
           sample = str_trim(str_replace(sample, "Sample : ", ""))) %>%
    separate(value,
             into = c("Mag", "Rad", "M", "Tc (ms)"), sep = ":+") %>%
    mutate(Mag = as.double(str_extract(Mag, "\\d+.\\d+")),
           Rad = as.double(str_extract(Rad, "\\d+.\\d+")),
           `Tc (ms)` = as.double(str_extract(`Tc (ms)`, "\\d+.\\d+"))) %>%
    separate(M,
             into = c("M", "species"), sep = "\\(+")  %>%
    mutate(M = as.double(str_extract(M, "\\d+.\\d+")),
           species = str_trim(str_replace(species, "\\)", ""))) %>%
    # count intervals (sec) and count ratios (count/sec)
    group_by(file_nm, species) %>%
    mutate(dt = min(diff(t)),
           # count rates
           Xt = N / dt ) %>%
    ungroup() %>%
    # add suffix to raw data for easy filtering later on
    rename_at(vars(N, t, Xt), ~paste0(., ".rw")) %>%

    # creation of unique ID within count which links the different ions
    mutate(ID = paste(file_nm, t.rw, sep = "/"))



  return(x.tb)

}

#' Validate directory of raw ion count data
#'
#' @param directory A character string
#'
#' @return A named vector of file names
#' @examples
#' read_validator("2018-01-19-GLENDON")
#'
#' @importFrom dplyr
#' @importFrom purrr
#' @importFrom stringr
#' @importFrom tibble
#' @importFrom tidyr
#'

read_validator <- function(directory){

  # argument class check
  stopifnot(is_character(directory))

  # extract txt files with count data blocks of each single point measurement
  l.c <- dir(directory,
             pattern = ".is_txt") %>%
           set_names()

  # length check of txt files
    if (any(map_dbl(l.c, ~
                    length(read_lines(paste0(directory, "/", .))
                           )) == 0)) {
      good <- map_dbl(l.c, ~
                        length(read_lines(paste0(directory, "/", .)))) > 0

      l.c <- l.c[good]

      warning("empty txt file removed")

    }

  # column content check of txt files
    if (any(map_dbl(l.c, ~
               nrow(read_tsv(paste0(directory, "/", .),
                             comment = "B",
                             skip = 1,
                             col_names = c("X1", "t", "N"),
                             col_types = cols(X1 = col_skip(),
                                              t = col_character(),
                                              N = col_character()),
                             n_max = 2) %>%
                              filter(t != "X", N  != "Y"))) == 0)){

       good <- map_dbl(l.c, ~
                         nrow(read_tsv(paste0(directory, "/", .),
                                       comment = "B",
                                       skip = 1,
                                       col_names = c("X1", "t", "N"),
                                       col_types = cols(X1 = col_skip(),
                                                        t = col_character(),
                                                        N = col_character()),
                                       n_max = 2) %>%
                                filter(t != "X", N  != "Y"))) > 0

       l.c <- l.c[good]

       warning("txt file contains empty columns")


    }

  return(l.c)
}


