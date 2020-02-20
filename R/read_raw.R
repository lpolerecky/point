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

# collecting metadata (stat file)
  tb.meta <- read_meta(directory)

#   x.meta <- left_join(
#                 x.meta ,
#
# # saving the metadata of measurement from the raw coun file
#                 map_df(l.c, ~(read_lines(paste0(directory, "/", .),
#                               skip_empty_rows = TRUE,
# # n max = 3 lines for the first species and 2 for each consecutive
#                               n_max = with(x.meta, (unique(n[file_nm == .])) *
#                                              length(unique(ions)) +
#                                              length(unique(ions)) * 2 + 1)
#                               ) %>%
#                   enframe(name = NULL) %>%
#                   filter(str_detect(value,"(^\\nB) ([^ ]+) ([^ ]+)"))),
#                   .id = "file_nm")  %>%
#                   separate(value,
#                            into = c("Mag", "Rad", "M", "Tc (ms)"), sep = ":+") %>%
#                   mutate(Mag = as.double(str_extract(Mag, "\\d+.\\d+")),
#                   Rad = as.double(str_extract(Rad, "\\d+.\\d+")),
#                   `Tc (ms)` = as.double(str_extract(`Tc (ms)`, "\\d+.\\d+"))) %>%
#                   separate(M,
#                            into = c("M", "species"), sep = "\\(+")  %>%
#                   mutate(M = as.double(str_extract(M, "\\d+.\\d+")),
#                          species = str_trim(str_replace(species, "\\)", ""))),
#
#
#               by = "file_nm")

# collecting measurement data and metadata
  tb.rw <- left_join(

# raw count data of measrument
      map_df(l.c,
             ~read_tsv(paste0(directory, "/", .),
                       col_names = c("X1", "t.rw", "N.rw"),
                       col_types = cols(X1 = col_skip(),
                                        t.rw = col_character(),
                                        N.rw = col_character()),
                        comment = "B",
                        skip = 1,
# n-max is n times number of species
                        n_max = with(tb.meta, (unique(n.rw[file.nm == .]) + 1) *
                                      length(unique(species.mt)))
                       ),
             .id = "file.nm") %>%
# remove old column headers
        filter(t.rw != "X", N.rw  != "Y") %>%
# coercion to numeric values
        mutate(t.rw = as.numeric(t.rw), N.rw = as.numeric(N.rw)) %>%
        group_by(file.nm) %>%
        mutate(num.mt = ntile(n = with(tb.meta, length(unique(species.mt))))) %>%
        ungroup(),


      tb.meta,

      by = c("file.nm", "num.mt")) %>%
# creation of unique ID within count which links the different ions
        mutate(ID = paste(file.nm, t.rw, sep = "/"))


}




read_validator <- function(directory){

# argument class check
  stopifnot(is_character(directory))

# extract txt files with count data blocks of each single point measurement
  l.c <- dir(directory,
             pattern = ".is_txt") %>%
           set_names()

# length check of txt files
    if (any(map_dbl(l.c, ~length(read_lines(paste0(directory, "/", .), n_max = 2)
                           )) == 0)) {

      good <- map_dbl(l.c, ~length(read_lines(paste0(directory, "/", .), n_max = 2))) > 0

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

         good <- map_dbl(l.c, ~nrow(read_tsv(paste0(directory, "/", .),
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


read_meta <- function(directory){

# NA aliases
  NA_aliases <- c("N/A", "none", "None") %>%
                  set_names( rep(NA_character_, length(.)), nm = .)

# extract stat files with diagnostics of the machine and statistics
  l.s <- dir(directory,
             pattern = ".stat$") %>%
# set names for subsequent storage
           set_names() %>%
# remove transect files
           discard(., str_detect(., "transect.stat"))

  min_n <- lapply(map(l.s, ~read_lines(paste0(directory, "/", .),
                                       skip_empty_rows = TRUE,
                                       n_max = 50)),
                  str_which, "#") %>%
            map(., 1) %>%
              flatten_dbl()

  max_n <- lapply(map(l.s, ~read_lines(paste0(directory, "/", .),
                                     skip_empty_rows = TRUE,
                                     n_max = 50)),
                str_which, "--") %>%
            map(., 1) %>%
            flatten_dbl()


# ions and MS setup
  l <- lst(a = l.s, b =  min_n + 3 , c = max_n - b + 1)
  f <- function(a, b, c) {

        read_table(paste0(directory, "/", a),
                   skip_empty_rows = TRUE, skip = b,
                   n_max = c,
                   col_names = c("num.mt", "species.mt", "mass.mt", "det.mt", "tc.mt (s)", "bfield.mt", "rad.mt", "X2", "X3", "X4"),
                   col_types = cols(
                                    num.mt = col_integer(),
                                    species.mt = col_character(),
                                    mass.mt = col_double(),
                                    det.mt = col_character(),
                                    `tc.mt (s)` = col_double(),
                                    bfield.mt = col_double(),
                                    rad.mt = col_double(),
                                    X2  = col_skip(),
                                    X3 = col_skip(),
                                    X4 = col_skip()
                   ))

  }


  tb.ion <- pmap_dfr(l, f, .id = "file.nm") %>%
              mutate(file.nm = str_sub(file.nm,
                                       end = (str_length(file.nm) - 5)))

# primary and secondary ion beam metadata

  tb.meas <- map2_df(l.s, min_n - 2, ~(read_lines(paste0(directory, "/", .x),
                                     skip_empty_rows = TRUE,
                                     n_max = .y)) %>%
                            list() %>%
                            enframe(name = NULL),
                   .id = "file.nm") %>%
              mutate(file.nm = str_sub(file.nm,
                                       end = (str_length(file.nm) - 5)))

# function to make metadate readible
  str_unfold <-  function(string) {

    string %>%
# create tibble from list
      enframe(name = NULL)  %>%
# remove date
      filter(!str_detect(value, pattern  = "\n(?=\t)")) %>%
      separate_rows(value, sep = "(/(?=[:blank:])) | =") %>%
      separate(value,
               into = c("variable", "value"),
               sep =":(?!\\\\)",
               extra = "merge") %>%
      mutate_all(str_trim) %>%
      distinct(., variable, .keep_all = TRUE) %>%
# convert NA aliases to NA
      mutate(value = recode(value, !!!NA_aliases),
# remove units behind numerics
             value =str_replace(value, "(?<=[:digit:]|[:blank:])(pA|um|%)", "")) %>%
      pivot_wider(names_from = variable, values_from = value) %>%
# separate date
      mutate(date = as.POSIXct(
                        str_replace_all(
                                        string[str_detect(string,
                                               pattern  = "\n(?=\t)")],
                                        "\\t|\\n", ""),
                       format = c("%d.%m.%y  %H:%M")
                       )) %>%
      rename(sample.nm = "CAMECA \\ ISOTOPES \\ Sample") %>%
# remove extra dot
      mutate(`Pre Sputtering Time (s)` = str_sub(`Pre Sputtering Time (s)`, end = -2))

  }

# named vector for renaming variables
  meta.nm  <- c(
    `presput.mt (s)` = "Pre Sputtering Time (s)",
    `bl_num.mt (n)` = "Block number",
    `meas_bl.mt (n)` = "Meas. per block",
    `width_hor.mt (V)` = "Width Horizontal(V)",
    `width_ver.mt (V)` =  "Vertical(V)",
    `prim_cur_start.mt (pA)` = "Primary Current before acq",
    `prim_cur_after.mt (pA)` = "after acq",
    `rast_com.mt (um)` = "Raster (um)",
    `blank_rast.mt (%)` = "Blanking"

  )

  tb.meas <- tb.meas %>%
               mutate(value = map(value, str_unfold)) %>%
               unnest(cols = c(value)) %>%
               select(c(file.nm, date, !!! meta.nm[meta.nm %in% c(colnames(tb.ion), colnames(.))])) %>%
               mutate_at(vars(one_of(names(meta.nm))), as.double) %>%
# add measurement number
               mutate(n.rw = `bl_num.mt (n)` * `meas_bl.mt (n)`)

# combine MS and beam metadata

  tb.meta <- left_join(tb.ion, tb.meas, by = "file.nm") %>%
               mutate(file.nm = paste0(file.nm, ".is_txt"))

  }


# function to read CAMECA output to validate point output (future work!!)
# read_test

