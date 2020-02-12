

read_IC <- function(directory){

  # extract txt files with count data blocks of each single point measuremen
  l.c <- dir(paste0(data.path, directory),
             pattern= ".is_txt") %>%
    # set names for subsequent storage
    set_names()

  # extract stat files with diagnostics of the machine and statistics
  # l.s <- dir(paste0(data.path, directory),
  #            pattern= ".stat$") %>%
  #   # set names for subsequent storage
  #   set_names()
  #
  # x.stat <- map_df(l.s, ~(read_lines(paste0(data.path, directory, "/", .),
  #                                    skip_empty_rows = TRUE,
  #                                    n_max = 22)) %>%
  #                    list() %>% enframe(),
  #                  .id = "file_nm") %>%
  #   mutate(file_nm = str_sub(file_nm, end = (str_length(file_nm) - 6)))

  # function for extracting metadata
  # num_ext <- function(string, regular) {
  #   if(str_detect(string, regular)){
  #     string %>%
  #       str_subset(. , regular) %>%
  #       str_split(. ,":|=") %>%
  #       flatten() %>%
  #       str_trim() %>%
  #       str_subset(., "^[0-9]") %>%
  #       as.numeric()} else {NA}
  # }

  # x.stat <- x.stat %>%
  #   mutate(
  #     Bl.n = map_dbl(value, num_ext, "Block number"),
  #     Bl.m = map_dbl(value, num_ext, "Meas. per block"),
  #     FC.b = map_dbl(value, num_ext, "FC Background before acq"),
  #     FC.a = map_dbl(value, num_ext, "FC Background after acq")
  #   )

  # collecting measurement data and metadata
  x.tb <- left_join(
    bind_cols(
      # raw count data of measrument
      map_df(l.c, ~read_tsv(paste0(data.path, directory, "/", .),
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
      map_df(l.c, ~(read_lines(paste0(data.path, directory, "/", .),
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
    map_df(l.c, ~read_lines(
      paste0(data.path, directory, "/", .),
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


  # empty directory list
  # l.c <- NULL
  # return(x.stat)
  return(x.tb)
  # x <-left_join(x.tb, x.stat, by = "file_nm")
  # return(x)
}

