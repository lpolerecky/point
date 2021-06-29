## Names for the data and machine parameter settings
names_cameca <- readr::read_csv(
  "data-raw/names_cameca.csv",
  col_types = readr::cols(.default = readr::col_character)
  )

usethis::use_data(names_cameca, overwrite = TRUE)
