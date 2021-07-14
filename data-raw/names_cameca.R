## Names for the data and machine parameter settings
names_cameca <- readr::read_csv(
  "data-raw/names_cameca.csv"
  )

usethis::use_data(names_cameca, overwrite = TRUE)
