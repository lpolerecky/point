names_diag <- tibble::tibble(
  argument = c(".execution", ".flag"),
  point = c("execution", "flag"),
  use = rep("diagnostics", 2)
)

usethis::use_data(names_diag, overwrite = TRUE)
