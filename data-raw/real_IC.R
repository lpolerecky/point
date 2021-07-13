## code to prepare `real_IC` dataset

#-------------------------------------------------------------------------------
# Point workflow
#-------------------------------------------------------------------------------

# load with meta
tb_rw <- read_IC(point_example("2018-01-19-GLENDON"), meta = TRUE)

real_IC <- cor_IC(tb_rw)

# Save small data-set for package
usethis::use_data(real_IC, overwrite = TRUE, compress = "xz")
