## code to prepare `DATASET` dataset goes here

# https://wwwrcamnl.wr.usgs.gov/isoig/isopubs/itchch2.html
reference_R <- tibble::tibble(
  isotope = c("2H", "3He","6Li", "11B", "13C", "15N",  "18O", "18O", "34S",
              "37Cl"),
  ratio = c("2H/1H", "3He/4He", "6Li/7Li", "11B/10B", "13C/12C", "15N/14N",
            "18O/16O", "18O/16O", "34S/32S", "37Cl/35Cl"),
  reference = c("VSMOW", "atmospheric He", "L-SVEC", "NBS 951", "VPDB",
                "atmospheric N2", "VSMOW", "VPDB", "CDT", "SMOC"),
  value = c(1.5575e-4, 1.3e-6, 8.32e-2, 4.04362, 1.1237e-2, 3.677e-3,
            2.0052e-3, 2.0672e-3, 4.5005e-2, 	0.324)
)

usethis::use_data(reference_R, overwrite = TRUE)
