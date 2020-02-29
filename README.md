
<!--  use the --webtex argument in the YAML to render equations -->

# Introduction to point

<!-- badges: start -->

<!-- badges: end -->

This projects was originally inspired by the lack of detailed insight in
the inner workings of the default software for the Cameca NanoSIMS50L
(Utrecht University). Hence this project has the objective of processing
raw count data into ion and isotope ratios of point-sourced
measurements, and to establish the internal and external precision of
individual measurements and complete series of measurements,
respectively. Access to raw ion count data is useful as it allows
detection of anomolous values associated with e.g. machine instability
or heterogeneity of the analysed sample. Such anomolous values can then
be omitted, or further analysed to delineate the source of variation.

## Installation

You can install the released version of point

<!-- ``` r -->

<!-- install.packages("point") -->

<!-- ``` -->

## Example1: internal precision of isotope ratios

This is an example of how Cameca NanoSIMS50L raw datafiles can be
extracted, processed and analysed for the <sup>13</sup>C/<sup>12</sup>C
isotope ratio (R). This produces a tibble with descriptive and
predictive poisson statistics (demarcated with an
![\\\\\\hat{}](https://latex.codecogs.com/png.latex?%5C%5C%5Chat%7B%7D
"\\\\\\hat{}")) of the ion count data. This can be done for single count
blocks in order to obtain internal precision.

``` r

# library(point)
devtools::load_all(".")
#> Loading point

# Use point_example() to access the examples bundled with this package in the
# inst/extdata directory. The examples directories are named:
# 2020-01-17-TREASURE and "2018-01-19-GLENDON"

# raw data containing 13C and 12C counts on carbonate
tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))

# processing raw ion count data
tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt, deadtime = 44, thr = 180)

# single ion descriptive an predictive statistics for all measured ions
tb.Xt <- stat_Xt(tb.pr, Xt.pr, N.pr, species.nm, sample.nm, file.nm)

# descriptive an predictive statistics for 13C/12C ratios
tb.R <- stat_R(tb.pr, Xt.pr, N.pr, species.nm, ion1 = "13C", ion2 = "12C", 
               sample.nm, file.nm, latex = TRUE)


# use the second element of the list to get nice variable names rendering in 
# Rmarkdown
knitr::kable(head(tb.R[[1]] %>% 
                    rename(sample = "sample.nm", file = "file.nm", R = "R.nm", 
                                   !!! tb.R[[2]])),
             format.args = list(digits = 5, format = "G", flag = "0")) 
```

| sample                           | file                     | ![n](https://latex.codecogs.com/png.latex?n "n") | ![\\bar{R}](https://latex.codecogs.com/png.latex?%5Cbar%7BR%7D "\\bar{R}") | ![r](https://latex.codecogs.com/png.latex?r "r") | ![s\_{R}](https://latex.codecogs.com/png.latex?s_%7BR%7D "s_{R}") | ![\\epsilon\_{R} (\\text{‰})](https://latex.codecogs.com/png.latex?%5Cepsilon_%7BR%7D%20%28%5Ctext%7B%E2%80%B0%7D%29 "\\epsilon_{R} (\\text{‰})") | ![s\_{\\bar{R}}](https://latex.codecogs.com/png.latex?s_%7B%5Cbar%7BR%7D%7D "s_{\\bar{R}}") | ![\\epsilon\_{\\bar{R}} (\\text{‰})](https://latex.codecogs.com/png.latex?%5Cepsilon_%7B%5Cbar%7BR%7D%7D%20%28%5Ctext%7B%E2%80%B0%7D%29 "\\epsilon_{\\bar{R}} (\\text{‰})") | ![\\hat{s}\_{R}](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_%7BR%7D "\\hat{s}_{R}") | ![\\hat{\\epsilon}\_{R} (\\text{‰})](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cepsilon%7D_%7BR%7D%20%28%5Ctext%7B%E2%80%B0%7D%29 "\\hat{\\epsilon}_{R} (\\text{‰})") | ![\\hat{s}\_{\\bar{R}}](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_%7B%5Cbar%7BR%7D%7D "\\hat{s}_{\\bar{R}}") | ![\\hat{\\epsilon}\_{\\bar{R}} (\\text{‰})](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cepsilon%7D_%7B%5Cbar%7BR%7D%7D%20%28%5Ctext%7B%E2%80%B0%7D%29 "\\hat{\\epsilon}_{\\bar{R}} (\\text{‰})") | ![\\chi^{2}](https://latex.codecogs.com/png.latex?%5Cchi%5E%7B2%7D "\\chi^{2}") | R                                                                                                                                                                                                                                                             |
| :------------------------------- | :----------------------- | -----------------------------------------------: | -------------------------------------------------------------------------: | -----------------------------------------------: | ----------------------------------------------------------------: | ------------------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------: | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------: | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | --------------------------------------------------------------------------------------------------------------------: | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Belemnite, Indium                | 2018-01-19-GLENDON\_1\_1 |                                             3900 |                                                                   0.010906 |                                          0.67778 |                                                         0.0010140 |                                                                                                                                            92.969 |                                                                                    1.62e-05 |                                                                                                                                                                      1.4887 |                                                                                   0.0010117 |                                                                                                                                                                      92.766 |                                                                                                              1.62e-05 |                                                                                                                                                                                                1.4854 |                                                                         1.00438 | ![^{13}](https://latex.codecogs.com/png.latex?%5E%7B13%7D "^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![^{12}](https://latex.codecogs.com/png.latex?%5E%7B12%7D "^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |
| Belemnite, Indium                | 2018-01-19-GLENDON\_1\_2 |                                             3900 |                                                                   0.010928 |                                          0.72679 |                                                         0.0007739 |                                                                                                                                            70.814 |                                                                                    1.24e-05 |                                                                                                                                                                      1.1339 |                                                                                   0.0007661 |                                                                                                                                                                      70.109 |                                                                                                              1.23e-05 |                                                                                                                                                                                                1.1226 |                                                                         1.02022 | ![^{13}](https://latex.codecogs.com/png.latex?%5E%7B13%7D "^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![^{12}](https://latex.codecogs.com/png.latex?%5E%7B12%7D "^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |
| Belemnite, Indium                | 2018-01-19-GLENDON\_1\_3 |                                             3900 |                                                                   0.011001 |                                          0.56653 |                                                         0.0007317 |                                                                                                                                            66.511 |                                                                                    1.17e-05 |                                                                                                                                                                      1.0650 |                                                                                   0.0007210 |                                                                                                                                                                      65.545 |                                                                                                              1.15e-05 |                                                                                                                                                                                                1.0496 |                                                                         1.02971 | ![^{13}](https://latex.codecogs.com/png.latex?%5E%7B13%7D "^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![^{12}](https://latex.codecogs.com/png.latex?%5E%7B12%7D "^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |
| W12, Indium (W12\_Indium\_1\_1)  | 2018-01-19-GLENDON\_2\_1 |                                             4000 |                                                                   0.010729 |                                          0.96154 |                                                         0.0007564 |                                                                                                                                            70.499 |                                                                                    1.20e-05 |                                                                                                                                                                      1.1147 |                                                                                   0.0007565 |                                                                                                                                                                      70.512 |                                                                                                              1.20e-05 |                                                                                                                                                                                                1.1149 |                                                                         0.99964 | ![^{13}](https://latex.codecogs.com/png.latex?%5E%7B13%7D "^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![^{12}](https://latex.codecogs.com/png.latex?%5E%7B12%7D "^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |
| W12, Indium (W12\_Indium\_1\_10) | 2018-01-19-GLENDON\_3\_1 |                                             4000 |                                                                   0.010868 |                                          0.95873 |                                                         0.0009103 |                                                                                                                                            83.757 |                                                                                    1.44e-05 |                                                                                                                                                                      1.3243 |                                                                                   0.0009124 |                                                                                                                                                                      83.953 |                                                                                                              1.44e-05 |                                                                                                                                                                                                1.3274 |                                                                         0.99532 | ![^{13}](https://latex.codecogs.com/png.latex?%5E%7B13%7D "^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![^{12}](https://latex.codecogs.com/png.latex?%5E%7B12%7D "^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |
| W12, Indium (W12\_Indium\_1\_11) | 2018-01-19-GLENDON\_3\_2 |                                             4000 |                                                                   0.010708 |                                          0.93340 |                                                         0.0007287 |                                                                                                                                            68.052 |                                                                                    1.15e-05 |                                                                                                                                                                      1.0760 |                                                                                   0.0007198 |                                                                                                                                                                      67.215 |                                                                                                              1.14e-05 |                                                                                                                                                                                                1.0628 |                                                                         1.02508 | ![^{13}](https://latex.codecogs.com/png.latex?%5E%7B13%7D "^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![^{12}](https://latex.codecogs.com/png.latex?%5E%7B12%7D "^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |

## Example 2: external precision of isotope ratios

To calculate the external reproducibility of isotope ratios one needs to
use the total ion count of one block and the block count rate. The
latter is equivalent to the mean ion count rate, which was already
calculated with the function `stat_Xt`.

``` r
# For this particular run a belemnite was used as reference material. 
tb.R.ext  <- stat_R(tb.Xt, 
                    Xt = M_Xt.pr, 
                    N = Ntot_Xt.pr, 
                    species = species.nm, 
                    ion1 = "13C", 
                    ion2 = "12C",  
                    sample.nm, 
                    latex = TRUE)

knitr::kable(head(tb.R.ext[[1]]  %>% 
                    filter(str_detect(sample.nm, "Belemnite")) %>% 
                    rename(sample = "sample.nm", R = "R.nm", 
                                   !!! tb.R.ext[[2]])),
              format.args = list(digits = 5, format = "G", flag = "0")) 
```

| sample            | ![n](https://latex.codecogs.com/png.latex?n "n") | ![\\bar{R}](https://latex.codecogs.com/png.latex?%5Cbar%7BR%7D "\\bar{R}") | ![r](https://latex.codecogs.com/png.latex?r "r") | ![s\_{R}](https://latex.codecogs.com/png.latex?s_%7BR%7D "s_{R}") | ![\\epsilon\_{R} (\\text{‰})](https://latex.codecogs.com/png.latex?%5Cepsilon_%7BR%7D%20%28%5Ctext%7B%E2%80%B0%7D%29 "\\epsilon_{R} (\\text{‰})") | ![s\_{\\bar{R}}](https://latex.codecogs.com/png.latex?s_%7B%5Cbar%7BR%7D%7D "s_{\\bar{R}}") | ![\\epsilon\_{\\bar{R}} (\\text{‰})](https://latex.codecogs.com/png.latex?%5Cepsilon_%7B%5Cbar%7BR%7D%7D%20%28%5Ctext%7B%E2%80%B0%7D%29 "\\epsilon_{\\bar{R}} (\\text{‰})") | ![\\hat{s}\_{R}](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_%7BR%7D "\\hat{s}_{R}") | ![\\hat{\\epsilon}\_{R} (\\text{‰})](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cepsilon%7D_%7BR%7D%20%28%5Ctext%7B%E2%80%B0%7D%29 "\\hat{\\epsilon}_{R} (\\text{‰})") | ![\\hat{s}\_{\\bar{R}}](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_%7B%5Cbar%7BR%7D%7D "\\hat{s}_{\\bar{R}}") | ![\\hat{\\epsilon}\_{\\bar{R}} (\\text{‰})](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cepsilon%7D_%7B%5Cbar%7BR%7D%7D%20%28%5Ctext%7B%E2%80%B0%7D%29 "\\hat{\\epsilon}_{\\bar{R}} (\\text{‰})") | ![\\chi^{2}](https://latex.codecogs.com/png.latex?%5Cchi%5E%7B2%7D "\\chi^{2}") | R                                                                                                                                                                                                                                                             |
| :---------------- | -----------------------------------------------: | -------------------------------------------------------------------------: | -----------------------------------------------: | ----------------------------------------------------------------: | ------------------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------: | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------: | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | --------------------------------------------------------------------------------------------------------------------: | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Belemnite, Indium |                                                3 |                                                                   0.010954 |                                          0.99995 |                                                          5.09e-05 |                                                                                                                                            4.6445 |                                                                                    2.94e-05 |                                                                                                                                                                      2.6815 |                                                                                    1.29e-05 |                                                                                                                                                                        1.18 |                                                                                                               7.5e-06 |                                                                                                                                                                                               0.68129 |                                                                          15.491 | ![^{13}](https://latex.codecogs.com/png.latex?%5E%7B13%7D "^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![^{12}](https://latex.codecogs.com/png.latex?%5E%7B12%7D "^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |

For more detailed information about reading raw ion count data use
`vignette("IC-read")`; processing ion count data `vignette("IC-read")`;
and for statistical test concerning ion count precision
`vignette("IC-precision")`.
