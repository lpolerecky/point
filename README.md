
<!-- README.md is generated from README.Rmd. Please edit that file -->

# point

<!-- badges: start -->

<!-- badges: end -->

This projects was originally inspired by the lack of detailed insight in
the inner workings of the default software for the Cameca NanoSIMS50L
(Utrecht University). Hence this project has the objective of processing
raw count data into ion and isotope ratios of point-sourced
measurements, and to establish the internal and external precision of
individual measurements and complete series of measurements,
respectively

## Installation

You can install the released version of point

<!-- ``` r -->

<!-- install.packages("point") -->

<!-- ``` -->

## Example

This is a basic example which shows you how to solve a common problem:

``` r

devtools::load_all(".")
#> Loading point
# Use system.file() to access the examples bundled with this package in the
# inst/extdata directory. The examples directories are named:
# 2020-01-17-TREASURE and "2018-01-19-GLENDON"

# raw data containing 13C and 12C counts on carbonate
tb.rw <- read_IC(system.file("extdata", "2018-01-19-GLENDON", package = "point"))

# processing raw ion count data
tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt, deadtime = 44)

# single ion descriptive an predictive statistics for all measured ions
tb.Xt <- stat_Xt(tb.pr, Xt.pr, N.pr, file.nm, species.nm)

# descriptive an predictive statistics for 13C/12C ratios
tb.R <- stat_R(tb.pr, Xt.pr, N.pr, ID = "ID", ion1 = "13C", ion2 = "12C",
               file.nm, species.nm, latex = TRUE)



knitr::kable(tb.R[[1]] %>% rename(!!tb.R[[2]]), escape = FALSE) 
```

| file.nm                          | R.nm    |    n | \(\\bar{x}\) |     \(r\) |   \(s_x\) | \(\epsilon_x\) | \(s_\bar{x}\) | \(\epsilon_\bar{x}\) | \(\hat{s}_x\) | \(\hat{\epsilon}_x\) | \(\hat{s}_\bar{x}\) | \(\hat{\epsilon}_\bar{x}\) | \(\chi^2\) |
| :------------------------------- | :------ | ---: | -----------: | --------: | --------: | -------------: | ------------: | -------------------: | ------------: | -------------------: | ------------------: | -------------------------: | ---------: |
| 2018-01-19-GLENDON\_1\_1.is\_txt | 13C/12C | 3900 |    0.0109721 | 0.6777798 | 0.0010201 |       92.96876 |      1.63e-05 |             1.488692 |     0.0010179 |             92.76863 |            1.63e-05 |                   1.485487 |  1.0043193 |
| 2018-01-19-GLENDON\_1\_2.is\_txt | 13C/12C | 3900 |    0.0109778 | 0.7267907 | 0.0007774 |       70.81376 |      1.24e-05 |             1.133928 |     0.0007697 |             70.11032 |            1.23e-05 |                   1.122664 |  1.0201674 |
| 2018-01-19-GLENDON\_1\_3.is\_txt | 13C/12C | 3900 |    0.0110027 | 0.5665313 | 0.0007318 |       66.51139 |      1.17e-05 |             1.065035 |     0.0007212 |             65.54477 |            1.15e-05 |                   1.049556 |  1.0297123 |
| 2018-01-19-GLENDON\_2\_1.is\_txt | 13C/12C | 4000 |    0.0107790 | 0.9615434 | 0.0007599 |       70.49903 |      1.20e-05 |             1.114688 |     0.0007601 |             70.51348 |            1.20e-05 |                   1.114916 |  0.9995903 |
| 2018-01-19-GLENDON\_2\_2.is\_txt | 13C/12C | 4000 |    0.0107809 | 0.9553736 | 0.0007530 |       69.84868 |      1.19e-05 |             1.104405 |     0.0007388 |             68.52923 |            1.17e-05 |                   1.083542 |  1.0388783 |
| 2018-01-19-GLENDON\_2\_3.is\_txt | 13C/12C | 4000 |    0.0107487 | 0.9623700 | 0.0007365 |       68.52177 |      1.16e-05 |             1.083424 |     0.0007327 |             68.16257 |            1.16e-05 |                   1.077745 |  1.0105675 |
| 2018-01-19-GLENDON\_2\_4.is\_txt | 13C/12C | 4000 |    0.0107693 | 0.9377646 | 0.0007194 |       66.79846 |      1.14e-05 |             1.056176 |     0.0007248 |             67.30553 |            1.15e-05 |                   1.064194 |  0.9849892 |
| 2018-01-19-GLENDON\_2\_5.is\_txt | 13C/12C | 4000 |    0.0108907 | 0.9715567 | 0.0007725 |       70.93615 |      1.22e-05 |             1.121599 |     0.0007746 |             71.12562 |            1.22e-05 |                   1.124595 |  0.9946793 |
| 2018-01-19-GLENDON\_2\_6.is\_txt | 13C/12C | 4000 |    0.0109161 | 0.9940903 | 0.0007044 |       64.52575 |      1.11e-05 |             1.020242 |     0.0007134 |             65.35719 |            1.13e-05 |                   1.033388 |  0.9747189 |
| 2018-01-19-GLENDON\_2\_7.is\_txt | 13C/12C | 4000 |    0.0110876 | 0.9713950 | 0.0008162 |       73.61425 |      1.29e-05 |             1.163944 |     0.0008101 |             73.06514 |            1.28e-05 |                   1.155261 |  1.0150874 |
| 2018-01-19-GLENDON\_2\_8.is\_txt | 13C/12C | 4000 |    0.0110901 | 0.9580354 | 0.0007566 |       68.22399 |      1.20e-05 |             1.078716 |     0.0007437 |             67.06142 |            1.18e-05 |                   1.060334 |  1.0349721 |
| 2018-01-19-GLENDON\_2\_9.is\_txt | 13C/12C | 4000 |    0.0110818 | 0.9642637 | 0.0007135 |       64.38828 |      1.13e-05 |             1.018068 |     0.0007174 |             64.73397 |            1.13e-05 |                   1.023534 |  0.9893481 |
| 2018-01-19-GLENDON\_3\_1.is\_txt | 13C/12C | 4000 |    0.0109466 | 0.9587252 | 0.0009169 |       83.75694 |      1.45e-05 |             1.324313 |     0.0009190 |             83.95670 |            1.45e-05 |                   1.327472 |  0.9952471 |
| 2018-01-19-GLENDON\_3\_2.is\_txt | 13C/12C | 4000 |    0.0107554 | 0.9334001 | 0.0007319 |       68.05230 |      1.16e-05 |             1.076001 |     0.0007229 |             67.21633 |            1.14e-05 |                   1.062783 |  1.0250287 |
| 2018-01-19-GLENDON\_3\_3.is\_txt | 13C/12C | 4000 |    0.0107471 | 0.9881505 | 0.0008391 |       78.07477 |      1.33e-05 |             1.234471 |     0.0008246 |             76.72916 |            1.30e-05 |                   1.213194 |  1.0353818 |
