
<!--  use the --webtex argument in the YAML to render equations -->

# Introduction to point

<!-- badges: start -->

<!-- badges: end -->

This project was originally inspired by the lack of detailed insight in
the inner workings of the default software for the *Cameca NanoSIMS50L*
(Utrecht University). Hence this project has the objective of processing
raw count data into ion and isotope ratios of point-sourced
measurements; and to establish the internal and external precision of,
respectively, individual analyses and complete series of analyses.
Access to raw ion count data is useful as it allows detection of
anomolous values associated with e.g. machine instability or
heterogeneity of the analysed sample. Upon detection, anomolous values
can be omitted or further analysed to delineate the source of variation.

## Installation

You can install the released version of point

``` r
# Install point rom GitHub:
# install.packages("devtools")
devtools::install_github("point")
```

## Credits

The construction of the R (R Core Team 2020) package *point* and
associated documentation was aided by the packages; *devtools* (Wickham
et al. 2019c), *roxygen2* (Wickham et al. 2019a), *testthat* (Wickham
2011), *knitr* (Xie 2015, 2020 ), *rmarkdown* (Xie et al. 2018, Allaire
et al. 2019), and the superb guidance in the book: *R packages:
organize, test, document, and share your code*, by Wickham (2015). In
addition, this package relies on a set of external packages from the
tidyverse universe, including: *dplyr* (Wickham et al. 2019b), *tidyr*
(Wickham & Henry 2019), *tibble* (Müller & Wickham 2019), *stringr*
(Wickham 2019), *readr* (Wickham et al. 2018), *magrittr* (Bache &
Wickham 2014), *ggplot2* (Wickham 2016), *rlang* (Henry & Wickham 2020),
and *purrr* (Henry & Wickham 2019) for internal functioning as well as
specialised statistics; *polyaAeppli* (Burden 2014).

## Usage

Load point with `library`.

``` r
library(point)
```

## The point workflow

A more detailed outline of the general point workflow is given in the
vignette
[IC-introduction](https://github.com/MartinSchobben/point/tree/master/vignettes/IC-introduction.html).

<img src="vignettes/workflow.png" width="100%" />

To read, process and analyse raw ion count data use the functions:

  - `read_IC`: raw ion count data
  - `cor_IC`: process ion count data
  - `stat_Xt`: analyse single ion count data
  - `stat_R`: analyse ion ratios

## Example 1: internal precision of isotope ratios

This is an example of how Cameca NanoSIMS50L raw datafiles can be
extracted, processed and analysed for the <sup>13</sup>C/<sup>12</sup>C
isotope ratio (![R](https://latex.codecogs.com/png.latex?R "R")). This
produces a tibble with descriptive and predictive poisson statistics
(demarcated with an
![\\\\\\hat{\\\\\\phantom{,}}](https://latex.codecogs.com/png.latex?%5C%5C%5Chat%7B%5C%5C%5Cphantom%7B%2C%7D%7D
"\\\\\\hat{\\\\\\phantom{,}}") ) of the ion count data. This can be done
for single count blocks in order to obtain internal
precision.

``` r
# Use point_example() to access the examples bundled with this package in the
# inst/extdata directory.

# raw data containing 13C and 12C counts on carbonate
tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))

# processing raw ion count data
tb.pr <- cor_IC(tb.rw, 
                N = N.rw, 
                t = t.rw, 
                Det = det_type.mt, 
                deadtime = 44, 
                thr_PHD = 50)

# descriptive an predictive statistics for 13C/12C ratios
tb.R <- stat_R(tb.pr, 
               Xt = Xt.pr, 
               N = N.pr, 
               species = species.nm, 
               ion1 = "13C", 
               ion2 = "12C", 
               sample.nm, 
               file.nm, 
               latex = TRUE)
```

| sample.nm                        | file.nm                   | R                                                                                                                                                                                                                                                                                                                                                   | ![n](https://latex.codecogs.com/png.latex?n "n") | ![\\bar{R}](https://latex.codecogs.com/png.latex?%5Cbar%7BR%7D "\\bar{R}") | ![s\_{R}](https://latex.codecogs.com/png.latex?s_%7BR%7D "s_{R}") | ![\\epsilon\_{R} \\,](https://latex.codecogs.com/png.latex?%5Cepsilon_%7BR%7D%20%5C%2C "\\epsilon_{R} \\,") (‰) | ![s\_{\\bar{R}}](https://latex.codecogs.com/png.latex?s_%7B%5Cbar%7BR%7D%7D "s_{\\bar{R}}") | ![\\epsilon\_{\\bar{R}} \\,](https://latex.codecogs.com/png.latex?%5Cepsilon_%7B%5Cbar%7BR%7D%7D%20%5C%2C "\\epsilon_{\\bar{R}} \\,") (‰) | ![\\hat{s}\_{R}](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_%7BR%7D "\\hat{s}_{R}") | ![\\hat{\\epsilon}\_{R} \\,](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cepsilon%7D_%7BR%7D%20%5C%2C "\\hat{\\epsilon}_{R} \\,") (‰) | ![\\hat{s}\_{\\bar{R}}](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_%7B%5Cbar%7BR%7D%7D "\\hat{s}_{\\bar{R}}") | ![\\hat{\\epsilon}\_{\\bar{R}} \\,](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cepsilon%7D_%7B%5Cbar%7BR%7D%7D%20%5C%2C "\\hat{\\epsilon}_{\\bar{R}} \\,") (‰) | ![\\chi^{2}](https://latex.codecogs.com/png.latex?%5Cchi%5E%7B2%7D "\\chi^{2}") |
| :------------------------------- | :------------------------ | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -----------------------------------------------: | -------------------------------------------------------------------------: | ----------------------------------------------------------------: | --------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------: | ----------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------: | ----------------------------------------------------------------------------------------------------------------------------------------: | --------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------: |
| Belemnite, Indium                | 2018-01-19-GLENDON\_1\_1  | ![\\phantom{,}^{13}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![\\phantom{,}^{12}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |                                             3900 |                                                                    0.01097 |                                                          0.001020 |                                                                                                            93.0 |                                                                                    1.63e-05 |                                                                                                                                      1.49 |                                                                                    0.001017 |                                                                                                                                      92.8 |                                                                                                              1.63e-05 |                                                                                                                                                                1.49 |                                                                           1.004 |
| Belemnite, Indium                | 2018-01-19-GLENDON\_1\_2  | ![\\phantom{,}^{13}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![\\phantom{,}^{12}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |                                             3900 |                                                                    0.01097 |                                                          0.000777 |                                                                                                            70.8 |                                                                                    1.24e-05 |                                                                                                                                      1.13 |                                                                                    0.000769 |                                                                                                                                      70.1 |                                                                                                              1.23e-05 |                                                                                                                                                                1.12 |                                                                           1.020 |
| Belemnite, Indium                | 2018-01-19-GLENDON\_1\_3  | ![\\phantom{,}^{13}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![\\phantom{,}^{12}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |                                             3900 |                                                                    0.01099 |                                                          0.000731 |                                                                                                            66.5 |                                                                                    1.17e-05 |                                                                                                                                      1.07 |                                                                                    0.000721 |                                                                                                                                      65.5 |                                                                                                              1.15e-05 |                                                                                                                                                                1.05 |                                                                           1.030 |
| W12, Indium (W12\_Indium\_1\_1)  | 2018-01-19-GLENDON\_2\_1  | ![\\phantom{,}^{13}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![\\phantom{,}^{12}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |                                             4000 |                                                                    0.01077 |                                                          0.000759 |                                                                                                            70.5 |                                                                                    1.20e-05 |                                                                                                                                      1.11 |                                                                                    0.000759 |                                                                                                                                      70.5 |                                                                                                              1.20e-05 |                                                                                                                                                                1.11 |                                                                           1.000 |
| W12, Indium (W12\_Indium\_1\_10) | 2018-01-19-GLENDON\_2\_10 | ![\\phantom{,}^{13}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![\\phantom{,}^{12}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |                                             4000 |                                                                    0.00503 |                                                          0.318830 |                                                                                                         63447.2 |                                                                                    5.04e-03 |                                                                                                                                   1003.19 |                                                                                    0.318614 |                                                                                                                                   63404.3 |                                                                                                              5.04e-03 |                                                                                                                                                             1002.51 |                                                                           1.001 |
| W12, Indium (W12\_Indium\_1\_10) | 2018-01-19-GLENDON\_3\_1  | ![\\phantom{,}^{13}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![\\phantom{,}^{12}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |                                             4000 |                                                                    0.01094 |                                                          0.000916 |                                                                                                            83.8 |                                                                                    1.45e-05 |                                                                                                                                      1.32 |                                                                                    0.000919 |                                                                                                                                      84.0 |                                                                                                              1.45e-05 |                                                                                                                                                                1.33 |                                                                           0.996 |

## Example 2: external precision of isotope ratios

To calculate the external reproducibility of isotope ratios one needs to
use the total ion count of one block and the block count rate. The
latter is equivalent to the mean ion count rate, which can be calculated
with the function `stat_Xt`.

``` r
# single ion descriptive an predictive statistics for all measured ions
tb.Xt <- stat_Xt(tb.pr, 
                 Xt = Xt.pr, 
                 N = N.pr, 
                 species = species.nm, 
                 sample.nm, 
                 file.nm,
                 latex = FALSE)

# For this particular run a belemnite was used as reference material. 
tb.R.ext  <- stat_R(tb.Xt, 
                    Xt = M_Xt.pr, 
                    N = Ntot_Xt.pr, 
                    species = species.nm, 
                    ion1 = "13C", 
                    ion2 = "12C",  
                    sample.nm, 
                    latex = TRUE)
```

| sample.nm         | R                                                                                                                                                                                                                                                                                                                                                   | ![n](https://latex.codecogs.com/png.latex?n "n") | ![\\bar{R}](https://latex.codecogs.com/png.latex?%5Cbar%7BR%7D "\\bar{R}") | ![s\_{R}](https://latex.codecogs.com/png.latex?s_%7BR%7D "s_{R}") | ![\\epsilon\_{R} \\,](https://latex.codecogs.com/png.latex?%5Cepsilon_%7BR%7D%20%5C%2C "\\epsilon_{R} \\,") (‰) | ![s\_{\\bar{R}}](https://latex.codecogs.com/png.latex?s_%7B%5Cbar%7BR%7D%7D "s_{\\bar{R}}") | ![\\epsilon\_{\\bar{R}} \\,](https://latex.codecogs.com/png.latex?%5Cepsilon_%7B%5Cbar%7BR%7D%7D%20%5C%2C "\\epsilon_{\\bar{R}} \\,") (‰) | ![\\hat{s}\_{R}](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_%7BR%7D "\\hat{s}_{R}") | ![\\hat{\\epsilon}\_{R} \\,](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cepsilon%7D_%7BR%7D%20%5C%2C "\\hat{\\epsilon}_{R} \\,") (‰) | ![\\hat{s}\_{\\bar{R}}](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_%7B%5Cbar%7BR%7D%7D "\\hat{s}_{\\bar{R}}") | ![\\hat{\\epsilon}\_{\\bar{R}} \\,](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cepsilon%7D_%7B%5Cbar%7BR%7D%7D%20%5C%2C "\\hat{\\epsilon}_{\\bar{R}} \\,") (‰) | ![\\chi^{2}](https://latex.codecogs.com/png.latex?%5Cchi%5E%7B2%7D "\\chi^{2}") |
| :---------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -----------------------------------------------: | -------------------------------------------------------------------------: | ----------------------------------------------------------------: | --------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------: | ----------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------: | ----------------------------------------------------------------------------------------------------------------------------------------: | --------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------------------------------------------------------------------------------------------: | ------------------------------------------------------------------------------: |
| Belemnite, Indium | ![\\phantom{,}^{13}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}")/![\\phantom{,}^{12}](https://latex.codecogs.com/png.latex?%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](https://latex.codecogs.com/png.latex?_%7B%7D "_{}") |                                                3 |                                                                      0.011 |                                                          1.55e-05 |                                                                                                            1.41 |                                                                                     8.9e-06 |                                                                                                                                     0.814 |                                                                                     1.3e-05 |                                                                                                                                      1.18 |                                                                                                               7.5e-06 |                                                                                                                                                               0.681 |                                                                            1.43 |

For more detailed
    information:

  - [IC-read](https://github.com/MartinSchobben/point/tree/master/vignettes/IC-read.html):
    reading raw ion count data
    use
  - [IC-process](https://github.com/MartinSchobben/point/tree/master/vignettes/IC-process.html):
    processing ion count
    data  
  - [IC-precision](https://github.com/MartinSchobben/point/tree/master/vignettes/IC-precision.html):
    statistics concerning ion count precision

# References

<div id="refs" class="references">

<div id="ref-rmarkdown1">

Allaire J, Xie Y, McPherson J, Luraschi J, Ushey K, Atkins A, Wickham H,
Cheng J, Chang W, Iannone R (2019) Rmarkdown: Dynamic documents for r.

</div>

<div id="ref-magrittr">

Bache SM, Wickham H (2014) Magrittr: A forward-pipe operator for r.

</div>

<div id="ref-polyaAeppli">

Burden C (2014) PolyaAeppli: Implementation of the polya-aeppli
distribution.

</div>

<div id="ref-purrr">

Henry L, Wickham H (2019) Purrr: Functional programming tools.

</div>

<div id="ref-rlang">

Henry L, Wickham H (2020) Rlang: Functions for base types and core r and
’tidyverse’ features.

</div>

<div id="ref-tibble">

Müller K, Wickham H (2019) Tibble: Simple data frames.

</div>

<div id="ref-rversion">

R Core Team (2020) R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.

</div>

<div id="ref-ggplot2">

Wickham H (2016) Ggplot2: Elegant graphics for data analysis.
Springer-Verlag New York.

</div>

<div id="ref-Wickham2015">

Wickham H (2015) R packages: Organize, test, document, and share your
code. O’Reilly Media, Inc.

</div>

<div id="ref-stringr">

Wickham H (2019) Stringr: Simple, consistent wrappers for common string
operations.

</div>

<div id="ref-testthat">

Wickham H (2011) Testthat: Get Started with Testing. The R Journal
3:5–10.

</div>

<div id="ref-roxygen2">

Wickham H, Danenberg P, Csárdi G, Eugster M (2019a) Roxygen2: In-line
documentation for r.

</div>

<div id="ref-dplyr">

Wickham H, François R, Henry L, Müller K (2019b) Dplyr: A grammar of
data manipulation.

</div>

<div id="ref-tidyr">

Wickham H, Henry L (2019) Tidyr: Tidy messy data.

</div>

<div id="ref-devtools">

Wickham H, Hester J, Chang W (2019c) Devtools: Tools to make developing
r packages easier.

</div>

<div id="ref-readr">

Wickham H, Hester J, Francois R (2018) Readr: Read rectangular text
data.

</div>

<div id="ref-knitr2">

Xie Y (2015) Dynamic documents with R and knitr, 2nd ed. Chapman;
Hall/CRC, Boca Raton, Florida.

</div>

<div id="ref-knitr1">

Xie Y (2020) Knitr: A general-purpose package for dynamic report
generation in r.

</div>

<div id="ref-rmarkdown2">

Xie Y, Allaire J, Grolemund G (2018) R markdown: The definitive guide.
Chapman; Hall/CRC, Boca Raton, Florida.

</div>

</div>
