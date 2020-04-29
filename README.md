
<!--  use the --webtex argument in the YAML to render equations -->

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

# Introduction to point

This project was originally inspired by the lack of detailed insight in
the inner workings of the default software for the *Cameca NanoSIMS50L*
(Utrecht University). Hence this project has the objective of processing
raw count data into ion and isotope ratios of point-sourced
measurements; and to establish the internal and external precision of,
respectively, individual analyses and complete series of analyses.
Access to raw ion count data is useful as it allows detection of
anomalous values associated with e.g. machine instability or
heterogeneity of the analysed sample. Upon detection, anomalous values
can be omitted or further analysed to delineate the source of variation.

## Credits

The construction of the R (R Core Team 2020) package *point* and
associated documentation was aided by the packages; *devtools* (Wickham,
Hester, and Chang 2019), *roxygen2* (Wickham, Danenberg, et al. 2019),
*testthat* (Wickham 2011), *knitr* (Xie 2020 , 2015), *rmarkdown*
(Allaire et al. 2019; Xie, Allaire, and Grolemund 2018), and the superb
guidance in the book: *R packages: organize, test, document, and share
your code*, by Wickham (2015). In addition, this package relies on a set
of external packages from the tidyverse universe, including: *dplyr*
(Wickham, François, et al. 2019), *tidyr* (Wickham and Henry 2019),
*tibble* (Müller and Wickham 2019), *stringr* (Wickham 2019), *readr*
(Wickham, Hester, and Francois 2018), *magrittr* (Bache and Wickham
2014), *ggplot2* (Wickham 2016), *rlang* (Henry and Wickham 2020), and
*purrr* (Henry and Wickham 2019) for internal functioning as well as
specialised statistics; *polyaAeppli* (Burden 2014).

## Installation

You can install the released version of point

``` r
# Install point from GitHub:
# install.packages("devtools")
devtools::install_github("point")
```

## Usage

Load point with `library`.

``` r
library(point)
```

## The point workflow

A more detailed outline of the general point workflow is given in the
vignette *IC-introduction* (`vignette("IC-introduction")`).

<img src="vignettes/workflow.png" width="100%" />

To read, process and analyse raw ion count data use the functions:

  - `read_IC`: raw ion count data
  - `cor_IC`: process ion count data
  - `stat_Xt`: analyse single ion count data
  - `stat_R`: analyse ion ratios

## Example 1: internal precision of isotope ratios

This is an example of how *Cameca NanoSIMS50L* raw data files can be
extracted, processed and analysed for the <sup>13</sup>C/<sup>12</sup>C
isotope ratio (![R](http://chart.apis.google.com/chart?cht=tx&chl=R
"R")). This produces a [tibble](https://tibble.tidyverse.org/) with
descriptive and predictive (Poisson) statistics (demarcated with an
![\\\\\\hat{\\\\\\phantom{,}}](http://chart.apis.google.com/chart?cht=tx&chl=%5C%5C%5Chat%7B%5C%5C%5Cphantom%7B%2C%7D%7D
"\\\\\\hat{\\\\\\phantom{,}}")) of the ion count data. This can be done
for single analysis in order to obtain internal precision.

``` r
library(point)
library(dplyr) # for data manipulation
library(stringr) # character string manipulation
```

``` r
# Use point_example() to access the examples bundled with this package in the
# inst/extdata directory.

# Raw data containing 13C and 12C counts on carbonate
tb.rw <- read_IC(point_example("2018-01-19-GLENDON"))

# Processing raw ion count data
tb.pr <- cor_IC(tb.rw, 
                N = N.rw, 
                t = t.rw, 
                Det = det_type.mt, 
                deadtime = 44, 
                thr_PHD = 50)

# Descriptive an predictive statistics for 13C/12C ratios
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

| sample.nm                        | file.nm                   | R                                                                                                                                                                                                                                                                                                                                                                                       | ![n](http://chart.apis.google.com/chart?cht=tx&chl=n "n") | ![\\bar{R}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cbar%7BR%7D "\\bar{R}") | ![s\_{R}](http://chart.apis.google.com/chart?cht=tx&chl=s_%7BR%7D "s_{R}") | ![\\epsilon\_{R} \\,](http://chart.apis.google.com/chart?cht=tx&chl=%5Cepsilon_%7BR%7D%20%5C%2C "\\epsilon_{R} \\,") (‰) | ![s\_{\\bar{R}}](http://chart.apis.google.com/chart?cht=tx&chl=s_%7B%5Cbar%7BR%7D%7D "s_{\\bar{R}}") | ![\\epsilon\_{\\bar{R}} \\,](http://chart.apis.google.com/chart?cht=tx&chl=%5Cepsilon_%7B%5Cbar%7BR%7D%7D%20%5C%2C "\\epsilon_{\\bar{R}} \\,") (‰) | ![\\hat{s}\_{R}](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7Bs%7D_%7BR%7D "\\hat{s}_{R}") | ![\\hat{\\epsilon}\_{R} \\,](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7B%5Cepsilon%7D_%7BR%7D%20%5C%2C "\\hat{\\epsilon}_{R} \\,") (‰) | ![\\hat{s}\_{\\bar{R}}](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7Bs%7D_%7B%5Cbar%7BR%7D%7D "\\hat{s}_{\\bar{R}}") | ![\\hat{\\epsilon}\_{\\bar{R}} \\,](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7B%5Cepsilon%7D_%7B%5Cbar%7BR%7D%7D%20%5C%2C "\\hat{\\epsilon}_{\\bar{R}} \\,") (‰) | ![\\chi^{2}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E%7B2%7D "\\chi^{2}") |
| :------------------------------- | :------------------------ | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------: | ----------------------------------------------------------------------------------: | -------------------------------------------------------------------------: | -----------------------------------------------------------------------------------------------------------------------: | ---------------------------------------------------------------------------------------------------: | -------------------------------------------------------------------------------------------------------------------------------------------------: | ---------------------------------------------------------------------------------------------------: | -------------------------------------------------------------------------------------------------------------------------------------------------: | -----------------------------------------------------------------------------------------------------------------------------: | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | ---------------------------------------------------------------------------------------: |
| Belemnite, Indium                | 2018-01-19-GLENDON\_1\_1  | ![\\phantom{,}^{13}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}")/![\\phantom{,}^{12}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}") |                                                      3900 |                                                                             0.01097 |                                                                   0.001020 |                                                                                                                     93.0 |                                                                                             1.63e-05 |                                                                                                                                               1.49 |                                                                                             0.001017 |                                                                                                                                               92.8 |                                                                                                                       1.63e-05 |                                                                                                                                                                         1.49 |                                                                                    1.004 |
| Belemnite, Indium                | 2018-01-19-GLENDON\_1\_2  | ![\\phantom{,}^{13}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}")/![\\phantom{,}^{12}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}") |                                                      3900 |                                                                             0.01097 |                                                                   0.000777 |                                                                                                                     70.8 |                                                                                             1.24e-05 |                                                                                                                                               1.13 |                                                                                             0.000769 |                                                                                                                                               70.1 |                                                                                                                       1.23e-05 |                                                                                                                                                                         1.12 |                                                                                    1.020 |
| Belemnite, Indium                | 2018-01-19-GLENDON\_1\_3  | ![\\phantom{,}^{13}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}")/![\\phantom{,}^{12}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}") |                                                      3900 |                                                                             0.01099 |                                                                   0.000731 |                                                                                                                     66.5 |                                                                                             1.17e-05 |                                                                                                                                               1.07 |                                                                                             0.000721 |                                                                                                                                               65.5 |                                                                                                                       1.15e-05 |                                                                                                                                                                         1.05 |                                                                                    1.030 |
| W12, Indium (W12\_Indium\_1\_1)  | 2018-01-19-GLENDON\_2\_1  | ![\\phantom{,}^{13}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}")/![\\phantom{,}^{12}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}") |                                                      4000 |                                                                             0.01077 |                                                                   0.000759 |                                                                                                                     70.5 |                                                                                             1.20e-05 |                                                                                                                                               1.11 |                                                                                             0.000759 |                                                                                                                                               70.5 |                                                                                                                       1.20e-05 |                                                                                                                                                                         1.11 |                                                                                    1.000 |
| W12, Indium (W12\_Indium\_1\_10) | 2018-01-19-GLENDON\_2\_10 | ![\\phantom{,}^{13}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}")/![\\phantom{,}^{12}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}") |                                                      4000 |                                                                             0.00503 |                                                                   0.318830 |                                                                                                                  63447.2 |                                                                                             5.04e-03 |                                                                                                                                            1003.19 |                                                                                             0.318614 |                                                                                                                                            63404.3 |                                                                                                                       5.04e-03 |                                                                                                                                                                      1002.51 |                                                                                    1.001 |
| W12, Indium (W12\_Indium\_1\_10) | 2018-01-19-GLENDON\_3\_1  | ![\\phantom{,}^{13}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}")/![\\phantom{,}^{12}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}") |                                                      4000 |                                                                             0.01094 |                                                                   0.000916 |                                                                                                                     83.8 |                                                                                             1.45e-05 |                                                                                                                                               1.32 |                                                                                             0.000919 |                                                                                                                                               84.0 |                                                                                                                       1.45e-05 |                                                                                                                                                                         1.33 |                                                                                    0.996 |

## Example 2: external precision of isotope ratios

To calculate the external reproducibility of isotope ratios one needs to
use the total ion counts and count rate of one analysis. The latter is
equivalent to the mean ion count rate, which can be calculated with the
function `stat_Xt`.

``` r
# Single ion descriptive an predictive statistics for all measured ions
tb.Xt <- stat_Xt(tb.pr, 
                 Xt = Xt.pr, 
                 N = N.pr, 
                 species = species.nm, 
                 sample.nm, 
                 file.nm,
                 latex = FALSE)

# For this particular study a belemnite was used as reference material. 
tb.R.ext  <- stat_R(tb.Xt, 
                    Xt = M_Xt.pr, 
                    N = Ntot_Xt.pr, 
                    species = species.nm, 
                    ion1 = "13C", 
                    ion2 = "12C",  
                    sample.nm, 
                    latex = TRUE)
```

| sample.nm         | R                                                                                                                                                                                                                                                                                                                                                                                       | ![n](http://chart.apis.google.com/chart?cht=tx&chl=n "n") | ![\\bar{R}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cbar%7BR%7D "\\bar{R}") | ![s\_{R}](http://chart.apis.google.com/chart?cht=tx&chl=s_%7BR%7D "s_{R}") | ![\\epsilon\_{R} \\,](http://chart.apis.google.com/chart?cht=tx&chl=%5Cepsilon_%7BR%7D%20%5C%2C "\\epsilon_{R} \\,") (‰) | ![s\_{\\bar{R}}](http://chart.apis.google.com/chart?cht=tx&chl=s_%7B%5Cbar%7BR%7D%7D "s_{\\bar{R}}") | ![\\epsilon\_{\\bar{R}} \\,](http://chart.apis.google.com/chart?cht=tx&chl=%5Cepsilon_%7B%5Cbar%7BR%7D%7D%20%5C%2C "\\epsilon_{\\bar{R}} \\,") (‰) | ![\\hat{s}\_{R}](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7Bs%7D_%7BR%7D "\\hat{s}_{R}") | ![\\hat{\\epsilon}\_{R} \\,](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7B%5Cepsilon%7D_%7BR%7D%20%5C%2C "\\hat{\\epsilon}_{R} \\,") (‰) | ![\\hat{s}\_{\\bar{R}}](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7Bs%7D_%7B%5Cbar%7BR%7D%7D "\\hat{s}_{\\bar{R}}") | ![\\hat{\\epsilon}\_{\\bar{R}} \\,](http://chart.apis.google.com/chart?cht=tx&chl=%5Chat%7B%5Cepsilon%7D_%7B%5Cbar%7BR%7D%7D%20%5C%2C "\\hat{\\epsilon}_{\\bar{R}} \\,") (‰) | ![\\chi^{2}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cchi%5E%7B2%7D "\\chi^{2}") |
| :---------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------: | ----------------------------------------------------------------------------------: | -------------------------------------------------------------------------: | -----------------------------------------------------------------------------------------------------------------------: | ---------------------------------------------------------------------------------------------------: | -------------------------------------------------------------------------------------------------------------------------------------------------: | ---------------------------------------------------------------------------------------------------: | -------------------------------------------------------------------------------------------------------------------------------------------------: | -----------------------------------------------------------------------------------------------------------------------------: | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | ---------------------------------------------------------------------------------------: |
| Belemnite, Indium | ![\\phantom{,}^{13}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B13%7D "\\phantom{,}^{13}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}")/![\\phantom{,}^{12}](http://chart.apis.google.com/chart?cht=tx&chl=%5Cphantom%7B%2C%7D%5E%7B12%7D "\\phantom{,}^{12}")C![\_{}](http://chart.apis.google.com/chart?cht=tx&chl=_%7B%7D "_{}") |                                                         3 |                                                                               0.011 |                                                                   1.55e-05 |                                                                                                                     1.41 |                                                                                              8.9e-06 |                                                                                                                                              0.814 |                                                                                              1.3e-05 |                                                                                                                                               1.18 |                                                                                                                        7.5e-06 |                                                                                                                                                                        0.681 |                                                                                     1.43 |

For more detailed information:

*IC-read* (`vignette("IC-read")`): reading raw ion count data use  
*IC-process* (`vignette("IC-process")`): processing ion count data  
*IC-precision* (`vignette("IC-precision")`): statistics concerning ion
count precision  
*IC-diagnostics* (`vignette("IC-diagnostics")`): diagnostics on internal
variation

# References

<div id="refs" class="references">

<div id="ref-rmarkdown1">

Allaire, JJ, Yihui Xie, Jonathan McPherson, Javier Luraschi, Kevin
Ushey, Aron Atkins, Hadley Wickham, Joe Cheng, Winston Chang, and
Richard Iannone. 2019. *Rmarkdown: Dynamic Documents for R*.
<https://github.com/rstudio/rmarkdown>.

</div>

<div id="ref-magrittr">

Bache, Stefan Milton, and Hadley Wickham. 2014. *Magrittr: A
Forward-Pipe Operator for R*.
<https://CRAN.R-project.org/package=magrittr>.

</div>

<div id="ref-polyaAeppli">

Burden, Conrad. 2014. *PolyaAeppli: Implementation of the Polya-Aeppli
Distribution*. <https://CRAN.R-project.org/package=polyaAeppli>.

</div>

<div id="ref-purrr">

Henry, Lionel, and Hadley Wickham. 2019. *Purrr: Functional Programming
Tools*. <https://CRAN.R-project.org/package=purrr>.

</div>

<div id="ref-rlang">

———. 2020. *Rlang: Functions for Base Types and Core R and ’Tidyverse’
Features*. <https://CRAN.R-project.org/package=rlang>.

</div>

<div id="ref-tibble">

Müller, Kirill, and Hadley Wickham. 2019. *Tibble: Simple Data Frames*.
<https://CRAN.R-project.org/package=tibble>.

</div>

<div id="ref-rversion">

R Core Team. 2020. *R: A Language and Environment for Statistical
Computing*. Vienna, Austria: R Foundation for Statistical Computing.
<https://www.R-project.org/>.

</div>

<div id="ref-testthat">

Wickham, Hadley. 2011. “Testthat: Get Started with Testing.” *The R
Journal* 3: 5–10.
<https://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf>.

</div>

<div id="ref-Wickham2015">

———. 2015. *R Packages: Organize, Test, Document, and Share Your Code*.
O’Reilly Media, Inc.

</div>

<div id="ref-ggplot2">

———. 2016. *Ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag New York. <https://ggplot2.tidyverse.org>.

</div>

<div id="ref-stringr">

———. 2019. *Stringr: Simple, Consistent Wrappers for Common String
Operations*. <https://CRAN.R-project.org/package=stringr>.

</div>

<div id="ref-roxygen2">

Wickham, Hadley, Peter Danenberg, Gábor Csárdi, and Manuel Eugster.
2019. *Roxygen2: In-Line Documentation for R*.
<https://CRAN.R-project.org/package=roxygen2>.

</div>

<div id="ref-dplyr">

Wickham, Hadley, Romain François, Lionel Henry, and Kirill Müller. 2019.
*Dplyr: A Grammar of Data Manipulation*.
<https://CRAN.R-project.org/package=dplyr>.

</div>

<div id="ref-tidyr">

Wickham, Hadley, and Lionel Henry. 2019. *Tidyr: Tidy Messy Data*.
<https://CRAN.R-project.org/package=tidyr>.

</div>

<div id="ref-devtools">

Wickham, Hadley, Jim Hester, and Winston Chang. 2019. *Devtools: Tools
to Make Developing R Packages Easier*.
<https://CRAN.R-project.org/package=devtools>.

</div>

<div id="ref-readr">

Wickham, Hadley, Jim Hester, and Romain Francois. 2018. *Readr: Read
Rectangular Text Data*. <https://CRAN.R-project.org/package=readr>.

</div>

<div id="ref-knitr2">

Xie, Yihui. 2015. *Dynamic Documents with R and Knitr*. 2nd ed. Boca
Raton, Florida: Chapman; Hall/CRC. <https://yihui.org/knitr/>.

</div>

<div id="ref-knitr1">

———. 2020. *Knitr: A General-Purpose Package for Dynamic Report
Generation in R*. <https://yihui.org/knitr/>.

</div>

<div id="ref-rmarkdown2">

Xie, Yihui, J.J. Allaire, and Garrett Grolemund. 2018. *R Markdown: The
Definitive Guide*. Boca Raton, Florida: Chapman; Hall/CRC.
<https://bookdown.org/yihui/rmarkdown>.

</div>

</div>
