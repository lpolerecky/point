
<!--  use the --webtex argument in the YAML to render equations -->

# Introduction to point

<!-- badges: start -->

<!-- badges: end -->

This projects was originally inspired by the lack of detailed insight in
the inner workings of the default software for the Cameca NanoSIMS50L
(Utrecht University). Hence this project has the objective of processing
raw count data into ion and isotope ratios of point-sourced
measurements; and to establish the internal and external precision of,
respectively, individual measurements and complete series of
measurements. Access to raw ion count data is useful as it allows
detection of anomolous values associated with e.g. machine instability
or heterogeneity of the analysed sample. Upon detection, anomolous
values can be omitted or further analysed to delineate the source of
variation.

## Installation

You can install the released version of point

``` r
# Install point rom GitHub:
# install.packages("devtools")
devtools::install_github("point")
```

## Usage

Load point with library

``` r
library(point)
```

To read, process and analyse raw ion coun data use the functions:

  - `read_IC`: raw ion count data
  - `cor_IC`: process ion count data
  - `stat_Xt`: analyse single ion count data
  - `stat_R`: analyse ion ratios

## Example 1: internal precision of isotope ratios

This is an example of how Cameca NanoSIMS50L raw datafiles can be
extracted, processed and analysed for the <sup>13</sup>C/<sup>12</sup>C
isotope ratio (\(R\)). This produces a tibble with descriptive and
predictive poisson statistics (demarcated with an
\(\\\hat{\\\phantom{,}}\) ) of the ion count data. This can be done for
single count blocks in order to obtain internal
precision.

``` r
# Use point_example() to access the examples bundled with this package in the
# inst/extdata directory. The examples directories are named:
# 2020-01-17-TREASURE and "2018-01-19-GLENDON"

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
               latex = TRUE, 
               output = "sum")

knitr::kable(head(tb.R),
             format.args = list(digits = 3, 
                                format = "G", 
                                flag = "0")) 
```

| sample.nm                        | file.nm                  | R                                                         | \(n\) | \(\bar{R}\) | \(s_{R}\) | \(\epsilon_{R} \,\) (‰) | \(s_{\bar{R}}\) | \(\epsilon_{\bar{R}} \,\) (‰) | \(\hat{s}_{R}\) | \(\hat{\epsilon}_{R} \,\) (‰) | \(\hat{s}_{\bar{R}}\) | \(\hat{\epsilon}_{\bar{R}} \,\) (‰) | \(\chi^{2}\) |
| :------------------------------- | :----------------------- | :-------------------------------------------------------- | ----: | ----------: | --------: | ----------------------: | --------------: | ----------------------------: | --------------: | ----------------------------: | --------------------: | ----------------------------------: | -----------: |
| Belemnite, Indium                | 2018-01-19-GLENDON\_1\_1 | \(\phantom{,}^{13}\)C\(_{}\)/\(\phantom{,}^{12}\)C\(_{}\) |  3900 |      0.0110 |  0.001020 |                    93.0 |        1.63e-05 |                          1.49 |        0.001018 |                          92.8 |              1.63e-05 |                                1.49 |        1.004 |
| Belemnite, Indium                | 2018-01-19-GLENDON\_1\_2 | \(\phantom{,}^{13}\)C\(_{}\)/\(\phantom{,}^{12}\)C\(_{}\) |  3900 |      0.0110 |  0.000777 |                    70.8 |        1.24e-05 |                          1.13 |        0.000770 |                          70.1 |              1.23e-05 |                                1.12 |        1.020 |
| Belemnite, Indium                | 2018-01-19-GLENDON\_1\_3 | \(\phantom{,}^{13}\)C\(_{}\)/\(\phantom{,}^{12}\)C\(_{}\) |  3900 |      0.0110 |  0.000732 |                    66.5 |        1.17e-05 |                          1.07 |        0.000721 |                          65.5 |              1.15e-05 |                                1.05 |        1.030 |
| W12, Indium (W12\_Indium\_1\_1)  | 2018-01-19-GLENDON\_2\_1 | \(\phantom{,}^{13}\)C\(_{}\)/\(\phantom{,}^{12}\)C\(_{}\) |  4000 |      0.0108 |  0.000760 |                    70.5 |        1.20e-05 |                          1.11 |        0.000760 |                          70.5 |              1.20e-05 |                                1.11 |        1.000 |
| W12, Indium (W12\_Indium\_1\_10) | 2018-01-19-GLENDON\_3\_1 | \(\phantom{,}^{13}\)C\(_{}\)/\(\phantom{,}^{12}\)C\(_{}\) |  4000 |      0.0109 |  0.000917 |                    83.8 |        1.45e-05 |                          1.32 |        0.000919 |                          84.0 |              1.45e-05 |                                1.33 |        0.995 |
| W12, Indium (W12\_Indium\_1\_11) | 2018-01-19-GLENDON\_3\_2 | \(\phantom{,}^{13}\)C\(_{}\)/\(\phantom{,}^{12}\)C\(_{}\) |  4000 |      0.0108 |  0.000732 |                    68.1 |        1.16e-05 |                          1.08 |        0.000723 |                          67.2 |              1.14e-05 |                                1.06 |        1.025 |

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
                 latex = FALSE,
                 output = "sum")

# For this particular run a belemnite was used as reference material. 
tb.R.ext  <- stat_R(tb.Xt, 
                    Xt = M_Xt.pr, 
                    N = Ntot_Xt.pr, 
                    species = species.nm, 
                    ion1 = "13C", 
                    ion2 = "12C",  
                    sample.nm, 
                    latex = TRUE,
                    output = "sum")

knitr::kable(head(tb.R.ext %>% 
                    filter(str_detect(sample.nm, "Belemnite"))),
             format.args = list(digits = 3, 
                                format = "G", 
                                flag = "0")) 
```

| sample.nm         | R                                                         | \(n\) | \(\bar{R}\) | \(s_{R}\) | \(\epsilon_{R} \,\) (‰) | \(s_{\bar{R}}\) | \(\epsilon_{\bar{R}} \,\) (‰) | \(\hat{s}_{R}\) | \(\hat{\epsilon}_{R} \,\) (‰) | \(\hat{s}_{\bar{R}}\) | \(\hat{\epsilon}_{\bar{R}} \,\) (‰) | \(\chi^{2}\) |
| :---------------- | :-------------------------------------------------------- | ----: | ----------: | --------: | ----------------------: | --------------: | ----------------------------: | --------------: | ----------------------------: | --------------------: | ----------------------------------: | -----------: |
| Belemnite, Indium | \(\phantom{,}^{13}\)C\(_{}\)/\(\phantom{,}^{12}\)C\(_{}\) |     3 |       0.011 |  1.71e-05 |                    1.56 |         9.9e-06 |                         0.898 |         1.3e-05 |                          1.18 |               7.5e-06 |                               0.681 |         1.74 |

For more detailed information:

  - `vignette("IC-read")`: reading raw ion count data use
  - `vignette("IC-process")`: processing ion count data  
  - `vignette("IC-precision")`: statistics concerning ion count
    precision
