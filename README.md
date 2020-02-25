
<!--  use the --webtex argument in the YAML to render equations -->

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

This is an example of how Cameca NanoSIMS50L raw datafiles can be
extracted, processed and analysed for the
![^{13}](https://latex.codecogs.com/png.latex?%5E%7B13%7D
"^{13}")C/![^{12}](https://latex.codecogs.com/png.latex?%5E%7B12%7D
"^{12}")C isotope ratio (R). This produces a tibble with descriptive and
predictive poisson statistics (demarcated with
![\\hat](https://latex.codecogs.com/png.latex?%5Chat "\\hat")) of the
ion count data.

``` r

# library(point)
devtools::load_all(".")
#> Loading point
# Use system.file() to access the examples bundled with this package in the
# inst/extdata directory. The examples directories are named:
# 2020-01-17-TREASURE and "2018-01-19-GLENDON"

# raw data containing 13C and 12C counts on carbonate
tb.rw <- read_IC(system.file("extdata", "2018-01-19-GLENDON", package = "point"))

# processing raw ion count data
tb.pr <- cor_IC(tb.rw, N.rw, t.rw, det_type.mt, deadtime = 44, disc.value = NULL)

# single ion descriptive an predictive statistics for all measured ions
tb.Xt <- stat_Xt(tb.pr, Xt.pr, N.pr, file.nm, species.nm)

# descriptive an predictive statistics for 13C/12C ratios
tb.R <- stat_R(tb.pr, Xt.pr, N.pr, ID = "ID", ion1 = "13C", ion2 = "12C",
               file.nm, species.nm, latex = TRUE)


# use the second element of the list to get nice variable names rendering in 
# Rmarkdown
knitr::kable(tb.R[[1]] , escape = FALSE, 
             col.names = c("file","R", names(tb.R[[2]])),
             format.args = list(digits = 5, format = "G", flag = "0")) 
```

| file                             | R       | ![n](https://latex.codecogs.com/png.latex?n "n") | ![\\bar{x}](https://latex.codecogs.com/png.latex?%5Cbar%7Bx%7D "\\bar{x}") | ![r](https://latex.codecogs.com/png.latex?r "r") | ![s\_x](https://latex.codecogs.com/png.latex?s_x "s_x") | ![\\epsilon\_x](https://latex.codecogs.com/png.latex?%5Cepsilon_x "\\epsilon_x") | ![s\_\\bar{x}](https://latex.codecogs.com/png.latex?s_%5Cbar%7Bx%7D "s_\\bar{x}") | ![\\epsilon\_\\bar{x}](https://latex.codecogs.com/png.latex?%5Cepsilon_%5Cbar%7Bx%7D "\\epsilon_\\bar{x}") | ![\\hat{s}\_x](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_x "\\hat{s}_x") | ![\\hat{\\epsilon}\_x](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cepsilon%7D_x "\\hat{\\epsilon}_x") | ![\\hat{s}\_\\bar{x}](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_%5Cbar%7Bx%7D "\\hat{s}_\\bar{x}") | ![\\hat{\\epsilon}\_\\bar{x}](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cepsilon%7D_%5Cbar%7Bx%7D "\\hat{\\epsilon}_\\bar{x}") | ![\\chi^2](https://latex.codecogs.com/png.latex?%5Cchi%5E2 "\\chi^2") |
| :------------------------------- | :------ | -----------------------------------------------: | -------------------------------------------------------------------------: | -----------------------------------------------: | ------------------------------------------------------: | -------------------------------------------------------------------------------: | --------------------------------------------------------------------------------: | ---------------------------------------------------------------------------------------------------------: | --------------------------------------------------------------------------------: | ---------------------------------------------------------------------------------------------------------: | ----------------------------------------------------------------------------------------------------------: | -----------------------------------------------------------------------------------------------------------------------------------: | --------------------------------------------------------------------: |
| 2018-01-19-GLENDON\_1\_1.is\_txt | 13C/12C |                                             3900 |                                                                   0.010972 |                                          0.67778 |                                               0.0010201 |                                                                           92.969 |                                                                          1.63e-05 |                                                                                                     1.4887 |                                                                         0.0010179 |                                                                                                     92.769 |                                                                                                    1.63e-05 |                                                                                                                               1.4855 |                                                               1.00432 |
| 2018-01-19-GLENDON\_1\_2.is\_txt | 13C/12C |                                             3900 |                                                                   0.010978 |                                          0.72679 |                                               0.0007774 |                                                                           70.814 |                                                                          1.24e-05 |                                                                                                     1.1339 |                                                                         0.0007697 |                                                                                                     70.110 |                                                                                                    1.23e-05 |                                                                                                                               1.1227 |                                                               1.02017 |
| 2018-01-19-GLENDON\_1\_3.is\_txt | 13C/12C |                                             3900 |                                                                   0.011003 |                                          0.56653 |                                               0.0007318 |                                                                           66.511 |                                                                          1.17e-05 |                                                                                                     1.0650 |                                                                         0.0007212 |                                                                                                     65.545 |                                                                                                    1.15e-05 |                                                                                                                               1.0496 |                                                               1.02971 |
| 2018-01-19-GLENDON\_2\_1.is\_txt | 13C/12C |                                             4000 |                                                                   0.010779 |                                          0.96154 |                                               0.0007599 |                                                                           70.499 |                                                                          1.20e-05 |                                                                                                     1.1147 |                                                                         0.0007601 |                                                                                                     70.513 |                                                                                                    1.20e-05 |                                                                                                                               1.1149 |                                                               0.99959 |
| 2018-01-19-GLENDON\_2\_2.is\_txt | 13C/12C |                                             4000 |                                                                   0.010781 |                                          0.95537 |                                               0.0007530 |                                                                           69.849 |                                                                          1.19e-05 |                                                                                                     1.1044 |                                                                         0.0007388 |                                                                                                     68.529 |                                                                                                    1.17e-05 |                                                                                                                               1.0835 |                                                               1.03888 |
| 2018-01-19-GLENDON\_2\_3.is\_txt | 13C/12C |                                             4000 |                                                                   0.010749 |                                          0.96237 |                                               0.0007365 |                                                                           68.522 |                                                                          1.16e-05 |                                                                                                     1.0834 |                                                                         0.0007327 |                                                                                                     68.163 |                                                                                                    1.16e-05 |                                                                                                                               1.0777 |                                                               1.01057 |
| 2018-01-19-GLENDON\_2\_4.is\_txt | 13C/12C |                                             4000 |                                                                   0.010769 |                                          0.93776 |                                               0.0007194 |                                                                           66.798 |                                                                          1.14e-05 |                                                                                                     1.0562 |                                                                         0.0007248 |                                                                                                     67.306 |                                                                                                    1.15e-05 |                                                                                                                               1.0642 |                                                               0.98499 |
| 2018-01-19-GLENDON\_2\_5.is\_txt | 13C/12C |                                             4000 |                                                                   0.010891 |                                          0.97156 |                                               0.0007725 |                                                                           70.936 |                                                                          1.22e-05 |                                                                                                     1.1216 |                                                                         0.0007746 |                                                                                                     71.126 |                                                                                                    1.22e-05 |                                                                                                                               1.1246 |                                                               0.99468 |
| 2018-01-19-GLENDON\_2\_6.is\_txt | 13C/12C |                                             4000 |                                                                   0.010916 |                                          0.99409 |                                               0.0007044 |                                                                           64.526 |                                                                          1.11e-05 |                                                                                                     1.0202 |                                                                         0.0007134 |                                                                                                     65.357 |                                                                                                    1.13e-05 |                                                                                                                               1.0334 |                                                               0.97472 |
| 2018-01-19-GLENDON\_2\_7.is\_txt | 13C/12C |                                             4000 |                                                                   0.011088 |                                          0.97140 |                                               0.0008162 |                                                                           73.614 |                                                                          1.29e-05 |                                                                                                     1.1639 |                                                                         0.0008101 |                                                                                                     73.065 |                                                                                                    1.28e-05 |                                                                                                                               1.1553 |                                                               1.01509 |
| 2018-01-19-GLENDON\_2\_8.is\_txt | 13C/12C |                                             4000 |                                                                   0.011090 |                                          0.95804 |                                               0.0007566 |                                                                           68.224 |                                                                          1.20e-05 |                                                                                                     1.0787 |                                                                         0.0007437 |                                                                                                     67.061 |                                                                                                    1.18e-05 |                                                                                                                               1.0603 |                                                               1.03497 |
| 2018-01-19-GLENDON\_2\_9.is\_txt | 13C/12C |                                             4000 |                                                                   0.011082 |                                          0.96426 |                                               0.0007135 |                                                                           64.388 |                                                                          1.13e-05 |                                                                                                     1.0181 |                                                                         0.0007174 |                                                                                                     64.734 |                                                                                                    1.13e-05 |                                                                                                                               1.0235 |                                                               0.98935 |
| 2018-01-19-GLENDON\_3\_1.is\_txt | 13C/12C |                                             4000 |                                                                   0.010947 |                                          0.95873 |                                               0.0009169 |                                                                           83.757 |                                                                          1.45e-05 |                                                                                                     1.3243 |                                                                         0.0009190 |                                                                                                     83.957 |                                                                                                    1.45e-05 |                                                                                                                               1.3275 |                                                               0.99525 |
| 2018-01-19-GLENDON\_3\_2.is\_txt | 13C/12C |                                             4000 |                                                                   0.010755 |                                          0.93340 |                                               0.0007319 |                                                                           68.052 |                                                                          1.16e-05 |                                                                                                     1.0760 |                                                                         0.0007229 |                                                                                                     67.216 |                                                                                                    1.14e-05 |                                                                                                                               1.0628 |                                                               1.02503 |
| 2018-01-19-GLENDON\_3\_3.is\_txt | 13C/12C |                                             4000 |                                                                   0.010747 |                                          0.98815 |                                               0.0008391 |                                                                           78.075 |                                                                          1.33e-05 |                                                                                                     1.2345 |                                                                         0.0008246 |                                                                                                     76.729 |                                                                                                    1.30e-05 |                                                                                                                               1.2132 |                                                               1.03538 |

# Internal precision of the SIMS analysis

To check for internal consistency in the analytical output of the
*NanoSIMS Cameca 50L*, descriptive statistics and statistical inferences
are compared with the results of predictive statistics of the raw counts
(![N\_{(i)}](https://latex.codecogs.com/png.latex?N_%7B%28i%29%7D
"N_{(i)}")) for single ions and isotopic ratios.

## Descriptive statistics

As a first step, the counts of a single count cycle
(![N\_i](https://latex.codecogs.com/png.latex?N_i "N_i")) are normalised
against the time it took to complete the cycle
(![0.541](https://latex.codecogs.com/png.latex?0.541 "0.541") s) to
account for differences in the count times for two different isotopes
during stable isotopic SIMS analysis. Hence, for the time period
(![t](https://latex.codecogs.com/png.latex?t "t")) over which an isotope
species ![a](https://latex.codecogs.com/png.latex?a "a") during
measurement ![i](https://latex.codecogs.com/png.latex?i "i")
accumulated, the count rate is given by

  
![X\_i^{a} = N\_i^{a} /
t\_i^{a}](https://latex.codecogs.com/png.latex?X_i%5E%7Ba%7D%20%3D%20N_i%5E%7Ba%7D%20%2F%20t_i%5E%7Ba%7D
"X_i^{a} = N_i^{a} / t_i^{a}")  
{\#eq:N.rate}

The sample mean
(![\\bar{x}](https://latex.codecogs.com/png.latex?%5Cbar%7Bx%7D
"\\bar{x}")) of isotope species
![a](https://latex.codecogs.com/png.latex?a "a") over a complete count
block is given by:

  
![ \\bar{x}\_a = \\frac{1}{n} \\sum\_{i=1}^{n} X\_i^a
](https://latex.codecogs.com/png.latex?%20%5Cbar%7Bx%7D_a%20%3D%20%20%5Cfrac%7B1%7D%7Bn%7D%20%5Csum_%7Bi%3D1%7D%5E%7Bn%7D%20X_i%5Ea%20
" \\bar{x}_a =  \\frac{1}{n} \\sum_{i=1}^{n} X_i^a ")  
{\#eq:M.rate}

To validate the internal consistency of the SIMS data, it is necessary
to define the internal precision of the count blocks. This can be done
with the standard deviation
(![s\_x](https://latex.codecogs.com/png.latex?s_x "s_x")), which gives
the spread of the sample, and the standard error of the mean
(![s\_{\\bar{x}}](https://latex.codecogs.com/png.latex?s_%7B%5Cbar%7Bx%7D%7D
"s_{\\bar{x}}")), which defines how well this sample mean approximates
the true population mean
(![\\mu](https://latex.codecogs.com/png.latex?%5Cmu "\\mu")). These
statistics rely on the assumption that the underlying probability
distribution follows a normal distribution.

### Standard deviation

The standard deviation for a limited sample of the population gives a
measure of how individual measurements are spread about the mean in one
count block, and is given by

  
![ s\_x = \\sqrt{\\sum\_{i=1}^{n} \\frac{(x\_{i}-\\bar{x})^2}{n-1}}
](https://latex.codecogs.com/png.latex?%20s_x%20%3D%20%5Csqrt%7B%5Csum_%7Bi%3D1%7D%5E%7Bn%7D%20%20%5Cfrac%7B%28x_%7Bi%7D-%5Cbar%7Bx%7D%29%5E2%7D%7Bn-1%7D%7D%20
" s_x = \\sqrt{\\sum_{i=1}^{n}  \\frac{(x_{i}-\\bar{x})^2}{n-1}} ")  
{\#eq:std}

where ![n](https://latex.codecogs.com/png.latex?n "n") is the number of
measurement cycles in the count block and
![x\_i](https://latex.codecogs.com/png.latex?x_i "x_i") is the
![i](https://latex.codecogs.com/png.latex?i "i")-th measurement cycle.
The number of measurements is subtracted with one (![n
- 1](https://latex.codecogs.com/png.latex?n%20-%201 "n - 1")) to show
that only ![n - 1](https://latex.codecogs.com/png.latex?n%20-%201
"n - 1") of the
![(x\_{i}-\\bar{x})^2](https://latex.codecogs.com/png.latex?%28x_%7Bi%7D-%5Cbar%7Bx%7D%29%5E2
"(x_{i}-\\bar{x})^2") are independent. The sample standard can inform
about the confidence whether a single measurements falls within a given
range of the sample mean value. One can also report the standard
deviation as a number relative to the sample mean by applying the
following transformation

  
![ e\_x = \\frac{s\_x}{\\bar{x}}
](https://latex.codecogs.com/png.latex?%20e_x%20%3D%20%5Cfrac%7Bs_x%7D%7B%5Cbar%7Bx%7D%7D%20
" e_x = \\frac{s_x}{\\bar{x}} ")  
{\#eq:rel.std}

### Standard error of the mean

The standard error of the mean
(![s\_{\\bar{x}}](https://latex.codecogs.com/png.latex?s_%7B%5Cbar%7Bx%7D%7D
"s_{\\bar{x}}")) provides a measure of how well the mean of a limited
sample (i.e., count block) approximates the actual population mean. This
measure can be used to gauge the precision of the count block with
![n](https://latex.codecogs.com/png.latex?n "n") measurement cycles.
This value is dependent on the number of measurements
(![n](https://latex.codecogs.com/png.latex?n "n")) and thus becomes
smaller with increasing measurement numbers (i.e.
![\\bar{x}](https://latex.codecogs.com/png.latex?%5Cbar%7Bx%7D
"\\bar{x}") becomes more precise). The standard error of the mean is
given by the following equation

  
![ s\_{\\bar{x}} = \\frac{s\_{x}}{\\sqrt{n}}
](https://latex.codecogs.com/png.latex?%20s_%7B%5Cbar%7Bx%7D%7D%20%3D%20%5Cfrac%7Bs_%7Bx%7D%7D%7B%5Csqrt%7Bn%7D%7D%20
" s_{\\bar{x}} = \\frac{s_{x}}{\\sqrt{n}} ")  
{\#eq:se}

This equation as well can be recast to obtain the relative standard
error of the mean
(![e\_{\\bar{x}}](https://latex.codecogs.com/png.latex?e_%7B%5Cbar%7Bx%7D%7D
"e_{\\bar{x}}")), where

  
![ e\_{\\bar{x}} = \\frac{s\_{\\bar{x}}}{\\bar{x}}
](https://latex.codecogs.com/png.latex?%20e_%7B%5Cbar%7Bx%7D%7D%20%3D%20%5Cfrac%7Bs_%7B%5Cbar%7Bx%7D%7D%7D%7B%5Cbar%7Bx%7D%7D%20
" e_{\\bar{x}} = \\frac{s_{\\bar{x}}}{\\bar{x}} ")  
{\#eq:rel.se}

### Error propagation for descriptive statistics of isotope ratios

The
![\\bar{x}\_{R}](https://latex.codecogs.com/png.latex?%5Cbar%7Bx%7D_%7BR%7D
"\\bar{x}_{R}") can be calculated from the mean values of the specific
ions of the complete count block.

  
![ \\bar{x}\_{R} = \\frac{\\frac{1}{n}\\sum\_{i = 1}^{n}
X\_i^{b}}{\\frac{1}{n}\\sum\_{i = 1}^{n} X\_i^{a}}
](https://latex.codecogs.com/png.latex?%20%5Cbar%7Bx%7D_%7BR%7D%20%3D%20%5Cfrac%7B%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7Bi%20%3D%201%7D%5E%7Bn%7D%20X_i%5E%7Bb%7D%7D%7B%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7Bi%20%3D%201%7D%5E%7Bn%7D%20X_i%5E%7Ba%7D%7D%20
" \\bar{x}_{R} = \\frac{\\frac{1}{n}\\sum_{i = 1}^{n} X_i^{b}}{\\frac{1}{n}\\sum_{i = 1}^{n} X_i^{a}} ")  
{\#eq:mean.R}

and this value can be considered as an estimate of the true isotopic
value (![\\mu\_R](https://latex.codecogs.com/png.latex?%5Cmu_R
"\\mu_R")). The uncertainties associated with the SIMS count rates of
the individual variables
![X^{b}](https://latex.codecogs.com/png.latex?X%5E%7Bb%7D "X^{b}") (e.g.
<sup>13</sup>C) and
![X^{a}](https://latex.codecogs.com/png.latex?X%5E%7Ba%7D "X^{a}") (e.g.
<sup>12</sup>C) need to be combined. This can be achieved by applying;
*The formula for exact propagation of error*

  
![s\_x^{2} \\approx \\sum\_{i = 1}^{n} \\left\[ \\left( \\frac{\\partial
F}{\\partial z\_i} \\right) s\_i^{2} \\right\] + 2 \\sum\_{j = 1}^{n}
\\sum\_{k = 1}^{n} \\left\[ \\left( \\frac{\\partial F}{\\partial z\_j}
\\right) \\left( \\frac{\\partial F}{\\partial z\_k} \\right) s\_j s\_k
\\right\]](https://latex.codecogs.com/png.latex?s_x%5E%7B2%7D%20%5Capprox%20%5Csum_%7Bi%20%3D%201%7D%5E%7Bn%7D%20%5Cleft%5B%20%5Cleft%28%20%5Cfrac%7B%5Cpartial%20F%7D%7B%5Cpartial%20z_i%7D%20%5Cright%29%20s_i%5E%7B2%7D%20%5Cright%5D%20%2B%202%20%5Csum_%7Bj%20%3D%201%7D%5E%7Bn%7D%20%5Csum_%7Bk%20%3D%201%7D%5E%7Bn%7D%20%5Cleft%5B%20%5Cleft%28%20%5Cfrac%7B%5Cpartial%20F%7D%7B%5Cpartial%20z_j%7D%20%5Cright%29%20%5Cleft%28%20%5Cfrac%7B%5Cpartial%20F%7D%7B%5Cpartial%20z_k%7D%20%5Cright%29%20s_j%20s_k%20%5Cright%5D
"s_x^{2} \\approx \\sum_{i = 1}^{n} \\left[ \\left( \\frac{\\partial F}{\\partial z_i} \\right) s_i^{2} \\right] + 2 \\sum_{j = 1}^{n} \\sum_{k = 1}^{n} \\left[ \\left( \\frac{\\partial F}{\\partial z_j} \\right) \\left( \\frac{\\partial F}{\\partial z_k} \\right) s_j s_k \\right]")  
{\#eq:er.prop}

, which ensures proper propagation of the error. In this formulation
![r\_{jk}](https://latex.codecogs.com/png.latex?r_%7Bjk%7D "r_{jk}")
stands for the correlation coefficient for the variables
![z\_j](https://latex.codecogs.com/png.latex?z_j "z_j") and
![z\_k](https://latex.codecogs.com/png.latex?z_k "z_k"), as defined by

  
![ r\_{jk} = \\frac{1}{\\left(n-1\\right) s\_j s\_k} \\sum\_{i=1}^n{
\\left\[ \\left(z\_{j}\\right)\_i - \\bar{z}\_j \\right\] \\left\[
\\left(z\_{k}\\right)\_i - \\bar{z}\_k \\right\]}
](https://latex.codecogs.com/png.latex?%20r_%7Bjk%7D%20%3D%20%5Cfrac%7B1%7D%7B%5Cleft%28n-1%5Cright%29%20s_j%20s_k%7D%20%5Csum_%7Bi%3D1%7D%5En%7B%20%5Cleft%5B%20%5Cleft%28z_%7Bj%7D%5Cright%29_i%20-%20%5Cbar%7Bz%7D_j%20%5Cright%5D%20%5Cleft%5B%20%5Cleft%28z_%7Bk%7D%5Cright%29_i%20-%20%5Cbar%7Bz%7D_k%20%5Cright%5D%7D%20
" r_{jk} = \\frac{1}{\\left(n-1\\right) s_j s_k} \\sum_{i=1}^n{ \\left[ \\left(z_{j}\\right)_i - \\bar{z}_j \\right] \\left[ \\left(z_{k}\\right)_i - \\bar{z}_k \\right]} ")  
{\#eq:corr}

yields an estimate for the sample correlation coefficient, where values
can range between ![-1](https://latex.codecogs.com/png.latex?-1 "-1")
and ![+1](https://latex.codecogs.com/png.latex?%2B1 "+1"), and thereby
recording a inverse or positive linear correlation between the
variables, and no correlation if
![r](https://latex.codecogs.com/png.latex?r "r") falls close to zero.
For this the `R base` function `cor()` was used, with the `method`
argument set to `"pearson"`.

Recasting \[@eq:er.prop\] for when
![\\mathrm{F}(...)](https://latex.codecogs.com/png.latex?%5Cmathrm%7BF%7D%28...%29
"\\mathrm{F}(...)") is ![R](https://latex.codecogs.com/png.latex?R "R"),
and with the variables
![\\bar{x}^{b}](https://latex.codecogs.com/png.latex?%5Cbar%7Bx%7D%5E%7Bb%7D
"\\bar{x}^{b}") (e.g. <sup>13</sup>C) and
![\\bar{x}^{a}](https://latex.codecogs.com/png.latex?%5Cbar%7Bx%7D%5E%7Ba%7D
"\\bar{x}^{a}") (e.g. <sup>12</sup>C), yields the following equation:

  
![ s\_{x\_{R}} = \\sqrt{ \\left( \\frac{ s\_{x\_{b}}}{\\bar{x}\_{b}}
\\right)^2 + \\left( \\frac{ s\_{x\_{a}}}{\\bar{x}\_{a}} \\right)^2 -
\\frac{2r\_{X^{b}X^{a}} s\_{x\_{b}}
s\_{x\_{a}}}{\\bar{x}\_{b}\\bar{x}\_{a}}} \\times \\bar{x}\_{R}
](https://latex.codecogs.com/png.latex?%20%20s_%7Bx_%7BR%7D%7D%20%3D%20%5Csqrt%7B%20%5Cleft%28%20%5Cfrac%7B%20s_%7Bx_%7Bb%7D%7D%7D%7B%5Cbar%7Bx%7D_%7Bb%7D%7D%20%5Cright%29%5E2%20%2B%20%5Cleft%28%20%5Cfrac%7B%20s_%7Bx_%7Ba%7D%7D%7D%7B%5Cbar%7Bx%7D_%7Ba%7D%7D%20%5Cright%29%5E2%20-%20%5Cfrac%7B2r_%7BX%5E%7Bb%7DX%5E%7Ba%7D%7D%20s_%7Bx_%7Bb%7D%7D%20s_%7Bx_%7Ba%7D%7D%7D%7B%5Cbar%7Bx%7D_%7Bb%7D%5Cbar%7Bx%7D_%7Ba%7D%7D%7D%20%5Ctimes%20%5Cbar%7Bx%7D_%7BR%7D%20
"  s_{x_{R}} = \\sqrt{ \\left( \\frac{ s_{x_{b}}}{\\bar{x}_{b}} \\right)^2 + \\left( \\frac{ s_{x_{a}}}{\\bar{x}_{a}} \\right)^2 - \\frac{2r_{X^{b}X^{a}} s_{x_{b}} s_{x_{a}}}{\\bar{x}_{b}\\bar{x}_{a}}} \\times \\bar{x}_{R} ")  
{\#eq:er.prop.ad}

The standard error of the mean isotope value
![s\_{\\bar{x}\_{R}}](https://latex.codecogs.com/png.latex?s_%7B%5Cbar%7Bx%7D_%7BR%7D%7D
"s_{\\bar{x}_{R}}") is obtained through diving
![s\_{x\_{R}}](https://latex.codecogs.com/png.latex?s_%7Bx_%7BR%7D%7D
"s_{x_{R}}") by
![\\sqrt(n)](https://latex.codecogs.com/png.latex?%5Csqrt%28n%29
"\\sqrt(n)"). Both the standard deviation and standard error of the mean
of the isotope value can be expressed as relative values in ‰ by
dividing them through the
![\\bar{x}\_{R}](https://latex.codecogs.com/png.latex?%5Cbar%7Bx%7D_%7BR%7D
"\\bar{x}_{R}") and multiplying by
![1.000](https://latex.codecogs.com/png.latex?1.000 "1.000").

## Predictive statistics (Poisson)

SIMS measurements have an inherent fundamental imprecision, which is
dictated by the random nature of secondary ion production. This restrict
the precision of the analysis to a certain analytical threshold. The
amplitude of this inherent variation can be gauged with Poisson
statistics. The Poisson distribution describes the likelihood of random
events occurring over a defined (and fixed) time-period. Further
conditions to be satisfied to validate the assumption of a Poisson
distribution is the observation that
![N](https://latex.codecogs.com/png.latex?N "N") should be able to occur
over a larger number of occasions and that the probability of the event
occurring at a particular occasions is limited but constant. In the case
of SIMS measurements ![N\_i](https://latex.codecogs.com/png.latex?N_i
"N_i") is the number of secondary ions counted by the detector during a
single measurement cycle.

### Predicted standard deviation

The predicted standard deviation of a whole count block is directly
related to the population mean of
![N\_{(i)}](https://latex.codecogs.com/png.latex?N_%7B%28i%29%7D
"N_{(i)}")
(![\\mu\_{N}](https://latex.codecogs.com/png.latex?%5Cmu_%7BN%7D
"\\mu_{N}")) by the equation;

  
![ \\sigma =
\\sqrt(\\mu\_{N})](https://latex.codecogs.com/png.latex?%20%5Csigma%20%3D%20%5Csqrt%28%5Cmu_%7BN%7D%29
" \\sigma = \\sqrt(\\mu_{N})")  
{\#eq:base.pois}

In this formulation the population mean of N
(![\\mu\_{N}](https://latex.codecogs.com/png.latex?%5Cmu_%7BN%7D
"\\mu_{N}")) can be substituted by the mean number of events
(i.e. secondary ion counts) per time unit, or
![\\bar{N}](https://latex.codecogs.com/png.latex?%5Cbar%7BN%7D
"\\bar{N}"). The predicted standard deviation can therefore be deduced
from the mean number of counts for that particular ion per count block,
as follows

  
![\\hat{s}\_{N} =
\\sqrt{\\bar{N}}](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_%7BN%7D%20%3D%20%5Csqrt%7B%5Cbar%7BN%7D%7D
"\\hat{s}_{N} = \\sqrt{\\bar{N}}")  
{\#eq:std.pois}

where:

  
![\\bar{N} =
\\frac{1}{n}\\sum\_{i=1}^{n}N\_i](https://latex.codecogs.com/png.latex?%5Cbar%7BN%7D%20%3D%20%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7Bi%3D1%7D%5E%7Bn%7DN_i
"\\bar{N} = \\frac{1}{n}\\sum_{i=1}^{n}N_i")  
{\#eq:mean.N}

In this formulation, the hat on
![\\hat{s}\_x](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_x
"\\hat{s}_x") denotes that the statistics is predictive, instead of
![s\_x](https://latex.codecogs.com/png.latex?s_x "s_x") which is an
observed value. The commonality of the two measures is, however that
they are a estimate of the true population
![\\sigma](https://latex.codecogs.com/png.latex?%5Csigma "\\sigma").

### Predicted standard error of the mean

In a similar fashion the standard error of the mean for Poisson
statistics depends again on the number of measurements
![n](https://latex.codecogs.com/png.latex?n "n"), and can be formulated
as follows

  
![\\hat{s}\_{\\bar{N}} = \\sqrt{\\left( \\frac{
\\bar{N}}{n}\\right)}](https://latex.codecogs.com/png.latex?%5Chat%7Bs%7D_%7B%5Cbar%7BN%7D%7D%20%3D%20%5Csqrt%7B%5Cleft%28%20%5Cfrac%7B%20%5Cbar%7BN%7D%7D%7Bn%7D%5Cright%29%7D
"\\hat{s}_{\\bar{N}} = \\sqrt{\\left( \\frac{ \\bar{N}}{n}\\right)}")  
{\#eq:se.pois}

### Error propagation for predictive statistics of isotope ratios

For SIMS isotope analysis we need to have at least two different count
blocks, so that we can get a count ratio, as defined by \[@eq:mean.R\],
and where ![X\_i](https://latex.codecogs.com/png.latex?X_i "X_i") is a
time normalised count, or count rate. Satisfying this assumption
provides us with count-rate ratio
![R](https://latex.codecogs.com/png.latex?R "R") for measurement
![i](https://latex.codecogs.com/png.latex?i "i") of the isotopes
![a](https://latex.codecogs.com/png.latex?a "a") and
![b](https://latex.codecogs.com/png.latex?b "b"), where we take a mean
![\\bar{x}\_{R}](https://latex.codecogs.com/png.latex?%5Cbar%7Bx%7D_%7BR%7D
"\\bar{x}_{R}") from the completed count block as our estimate of the
true isotope value
![\\mu\_R](https://latex.codecogs.com/png.latex?%5Cmu_R "\\mu_R"). As
the predicted
![\\hat{S}\_x](https://latex.codecogs.com/png.latex?%5Chat%7BS%7D_x
"\\hat{S}_x") can be calculated for single ions, this should also mean
that the uncertainty in the isotope measurement can be predicted
(![\\hat{S}\_R](https://latex.codecogs.com/png.latex?%5Chat%7BS%7D_R
"\\hat{S}_R")). And, again this requires proper error propagation to
incorporate the cumulative errors on the counts of both isotopes;
![N^{a}](https://latex.codecogs.com/png.latex?N%5E%7Ba%7D "N^{a}") and
![N^{b}](https://latex.codecogs.com/png.latex?N%5E%7Bb%7D "N^{b}"), over
one count block \[@Fitzsimons2000a\]. Since the count-rate ratio
![R](https://latex.codecogs.com/png.latex?R "R") is a linear function of
the count ratio, it is possible to use the standard deviation of the
count ratio
![\\hat{S}\_{N^{b}/N^{a}}](https://latex.codecogs.com/png.latex?%5Chat%7BS%7D_%7BN%5E%7Bb%7D%2FN%5E%7Ba%7D%7D
"\\hat{S}_{N^{b}/N^{a}}") instead of
![\\hat{S}\_{R}](https://latex.codecogs.com/png.latex?%5Chat%7BS%7D_%7BR%7D
"\\hat{S}_{R}"), following that:

  
![ \\hat{S}\_{R} \\approx \\left(\\frac{t^{a}}{t^{b}} \\right)
\\hat{S}\_{N^{b}/N^{a}}](https://latex.codecogs.com/png.latex?%20%5Chat%7BS%7D_%7BR%7D%20%5Capprox%20%5Cleft%28%5Cfrac%7Bt%5E%7Ba%7D%7D%7Bt%5E%7Bb%7D%7D%20%5Cright%29%20%5Chat%7BS%7D_%7BN%5E%7Bb%7D%2FN%5E%7Ba%7D%7D
" \\hat{S}_{R} \\approx \\left(\\frac{t^{a}}{t^{b}} \\right) \\hat{S}_{N^{b}/N^{a}}")  
{\#eq:N\_R}

This provides the possibility to express
![\\hat{S}\_{N^{b}/N^{a}}](https://latex.codecogs.com/png.latex?%5Chat%7BS%7D_%7BN%5E%7Bb%7D%2FN%5E%7Ba%7D%7D
"\\hat{S}_{N^{b}/N^{a}}") in terms of the standard deviations of the
individual counts, and by using \[@eq:er.prop\], yielding;

  
![ \\hat{S}\_{N^{b}/N^{a}} \\approx \\sqrt{ \\left(
\\frac{\\hat{S}\_{N^{b}}}{N^{a}} \\right) + \\left(
\\frac{\\hat{S}\_{N^{b}}}{N^{a}} \\right) - \\frac{2r\_{N^{b}N^{a}}
s\_{N^{b}} s\_{N^{a}}}{N^{b}N^{a}} }\\times \\frac{N^{b}}{N^{a}}
](https://latex.codecogs.com/png.latex?%20%5Chat%7BS%7D_%7BN%5E%7Bb%7D%2FN%5E%7Ba%7D%7D%20%5Capprox%20%5Csqrt%7B%20%5Cleft%28%20%5Cfrac%7B%5Chat%7BS%7D_%7BN%5E%7Bb%7D%7D%7D%7BN%5E%7Ba%7D%7D%20%5Cright%29%20%2B%20%5Cleft%28%20%5Cfrac%7B%5Chat%7BS%7D_%7BN%5E%7Bb%7D%7D%7D%7BN%5E%7Ba%7D%7D%20%5Cright%29%20%20-%20%5Cfrac%7B2r_%7BN%5E%7Bb%7DN%5E%7Ba%7D%7D%20s_%7BN%5E%7Bb%7D%7D%20s_%7BN%5E%7Ba%7D%7D%7D%7BN%5E%7Bb%7DN%5E%7Ba%7D%7D%20%7D%5Ctimes%20%5Cfrac%7BN%5E%7Bb%7D%7D%7BN%5E%7Ba%7D%7D%20
" \\hat{S}_{N^{b}/N^{a}} \\approx \\sqrt{ \\left( \\frac{\\hat{S}_{N^{b}}}{N^{a}} \\right) + \\left( \\frac{\\hat{S}_{N^{b}}}{N^{a}} \\right)  - \\frac{2r_{N^{b}N^{a}} s_{N^{b}} s_{N^{a}}}{N^{b}N^{a}} }\\times \\frac{N^{b}}{N^{a}} ")  
{\#eq:std.pois.R1}

As the both count statistics are independent, the
![r](https://latex.codecogs.com/png.latex?r "r") becomes zero. The
predicted standard deviations for
![N^{b}](https://latex.codecogs.com/png.latex?N%5E%7Bb%7D "N^{b}") and
![N^{a}](https://latex.codecogs.com/png.latex?N%5E%7Ba%7D "N^{a}") can
be approximated by the population mean, according to \[@eq:std.pois\],
thereby transforming \[@eq:std.pois.R1\] to

  
![ \\hat{S}\_{N^{b}/N^{a}} \\approx \\sqrt{\\frac{1}{ \\bar{N}^{b}} +
\\frac{1}{ \\bar{N}^{a}}} \\times \\bar{N}^{b}/\\bar{N}^{a}
](https://latex.codecogs.com/png.latex?%20%5Chat%7BS%7D_%7BN%5E%7Bb%7D%2FN%5E%7Ba%7D%7D%20%5Capprox%20%5Csqrt%7B%5Cfrac%7B1%7D%7B%20%5Cbar%7BN%7D%5E%7Bb%7D%7D%20%2B%20%5Cfrac%7B1%7D%7B%20%5Cbar%7BN%7D%5E%7Ba%7D%7D%7D%20%5Ctimes%20%20%5Cbar%7BN%7D%5E%7Bb%7D%2F%5Cbar%7BN%7D%5E%7Ba%7D%20
" \\hat{S}_{N^{b}/N^{a}} \\approx \\sqrt{\\frac{1}{ \\bar{N}^{b}} + \\frac{1}{ \\bar{N}^{a}}} \\times  \\bar{N}^{b}/\\bar{N}^{a} ")  
{\#eq:std.pois.R2}

in which we can substitute \[@eq:N\_R\] to obtain

  
![ \\hat{S}\_{x} \\approx \\sqrt{\\frac{1}{ \\bar{N}^{b}} + \\frac{1}{
\\bar{N}^{a}}} \\times \\frac{\\bar{N}^{b}}{\\bar{N}^{a}} \\left(
\\frac{t^{a}}{t^{b}} \\right)
](https://latex.codecogs.com/png.latex?%20%5Chat%7BS%7D_%7Bx%7D%20%5Capprox%20%5Csqrt%7B%5Cfrac%7B1%7D%7B%20%5Cbar%7BN%7D%5E%7Bb%7D%7D%20%2B%20%5Cfrac%7B1%7D%7B%20%5Cbar%7BN%7D%5E%7Ba%7D%7D%7D%20%5Ctimes%20%20%5Cfrac%7B%5Cbar%7BN%7D%5E%7Bb%7D%7D%7B%5Cbar%7BN%7D%5E%7Ba%7D%7D%20%5Cleft%28%20%5Cfrac%7Bt%5E%7Ba%7D%7D%7Bt%5E%7Bb%7D%7D%20%5Cright%29%20
" \\hat{S}_{x} \\approx \\sqrt{\\frac{1}{ \\bar{N}^{b}} + \\frac{1}{ \\bar{N}^{a}}} \\times  \\frac{\\bar{N}^{b}}{\\bar{N}^{a}} \\left( \\frac{t^{a}}{t^{b}} \\right) ")  
{\#eq:std.pois.R3}

, which is equivalent to

  
![ \\hat{S}\_{N^{b}/N^{a}} \\approx \\sqrt{\\frac{1}{ \\bar{N}^{b}} +
\\frac{1}{ \\bar{N}^{a}}} \\times
\\bar{x}\_{R}](https://latex.codecogs.com/png.latex?%20%5Chat%7BS%7D_%7BN%5E%7Bb%7D%2FN%5E%7Ba%7D%7D%20%5Capprox%20%5Csqrt%7B%5Cfrac%7B1%7D%7B%20%5Cbar%7BN%7D%5E%7Bb%7D%7D%20%2B%20%5Cfrac%7B1%7D%7B%20%5Cbar%7BN%7D%5E%7Ba%7D%7D%7D%20%5Ctimes%20%20%5Cbar%7Bx%7D_%7BR%7D
" \\hat{S}_{N^{b}/N^{a}} \\approx \\sqrt{\\frac{1}{ \\bar{N}^{b}} + \\frac{1}{ \\bar{N}^{a}}} \\times  \\bar{x}_{R}")  
{\#eq:std.pois.R4}

In \[@eq:std.pois.R4\], we can substitute \[@eq:mean.N\] for
![\\bar{N}^{b}](https://latex.codecogs.com/png.latex?%5Cbar%7BN%7D%5E%7Bb%7D
"\\bar{N}^{b}") and
![\\bar{N}^{b}](https://latex.codecogs.com/png.latex?%5Cbar%7BN%7D%5E%7Bb%7D
"\\bar{N}^{b}"), respectively.

  
![ \\hat{s}\_{x\_{R}} = 
\\sqrt{ 
\\left( 
\\frac{1}{\\sum\_{i = 1}^{n}{N\_i^a}} \\right) + 
\\left( 
\\frac{1}{\\sum\_{i = 1}^{n}{N\_i^b}} \\right)} \\times \\bar{x}\_R
\\sqrt{n}
](https://latex.codecogs.com/png.latex?%20%5Chat%7Bs%7D_%7Bx_%7BR%7D%7D%20%3D%20%0A%20%20%20%20%20%20%20%20%5Csqrt%7B%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cleft%28%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cfrac%7B1%7D%7B%5Csum_%7Bi%20%3D%201%7D%5E%7Bn%7D%7BN_i%5Ea%7D%7D%20%5Cright%29%20%20%2B%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cleft%28%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cfrac%7B1%7D%7B%5Csum_%7Bi%20%3D%201%7D%5E%7Bn%7D%7BN_i%5Eb%7D%7D%20%5Cright%29%7D%20%5Ctimes%20%5Cbar%7Bx%7D_R%20%5Csqrt%7Bn%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20
" \\hat{s}_{x_{R}} = 
        \\sqrt{ 
            \\left( 
                \\frac{1}{\\sum_{i = 1}^{n}{N_i^a}} \\right)  + 
             \\left( 
                 \\frac{1}{\\sum_{i = 1}^{n}{N_i^b}} \\right)} \\times \\bar{x}_R \\sqrt{n}
                ")  
{\#eq:std.pois.R5}

The predicted standard error of the mean of a repeated set of
measurements in one count block is then:

  
![ \\hat{s}\_{\\bar{x}\_{R}} = 
\\sqrt{ 
\\left( 
\\frac{1}{\\sum\_{i = 1}^{n}{N\_i^a}} \\right) + 
\\left( 
\\frac{1}{\\sum\_{i = 1}^{n}{N\_i^b}} \\right)} \\times
\\bar{x}\_R](https://latex.codecogs.com/png.latex?%20%5Chat%7Bs%7D_%7B%5Cbar%7Bx%7D_%7BR%7D%7D%20%3D%20%0A%20%20%20%20%20%20%20%20%5Csqrt%7B%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cleft%28%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cfrac%7B1%7D%7B%5Csum_%7Bi%20%3D%201%7D%5E%7Bn%7D%7BN_i%5Ea%7D%7D%20%5Cright%29%20%20%2B%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cleft%28%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Cfrac%7B1%7D%7B%5Csum_%7Bi%20%3D%201%7D%5E%7Bn%7D%7BN_i%5Eb%7D%7D%20%5Cright%29%7D%20%20%5Ctimes%20%20%5Cbar%7Bx%7D_R
" \\hat{s}_{\\bar{x}_{R}} = 
        \\sqrt{ 
            \\left( 
                \\frac{1}{\\sum_{i = 1}^{n}{N_i^a}} \\right)  + 
             \\left( 
                 \\frac{1}{\\sum_{i = 1}^{n}{N_i^b}} \\right)}  \\times  \\bar{x}_R")  
{\#eq:std.pois.R6}

The latter to measures can be expressed as relative uncertainties in ‰,
following the same transformation as for the descriptive statistics.

## Verdict predictive and descriptive statistics

The ![\\chi^2](https://latex.codecogs.com/png.latex?%5Cchi%5E2
"\\chi^2") is used how well the machine approximated the theoretical
precision for isotope analysis:

  
![\\chi^2 = \\left( \\frac{s\_{\\bar{x}\_R}} {\\hat{s}\_{\\bar{x}\_R}}
\\right)^2
](https://latex.codecogs.com/png.latex?%5Cchi%5E2%20%3D%20%5Cleft%28%20%5Cfrac%7Bs_%7B%5Cbar%7Bx%7D_R%7D%7D%20%7B%5Chat%7Bs%7D_%7B%5Cbar%7Bx%7D_R%7D%7D%20%5Cright%29%5E2%20%20%20
"\\chi^2 = \\left( \\frac{s_{\\bar{x}_R}} {\\hat{s}_{\\bar{x}_R}} \\right)^2   ")
