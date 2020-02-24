
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



kableExtra::kable(tb.R[[1]] %>% rename(!!!tb.R[[2]]), escape = FALSE)
```

<table>

<thead>

<tr>

<th style="text-align:left;">

file.nm

</th>

<th style="text-align:left;">

R.nm

</th>

<th style="text-align:right;">

\[n\]

</th>

<th style="text-align:right;">

\[\bar{x}\]

</th>

<th style="text-align:right;">

\[r\]

</th>

<th style="text-align:right;">

\[s_x\]

</th>

<th style="text-align:right;">

\[\epsilon_x\]

</th>

<th style="text-align:right;">

\[s_\bar{x}\]

</th>

<th style="text-align:right;">

\[\epsilon_\bar{x}\]

</th>

<th style="text-align:right;">

\[\hat{s}_x\]

</th>

<th style="text-align:right;">

\[\hat{\epsilon}_x\]

</th>

<th style="text-align:right;">

\[\hat{s}_\bar{x}\]

</th>

<th style="text-align:right;">

\[\hat{\epsilon}_\bar{x}\]

</th>

<th style="text-align:right;">

\[\chi^2\]

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_1\_1.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

3900

</td>

<td style="text-align:right;">

0.0109721

</td>

<td style="text-align:right;">

0.6777798

</td>

<td style="text-align:right;">

0.0010201

</td>

<td style="text-align:right;">

92.96876

</td>

<td style="text-align:right;">

1.63e-05

</td>

<td style="text-align:right;">

1.488692

</td>

<td style="text-align:right;">

0.0010179

</td>

<td style="text-align:right;">

92.76863

</td>

<td style="text-align:right;">

1.63e-05

</td>

<td style="text-align:right;">

1.485487

</td>

<td style="text-align:right;">

1.0043193

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_1\_2.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

3900

</td>

<td style="text-align:right;">

0.0109778

</td>

<td style="text-align:right;">

0.7267907

</td>

<td style="text-align:right;">

0.0007774

</td>

<td style="text-align:right;">

70.81376

</td>

<td style="text-align:right;">

1.24e-05

</td>

<td style="text-align:right;">

1.133928

</td>

<td style="text-align:right;">

0.0007697

</td>

<td style="text-align:right;">

70.11032

</td>

<td style="text-align:right;">

1.23e-05

</td>

<td style="text-align:right;">

1.122664

</td>

<td style="text-align:right;">

1.0201674

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_1\_3.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

3900

</td>

<td style="text-align:right;">

0.0110027

</td>

<td style="text-align:right;">

0.5665313

</td>

<td style="text-align:right;">

0.0007318

</td>

<td style="text-align:right;">

66.51139

</td>

<td style="text-align:right;">

1.17e-05

</td>

<td style="text-align:right;">

1.065035

</td>

<td style="text-align:right;">

0.0007212

</td>

<td style="text-align:right;">

65.54477

</td>

<td style="text-align:right;">

1.15e-05

</td>

<td style="text-align:right;">

1.049556

</td>

<td style="text-align:right;">

1.0297123

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_2\_1.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

4000

</td>

<td style="text-align:right;">

0.0107790

</td>

<td style="text-align:right;">

0.9615434

</td>

<td style="text-align:right;">

0.0007599

</td>

<td style="text-align:right;">

70.49903

</td>

<td style="text-align:right;">

1.20e-05

</td>

<td style="text-align:right;">

1.114688

</td>

<td style="text-align:right;">

0.0007601

</td>

<td style="text-align:right;">

70.51348

</td>

<td style="text-align:right;">

1.20e-05

</td>

<td style="text-align:right;">

1.114916

</td>

<td style="text-align:right;">

0.9995903

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_2\_2.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

4000

</td>

<td style="text-align:right;">

0.0107809

</td>

<td style="text-align:right;">

0.9553736

</td>

<td style="text-align:right;">

0.0007530

</td>

<td style="text-align:right;">

69.84868

</td>

<td style="text-align:right;">

1.19e-05

</td>

<td style="text-align:right;">

1.104405

</td>

<td style="text-align:right;">

0.0007388

</td>

<td style="text-align:right;">

68.52923

</td>

<td style="text-align:right;">

1.17e-05

</td>

<td style="text-align:right;">

1.083542

</td>

<td style="text-align:right;">

1.0388783

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_2\_3.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

4000

</td>

<td style="text-align:right;">

0.0107487

</td>

<td style="text-align:right;">

0.9623700

</td>

<td style="text-align:right;">

0.0007365

</td>

<td style="text-align:right;">

68.52177

</td>

<td style="text-align:right;">

1.16e-05

</td>

<td style="text-align:right;">

1.083424

</td>

<td style="text-align:right;">

0.0007327

</td>

<td style="text-align:right;">

68.16257

</td>

<td style="text-align:right;">

1.16e-05

</td>

<td style="text-align:right;">

1.077745

</td>

<td style="text-align:right;">

1.0105675

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_2\_4.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

4000

</td>

<td style="text-align:right;">

0.0107693

</td>

<td style="text-align:right;">

0.9377646

</td>

<td style="text-align:right;">

0.0007194

</td>

<td style="text-align:right;">

66.79846

</td>

<td style="text-align:right;">

1.14e-05

</td>

<td style="text-align:right;">

1.056176

</td>

<td style="text-align:right;">

0.0007248

</td>

<td style="text-align:right;">

67.30553

</td>

<td style="text-align:right;">

1.15e-05

</td>

<td style="text-align:right;">

1.064194

</td>

<td style="text-align:right;">

0.9849892

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_2\_5.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

4000

</td>

<td style="text-align:right;">

0.0108907

</td>

<td style="text-align:right;">

0.9715567

</td>

<td style="text-align:right;">

0.0007725

</td>

<td style="text-align:right;">

70.93615

</td>

<td style="text-align:right;">

1.22e-05

</td>

<td style="text-align:right;">

1.121599

</td>

<td style="text-align:right;">

0.0007746

</td>

<td style="text-align:right;">

71.12562

</td>

<td style="text-align:right;">

1.22e-05

</td>

<td style="text-align:right;">

1.124595

</td>

<td style="text-align:right;">

0.9946793

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_2\_6.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

4000

</td>

<td style="text-align:right;">

0.0109161

</td>

<td style="text-align:right;">

0.9940903

</td>

<td style="text-align:right;">

0.0007044

</td>

<td style="text-align:right;">

64.52575

</td>

<td style="text-align:right;">

1.11e-05

</td>

<td style="text-align:right;">

1.020242

</td>

<td style="text-align:right;">

0.0007134

</td>

<td style="text-align:right;">

65.35719

</td>

<td style="text-align:right;">

1.13e-05

</td>

<td style="text-align:right;">

1.033388

</td>

<td style="text-align:right;">

0.9747189

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_2\_7.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

4000

</td>

<td style="text-align:right;">

0.0110876

</td>

<td style="text-align:right;">

0.9713950

</td>

<td style="text-align:right;">

0.0008162

</td>

<td style="text-align:right;">

73.61425

</td>

<td style="text-align:right;">

1.29e-05

</td>

<td style="text-align:right;">

1.163944

</td>

<td style="text-align:right;">

0.0008101

</td>

<td style="text-align:right;">

73.06514

</td>

<td style="text-align:right;">

1.28e-05

</td>

<td style="text-align:right;">

1.155261

</td>

<td style="text-align:right;">

1.0150874

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_2\_8.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

4000

</td>

<td style="text-align:right;">

0.0110901

</td>

<td style="text-align:right;">

0.9580354

</td>

<td style="text-align:right;">

0.0007566

</td>

<td style="text-align:right;">

68.22399

</td>

<td style="text-align:right;">

1.20e-05

</td>

<td style="text-align:right;">

1.078716

</td>

<td style="text-align:right;">

0.0007437

</td>

<td style="text-align:right;">

67.06142

</td>

<td style="text-align:right;">

1.18e-05

</td>

<td style="text-align:right;">

1.060334

</td>

<td style="text-align:right;">

1.0349721

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_2\_9.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

4000

</td>

<td style="text-align:right;">

0.0110818

</td>

<td style="text-align:right;">

0.9642637

</td>

<td style="text-align:right;">

0.0007135

</td>

<td style="text-align:right;">

64.38828

</td>

<td style="text-align:right;">

1.13e-05

</td>

<td style="text-align:right;">

1.018068

</td>

<td style="text-align:right;">

0.0007174

</td>

<td style="text-align:right;">

64.73397

</td>

<td style="text-align:right;">

1.13e-05

</td>

<td style="text-align:right;">

1.023534

</td>

<td style="text-align:right;">

0.9893481

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_3\_1.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

4000

</td>

<td style="text-align:right;">

0.0109466

</td>

<td style="text-align:right;">

0.9587252

</td>

<td style="text-align:right;">

0.0009169

</td>

<td style="text-align:right;">

83.75694

</td>

<td style="text-align:right;">

1.45e-05

</td>

<td style="text-align:right;">

1.324313

</td>

<td style="text-align:right;">

0.0009190

</td>

<td style="text-align:right;">

83.95670

</td>

<td style="text-align:right;">

1.45e-05

</td>

<td style="text-align:right;">

1.327472

</td>

<td style="text-align:right;">

0.9952471

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_3\_2.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

4000

</td>

<td style="text-align:right;">

0.0107554

</td>

<td style="text-align:right;">

0.9334001

</td>

<td style="text-align:right;">

0.0007319

</td>

<td style="text-align:right;">

68.05230

</td>

<td style="text-align:right;">

1.16e-05

</td>

<td style="text-align:right;">

1.076001

</td>

<td style="text-align:right;">

0.0007229

</td>

<td style="text-align:right;">

67.21633

</td>

<td style="text-align:right;">

1.14e-05

</td>

<td style="text-align:right;">

1.062783

</td>

<td style="text-align:right;">

1.0250287

</td>

</tr>

<tr>

<td style="text-align:left;">

2018-01-19-GLENDON\_3\_3.is\_txt

</td>

<td style="text-align:left;">

13C/12C

</td>

<td style="text-align:right;">

4000

</td>

<td style="text-align:right;">

0.0107471

</td>

<td style="text-align:right;">

0.9881505

</td>

<td style="text-align:right;">

0.0008391

</td>

<td style="text-align:right;">

78.07477

</td>

<td style="text-align:right;">

1.33e-05

</td>

<td style="text-align:right;">

1.234471

</td>

<td style="text-align:right;">

0.0008246

</td>

<td style="text-align:right;">

76.72916

</td>

<td style="text-align:right;">

1.30e-05

</td>

<td style="text-align:right;">

1.213194

</td>

<td style="text-align:right;">

1.0353818

</td>

</tr>

</tbody>

</table>
