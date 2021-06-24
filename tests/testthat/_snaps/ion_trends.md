# consistency of diagnostics for single ion trends

    Code
      predict_ionize(tb_0, file.nm, .plot = FALSE)
    Output
      # A tibble: 23,400 x 11
         file.nm  species.nm  t.nm sample.nm bl.nm  Xt.pr  N.pr Xt.nlt M_Xt.nlt  Xt.l0
         <chr>    <chr>      <dbl> <chr>     <int>  <dbl> <dbl>  <dbl>    <dbl>  <dbl>
       1 2018-01~ 12C         0.54 Belemnit~     1 34499. 12040 35215.   30345. 29628.
       2 2018-01~ 12C         1.08 Belemnit~     1 34241. 11950 35216.   30345. 29370.
       3 2018-01~ 12C         1.62 Belemnit~     1 34670. 12100 35216.   30345. 29799.
       4 2018-01~ 12C         2.16 Belemnit~     1 34229. 11946 35217.   30345. 29358.
       5 2018-01~ 12C         2.7  Belemnit~     1 34894. 12178 35217.   30345. 30022.
       6 2018-01~ 12C         3.24 Belemnit~     1 34711. 12114 35218.   30345. 29838.
       7 2018-01~ 12C         3.78 Belemnit~     1 34805. 12147 35219.   30345. 29932.
       8 2018-01~ 12C         4.32 Belemnit~     1 34648. 12092 35219.   30345. 29774.
       9 2018-01~ 12C         4.86 Belemnit~     1 34453. 12024 35220.   30345. 29578.
      10 2018-01~ 12C         5.4  Belemnit~     1 34513. 12045 35220.   30345. 29638.
      # ... with 23,390 more rows, and 1 more variable: N.l0 <dbl>

