# consistency of diagnostics for single ion trends

    Code
      predict_ionize(tb_0, file.nm, .plot = FALSE)
    Output
      # A tibble: 23,400 x 11
         file.nm  species.nm  t.nm sample.nm bl.nm  Xt.pr  N.pr Xt.nlt M_Xt.nlt  Xt.l0
         <chr>    <chr>      <dbl> <chr>     <int>  <dbl> <dbl>  <dbl>    <dbl>  <dbl>
       1 2018-01~ 12C         0.54 Belemnit~     1 34460. 12040 35176.   30311. 29595.
       2 2018-01~ 12C         1.08 Belemnit~     1 34202. 11950 35176.   30311. 29337.
       3 2018-01~ 12C         1.62 Belemnit~     1 34632. 12100 35177.   30311. 29766.
       4 2018-01~ 12C         2.16 Belemnit~     1 34191. 11946 35177.   30311. 29325.
       5 2018-01~ 12C         2.7  Belemnit~     1 34855. 12178 35178.   30311. 29988.
       6 2018-01~ 12C         3.24 Belemnit~     1 34672. 12114 35179.   30311. 29804.
       7 2018-01~ 12C         3.78 Belemnit~     1 34766. 12147 35179.   30311. 29898.
       8 2018-01~ 12C         4.32 Belemnit~     1 34609. 12092 35180.   30311. 29740.
       9 2018-01~ 12C         4.86 Belemnit~     1 34414. 12024 35180.   30311. 29545.
      10 2018-01~ 12C         5.4  Belemnit~     1 34474. 12045 35181.   30311. 29605.
      # ... with 23,390 more rows, and 1 more variable: N.l0 <dbl>

