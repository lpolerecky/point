# consistency of the zeroCt and cov_R

    Code
      zeroCt(tb_pr, "13C", "12C", file.nm)
    Output
      # A tibble: 23,400 x 7
         file.nm                 t.nm species.nm sample.nm         bl.nm  Xt.pr  N.pr
         <chr>                  <dbl> <chr>      <chr>             <int>  <dbl> <dbl>
       1 2018-01-19-GLENDON_1_1  0.54 12C        Belemnite, Indium     1 34499. 12040
       2 2018-01-19-GLENDON_1_1  1.08 12C        Belemnite, Indium     1 34241. 11950
       3 2018-01-19-GLENDON_1_1  1.62 12C        Belemnite, Indium     1 34670. 12100
       4 2018-01-19-GLENDON_1_1  2.16 12C        Belemnite, Indium     1 34229. 11946
       5 2018-01-19-GLENDON_1_1  2.7  12C        Belemnite, Indium     1 34894. 12178
       6 2018-01-19-GLENDON_1_1  3.24 12C        Belemnite, Indium     1 34711. 12114
       7 2018-01-19-GLENDON_1_1  3.78 12C        Belemnite, Indium     1 34805. 12147
       8 2018-01-19-GLENDON_1_1  4.32 12C        Belemnite, Indium     1 34648. 12092
       9 2018-01-19-GLENDON_1_1  4.86 12C        Belemnite, Indium     1 34453. 12024
      10 2018-01-19-GLENDON_1_1  5.4  12C        Belemnite, Indium     1 34513. 12045
      # ... with 23,390 more rows

---

    Code
      cov_R(tb_pr, c("13C", "12C"), file.nm)
    Output
      # A tibble: 11,700 x 10
         file.nm      t.nm sample.nm.12C  sample.nm.13C  bl.nm.12C bl.nm.13C Xt.pr.12C
         <chr>       <dbl> <chr>          <chr>              <int>     <int>     <dbl>
       1 2018-01-19~  0.54 Belemnite, In~ Belemnite, In~         1         1    34499.
       2 2018-01-19~  1.08 Belemnite, In~ Belemnite, In~         1         1    34241.
       3 2018-01-19~  1.62 Belemnite, In~ Belemnite, In~         1         1    34670.
       4 2018-01-19~  2.16 Belemnite, In~ Belemnite, In~         1         1    34229.
       5 2018-01-19~  2.7  Belemnite, In~ Belemnite, In~         1         1    34894.
       6 2018-01-19~  3.24 Belemnite, In~ Belemnite, In~         1         1    34711.
       7 2018-01-19~  3.78 Belemnite, In~ Belemnite, In~         1         1    34805.
       8 2018-01-19~  4.32 Belemnite, In~ Belemnite, In~         1         1    34648.
       9 2018-01-19~  4.86 Belemnite, In~ Belemnite, In~         1         1    34453.
      10 2018-01-19~  5.4  Belemnite, In~ Belemnite, In~         1         1    34513.
      # ... with 11,690 more rows, and 3 more variables: Xt.pr.13C <dbl>,
      #   N.pr.12C <dbl>, N.pr.13C <dbl>

