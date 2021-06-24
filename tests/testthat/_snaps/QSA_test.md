# consistency of the QSA test

    Code
      QSA_test(tb_pr, "13C", "12C", file.nm)
    Output
      # A tibble: 11,700 x 16
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
      # ... with 11,690 more rows, and 9 more variables: Xt.pr.13C <dbl>,
      #   N.pr.12C <dbl>, N.pr.13C <dbl>, R_Xt.pr <dbl>, alpha_Xt.pr.12C <dbl>,
      #   beta_Xt.pr.12C <dbl>, t_Xt.pr.12C <dbl>, p_Xt.pr.12C <dbl>,
      #   delta_Xt.pr.12C <dbl>

