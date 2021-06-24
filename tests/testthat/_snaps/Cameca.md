# consistency of Cameca diagnostics on internal dataset

    Code
      Cameca(tb_R, "13C", "12C", file.nm, .output = "flag")
    Output
      # A tibble: 11,700 x 25
         file.nm     sample.nm   t.nm bl.nm.12C bl.nm.13C Xt.pr.12C Xt.pr.13C N.pr.12C
         <chr>       <chr>      <dbl>     <int>     <int>     <dbl>     <dbl>    <dbl>
       1 2018-01-19~ Belemnite~  0.54         1         1    34499.      358.    12040
       2 2018-01-19~ Belemnite~  1.08         1         1    34241.      341.    11950
       3 2018-01-19~ Belemnite~  1.62         1         1    34670.      395.    12100
       4 2018-01-19~ Belemnite~  2.16         1         1    34229.      367.    11946
       5 2018-01-19~ Belemnite~  2.7          1         1    34894.      378.    12178
       6 2018-01-19~ Belemnite~  3.24         1         1    34711.      370.    12114
       7 2018-01-19~ Belemnite~  3.78         1         1    34805.      370.    12147
       8 2018-01-19~ Belemnite~  4.32         1         1    34648.      367.    12092
       9 2018-01-19~ Belemnite~  4.86         1         1    34453.      407.    12024
      10 2018-01-19~ Belemnite~  5.4          1         1    34513.      410.    12045
      # ... with 11,690 more rows, and 17 more variables: N.pr.13C <dbl>,
      #   n_R_t.nm <int>, M_R_Xt.pr <dbl>, S_R_Xt.pr <dbl>, RS_R_Xt.pr <dbl>,
      #   SeM_R_Xt.pr <dbl>, RSeM_R_Xt.pr <dbl>, hat_S_R_N.pr <dbl>,
      #   hat_RS_R_N.pr <dbl>, hat_SeM_R_N.pr <dbl>, hat_RSeM_R_N.pr <dbl>,
      #   chi2_R_N.pr <dbl>, R_Xt.pr <dbl>, ratio.nm <chr>, hat_s_R_Xt.pr <dbl>,
      #   hat_Xt.pr.13C <dbl>, flag <fct>

