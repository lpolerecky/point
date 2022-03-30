# Cameca diagnostics on internal dataset are consistent

    Code
      Cameca(tb_R, "13C", "12C", file.nm, .X = Xt.pr, .N = N.pr, .species = species.nm,
        .t = t.nm, .output = "flag")
    Output
      # A tibble: 11,700 x 25
         file.nm       t.nm sample.nm bl.nm.12C bl.nm.13C Xt.pr.12C Xt.pr.13C N.pr.12C
         <chr>        <dbl> <chr>         <int>     <int>     <dbl>     <dbl>    <dbl>
       1 2018-01-19-~  0.54 Belemnit~         1         1    34460.      358.    12040
       2 2018-01-19-~  1.08 Belemnit~         1         1    34202.      341.    11950
       3 2018-01-19-~  1.62 Belemnit~         1         1    34632.      395.    12100
       4 2018-01-19-~  2.16 Belemnit~         1         1    34191.      366.    11946
       5 2018-01-19-~  2.7  Belemnit~         1         1    34855.      378.    12178
       6 2018-01-19-~  3.24 Belemnit~         1         1    34672.      369.    12114
       7 2018-01-19-~  3.78 Belemnit~         1         1    34766.      369.    12147
       8 2018-01-19-~  4.32 Belemnit~         1         1    34609.      366.    12092
       9 2018-01-19-~  4.86 Belemnit~         1         1    34414.      406.    12024
      10 2018-01-19-~  5.4  Belemnit~         1         1    34474.      409.    12045
      # ... with 11,690 more rows, and 17 more variables: N.pr.13C <dbl>,
      #   n_R_t.nm <int>, M_R_Xt.pr <dbl>, S_R_Xt.pr <dbl>, RS_R_Xt.pr <dbl>,
      #   SeM_R_Xt.pr <dbl>, RSeM_R_Xt.pr <dbl>, hat_S_R_N.pr <dbl>,
      #   hat_RS_R_N.pr <dbl>, hat_SeM_R_N.pr <dbl>, hat_RSeM_R_N.pr <dbl>,
      #   chi2_R_N.pr <dbl>, R_Xt.pr <dbl>, ratio.nm <chr>, hat_s_R_Xt.pr <dbl>,
      #   hat_Xt.pr.13C <dbl>, flag <fct>

