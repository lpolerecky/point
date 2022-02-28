# consistency of the QSA test

    Code
      QSA_test(real_IC, "13C", "12C", file.nm)
    Output
      # A tibble: 11,700 x 16
         file.nm        t.nm sample.nm.12C sample.nm.13C bl.nm.12C bl.nm.13C Xt.pr.12C
         <chr>         <dbl> <chr>         <chr>             <int>     <int>     <dbl>
       1 2018-01-19-G~  0.54 Belemnite,In~ Belemnite,In~         1         1    34460.
       2 2018-01-19-G~  1.08 Belemnite,In~ Belemnite,In~         1         1    34202.
       3 2018-01-19-G~  1.62 Belemnite,In~ Belemnite,In~         1         1    34632.
       4 2018-01-19-G~  2.16 Belemnite,In~ Belemnite,In~         1         1    34191.
       5 2018-01-19-G~  2.7  Belemnite,In~ Belemnite,In~         1         1    34855.
       6 2018-01-19-G~  3.24 Belemnite,In~ Belemnite,In~         1         1    34672.
       7 2018-01-19-G~  3.78 Belemnite,In~ Belemnite,In~         1         1    34766.
       8 2018-01-19-G~  4.32 Belemnite,In~ Belemnite,In~         1         1    34609.
       9 2018-01-19-G~  4.86 Belemnite,In~ Belemnite,In~         1         1    34414.
      10 2018-01-19-G~  5.4  Belemnite,In~ Belemnite,In~         1         1    34474.
      # ... with 11,690 more rows, and 9 more variables: Xt.pr.13C <dbl>,
      #   N.pr.12C <dbl>, N.pr.13C <dbl>, R_Xt.pr <dbl>, alpha_Xt.pr.12C <dbl>,
      #   beta_Xt.pr.12C <dbl>, t_Xt.pr.12C <dbl>, p_Xt.pr.12C <dbl>,
      #   delta_Xt.pr.12C <dbl>

---

    Code
      QSA_test(real_IC, "13C", "12C", sample.nm, file.nm, .nest = file.nm)
    Output
      # A tibble: 11,700 x 20
         file.nm      sample.nm  t.nm bl.nm.12C bl.nm.13C Xt.pr.12C Xt.pr.13C N.pr.12C
         <chr>        <chr>     <dbl>     <int>     <int>     <dbl>     <dbl>    <dbl>
       1 2018-01-19-~ Belemnit~  0.54         1         1    34460.      358.    12040
       2 2018-01-19-~ Belemnit~  1.08         1         1    34202.      341.    11950
       3 2018-01-19-~ Belemnit~  1.62         1         1    34632.      395.    12100
       4 2018-01-19-~ Belemnit~  2.16         1         1    34191.      366.    11946
       5 2018-01-19-~ Belemnit~  2.7          1         1    34855.      378.    12178
       6 2018-01-19-~ Belemnit~  3.24         1         1    34672.      369.    12114
       7 2018-01-19-~ Belemnit~  3.78         1         1    34766.      369.    12147
       8 2018-01-19-~ Belemnit~  4.32         1         1    34609.      366.    12092
       9 2018-01-19-~ Belemnit~  4.86         1         1    34414.      406.    12024
      10 2018-01-19-~ Belemnit~  5.4          1         1    34474.      409.    12045
      # ... with 11,690 more rows, and 12 more variables: N.pr.13C <dbl>,
      #   R_Xt.pr <dbl>, alpha_Xt.pr.12C <dbl>, beta_Xt.pr.12C <dbl>,
      #   t_Xt.pr.12C <dbl>, p_Xt.pr.12C <dbl>, delta_Xt.pr.12C <dbl>,
      #   alpha_file.nm <dbl>, beta_file.nm <dbl>, t_file.nm <dbl>, p_file.nm <dbl>,
      #   delta_file.nm <dbl>

