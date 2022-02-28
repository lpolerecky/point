# consistency of the zeroCt and cov_R

    Code
      zeroCt(real_IC, "13C", "12C", file.nm)
    Output
      # A tibble: 23,400 x 7
         file.nm                 t.nm species.nm sample.nm        bl.nm  Xt.pr  N.pr
         <chr>                  <dbl> <chr>      <chr>            <int>  <dbl> <dbl>
       1 2018-01-19-GLENDON_1_1  0.54 12C        Belemnite,Indium     1 34460. 12040
       2 2018-01-19-GLENDON_1_1  1.08 12C        Belemnite,Indium     1 34202. 11950
       3 2018-01-19-GLENDON_1_1  1.62 12C        Belemnite,Indium     1 34632. 12100
       4 2018-01-19-GLENDON_1_1  2.16 12C        Belemnite,Indium     1 34191. 11946
       5 2018-01-19-GLENDON_1_1  2.7  12C        Belemnite,Indium     1 34855. 12178
       6 2018-01-19-GLENDON_1_1  3.24 12C        Belemnite,Indium     1 34672. 12114
       7 2018-01-19-GLENDON_1_1  3.78 12C        Belemnite,Indium     1 34766. 12147
       8 2018-01-19-GLENDON_1_1  4.32 12C        Belemnite,Indium     1 34609. 12092
       9 2018-01-19-GLENDON_1_1  4.86 12C        Belemnite,Indium     1 34414. 12024
      10 2018-01-19-GLENDON_1_1  5.4  12C        Belemnite,Indium     1 34474. 12045
      # ... with 23,390 more rows

---

    Code
      cov_R(real_IC, c("13C", "12C"), file.nm)
    Output
      # A tibble: 11,700 x 10
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
      # ... with 11,690 more rows, and 3 more variables: Xt.pr.13C <dbl>,
      #   N.pr.12C <dbl>, N.pr.13C <dbl>

# convert to wide format ion count data works

    Code
      cov_R(tb_pr, c("13C", "12C"), file.nm)
    Output
      # A tibble: 11,700 x 6
         file.nm                 t.nm Xt.pr.12C Xt.pr.13C N.pr.12C N.pr.13C
         <chr>                  <dbl>     <dbl>     <dbl>    <dbl>    <dbl>
       1 2018-01-19-GLENDON_1_1  0.54    34460.      358.    12040      125
       2 2018-01-19-GLENDON_1_1  1.08    34202.      341.    11950      119
       3 2018-01-19-GLENDON_1_1  1.62    34632.      395.    12100      138
       4 2018-01-19-GLENDON_1_1  2.16    34191.      366.    11946      128
       5 2018-01-19-GLENDON_1_1  2.7     34855.      378.    12178      132
       6 2018-01-19-GLENDON_1_1  3.24    34672.      369.    12114      129
       7 2018-01-19-GLENDON_1_1  3.78    34766.      369.    12147      129
       8 2018-01-19-GLENDON_1_1  4.32    34609.      366.    12092      128
       9 2018-01-19-GLENDON_1_1  4.86    34414.      406.    12024      142
      10 2018-01-19-GLENDON_1_1  5.4     34474.      409.    12045      143
      # ... with 11,690 more rows

