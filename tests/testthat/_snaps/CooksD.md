# consistency of residual diagnostics on internal dataset

    Code
      CooksD(tb_R, "13C", "12C", file.nm, .output = "flag")
    Output
      # A tibble: 11,700 x 5
         file.nm                 t.nm hat_Xt.pr.13C     CooksD flag     
         <chr>                  <dbl>         <dbl>      <dbl> <fct>    
       1 2018-01-19-GLENDON_1_1  0.54          379. 0.000112   confluent
       2 2018-01-19-GLENDON_1_1  1.08          376. 0.000322   confluent
       3 2018-01-19-GLENDON_1_1  1.62          381. 0.0000564  confluent
       4 2018-01-19-GLENDON_1_1  2.16          376. 0.0000218  confluent
       5 2018-01-19-GLENDON_1_1  2.7           383. 0.00000647 confluent
       6 2018-01-19-GLENDON_1_1  3.24          381. 0.0000349  confluent
       7 2018-01-19-GLENDON_1_1  3.78          382. 0.0000415  confluent
       8 2018-01-19-GLENDON_1_1  4.32          380. 0.0000493  confluent
       9 2018-01-19-GLENDON_1_1  4.86          378. 0.000213   confluent
      10 2018-01-19-GLENDON_1_1  5.4           379. 0.000247   confluent
      # ... with 11,690 more rows

---

    Code
      Rm(tb_R, "13C", "12C", file.nm, .output = "flag")
    Output
      # A tibble: 11,700 x 5
         file.nm                 t.nm hat_Xt.pr.13C  studE flag     
         <chr>                  <dbl>         <dbl>  <dbl> <fct>    
       1 2018-01-19-GLENDON_1_1  0.54          379. -0.623 confluent
       2 2018-01-19-GLENDON_1_1  1.08          376. -1.06  confluent
       3 2018-01-19-GLENDON_1_1  1.62          381.  0.441 confluent
       4 2018-01-19-GLENDON_1_1  2.16          376. -0.276 confluent
       5 2018-01-19-GLENDON_1_1  2.7           383. -0.149 confluent
       6 2018-01-19-GLENDON_1_1  3.24          381. -0.347 confluent
       7 2018-01-19-GLENDON_1_1  3.78          382. -0.378 confluent
       8 2018-01-19-GLENDON_1_1  4.32          380. -0.412 confluent
       9 2018-01-19-GLENDON_1_1  4.86          378.  0.860 confluent
      10 2018-01-19-GLENDON_1_1  5.4           379.  0.925 confluent
      # ... with 11,690 more rows

---

    Code
      CV(tb_R, "13C", "12C", file.nm, .output = "flag")
    Output
      # A tibble: 11,700 x 5
         file.nm                 t.nm hat_Xt.pr.13C  studE flag     
         <chr>                  <dbl>         <dbl>  <dbl> <fct>    
       1 2018-01-19-GLENDON_1_1  0.54          379. -0.623 confluent
       2 2018-01-19-GLENDON_1_1  1.08          376. -1.06  confluent
       3 2018-01-19-GLENDON_1_1  1.62          381.  0.441 confluent
       4 2018-01-19-GLENDON_1_1  2.16          376. -0.276 confluent
       5 2018-01-19-GLENDON_1_1  2.7           383. -0.149 confluent
       6 2018-01-19-GLENDON_1_1  3.24          381. -0.347 confluent
       7 2018-01-19-GLENDON_1_1  3.78          382. -0.378 confluent
       8 2018-01-19-GLENDON_1_1  4.32          380. -0.412 confluent
       9 2018-01-19-GLENDON_1_1  4.86          378.  0.860 confluent
      10 2018-01-19-GLENDON_1_1  5.4           379.  0.925 confluent
      # ... with 11,690 more rows

---

    Code
      QQ(tb_R, "13C", "12C", file.nm, .output = "flag")
    Output
      # A tibble: 11,700 x 8
         file.nm                 t.nm    RQ    TQ    QE hat_RQ hat_e_RQ flag     
         <chr>                  <dbl> <dbl> <dbl> <dbl>  <dbl>    <dbl> <fct>    
       1 2018-01-19-GLENDON_1_1  0.54 -3.08 -3.65 0.570  -3.64    0.353 confluent
       2 2018-01-19-GLENDON_1_1  1.08 -2.99 -3.36 0.366  -3.35    0.220 confluent
       3 2018-01-19-GLENDON_1_1  1.62 -2.93 -3.21 0.282  -3.21    0.177 confluent
       4 2018-01-19-GLENDON_1_1  2.16 -2.87 -3.12 0.242  -3.11    0.154 confluent
       5 2018-01-19-GLENDON_1_1  2.7  -2.81 -3.04 0.229  -3.04    0.139 confluent
       6 2018-01-19-GLENDON_1_1  3.24 -2.74 -2.98 0.240  -2.97    0.128 confluent
       7 2018-01-19-GLENDON_1_1  3.78 -2.66 -2.93 0.265  -2.92    0.119 divergent
       8 2018-01-19-GLENDON_1_1  4.32 -2.63 -2.88 0.251  -2.88    0.113 divergent
       9 2018-01-19-GLENDON_1_1  4.86 -2.59 -2.85 0.260  -2.84    0.107 divergent
      10 2018-01-19-GLENDON_1_1  5.4  -2.56 -2.81 0.254  -2.80    0.102 divergent
      # ... with 11,690 more rows

---

    Code
      norm_E(tb_R, "13C", "12C", file.nm, .output = "flag")
    Output
      # A tibble: 11,700 x 6
         file.nm                 t.nm  studE   hat_Xi     CooksD flag     
         <chr>                  <dbl>  <dbl>    <dbl>      <dbl> <fct>    
       1 2018-01-19-GLENDON_1_1  0.54 -0.623 0.000289 0.000112   confluent
       2 2018-01-19-GLENDON_1_1  1.08 -1.06  0.000286 0.000322   confluent
       3 2018-01-19-GLENDON_1_1  1.62  0.441 0.000290 0.0000564  confluent
       4 2018-01-19-GLENDON_1_1  2.16 -0.276 0.000286 0.0000218  confluent
       5 2018-01-19-GLENDON_1_1  2.7  -0.149 0.000292 0.00000647 confluent
       6 2018-01-19-GLENDON_1_1  3.24 -0.347 0.000290 0.0000349  confluent
       7 2018-01-19-GLENDON_1_1  3.78 -0.378 0.000291 0.0000415  confluent
       8 2018-01-19-GLENDON_1_1  4.32 -0.412 0.000290 0.0000493  confluent
       9 2018-01-19-GLENDON_1_1  4.86  0.860 0.000288 0.000213   confluent
      10 2018-01-19-GLENDON_1_1  5.4   0.925 0.000289 0.000247   confluent
      # ... with 11,690 more rows

---

    Code
      IR(tb_R, "13C", "12C", file.nm, .output = "flag")
    Output
      # A tibble: 105 x 5
         file.nm                  lag       acf  e_acf flag     
         <chr>                  <dbl>     <dbl>  <dbl> <fct>    
       1 2018-01-19-GLENDON_1_1     1 -0.00222  0.0314 confluent
       2 2018-01-19-GLENDON_1_1     2  0.00421  0.0314 confluent
       3 2018-01-19-GLENDON_1_1     3 -0.00344  0.0314 confluent
       4 2018-01-19-GLENDON_1_1     4  0.00351  0.0314 confluent
       5 2018-01-19-GLENDON_1_1     5 -0.00730  0.0314 confluent
       6 2018-01-19-GLENDON_1_1     6  0.000603 0.0314 confluent
       7 2018-01-19-GLENDON_1_1     7 -0.0214   0.0314 confluent
       8 2018-01-19-GLENDON_1_1     8  0.00584  0.0314 confluent
       9 2018-01-19-GLENDON_1_1     9 -0.00431  0.0314 confluent
      10 2018-01-19-GLENDON_1_1    10 -0.0133   0.0314 confluent
      # ... with 95 more rows

