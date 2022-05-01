# diagnostics wrapper on synthetic data is consistent

    Code
      diag_R(simu_IC, "13C", "12C", type.nm, spot.nm)
    Output
      # A tibble: 9 x 7
        execution type.nm    spot.nm ratio.nm M_R_Xt.pr F_R_Xt.pr p_R_Xt.pr
            <dbl> <chr>        <int> <chr>        <dbl>     <dbl>     <dbl>
      1         1 asymmetric       1 13C/12C     0.0111    63.2    1.28e-39
      2         1 asymmetric       2 13C/12C     0.0111    74.9    1.09e-46
      3         1 asymmetric       3 13C/12C     0.0111    71.9    6.70e-45
      4         1 ideal            1 13C/12C     0.0112     0.267  8.49e- 1
      5         1 ideal            2 13C/12C     0.0112     0.219  8.83e- 1
      6         1 ideal            3 13C/12C     0.0112     1.73   1.59e- 1
      7         1 symmetric        1 13C/12C     0.0110   105.     2.00e-64
      8         1 symmetric        2 13C/12C     0.0110   127.     2.81e-77
      9         1 symmetric        3 13C/12C     0.0110   122.     1.42e-74

---

    Code
      diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .method = "CV")
    Output
      # A tibble: 9 x 3
        type.nm    spot.nm hyp                  
        <chr>        <int> <chr>                
      1 symmetric        1 H0 (homoskedasticity)
      2 symmetric        2 H0 (homoskedasticity)
      3 symmetric        3 H0 (homoskedasticity)
      4 asymmetric       1 H0 (homoskedasticity)
      5 asymmetric       2 H0 (homoskedasticity)
      6 asymmetric       3 H0 (homoskedasticity)
      7 ideal            1 H0 (homoskedasticity)
      8 ideal            2 H0 (homoskedasticity)
      9 ideal            3 H0 (homoskedasticity)

---

    Code
      diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .output = "complete")
    Output
      # A tibble: 27,000 x 56
         execution type.nm    spot.nm  t.nm trend.nm.13C trend.nm.12C base.nm.13C
             <dbl> <chr>        <int> <int>        <dbl>        <dbl>       <dbl>
       1         1 asymmetric       1     1          120          120           0
       2         1 asymmetric       1     2          120          120           0
       3         1 asymmetric       1     3          120          120           0
       4         1 asymmetric       1     4          120          120           0
       5         1 asymmetric       1     5          120          120           0
       6         1 asymmetric       1     6          120          120           0
       7         1 asymmetric       1     7          120          120           0
       8         1 asymmetric       1     8          120          120           0
       9         1 asymmetric       1     9          120          120           0
      10         1 asymmetric       1    10          120          120           0
      # ... with 26,990 more rows, and 49 more variables: base.nm.12C <dbl>,
      #   force.nm.13C <dbl>, force.nm.12C <dbl>, bl.nm.13C <int>, bl.nm.12C <int>,
      #   n.rw.13C <dbl>, n.rw.12C <dbl>, N.pr.13C <dbl>, N.pr.12C <dbl>,
      #   Xt.pr.13C <dbl>, Xt.pr.12C <dbl>, n_t.nm.13C <int>, n_t.nm.12C <int>,
      #   tot_N.pr.13C <dbl>, tot_N.pr.12C <dbl>, M_Xt.pr.13C <dbl>,
      #   M_Xt.pr.12C <dbl>, S_Xt.pr.13C <dbl>, S_Xt.pr.12C <dbl>,
      #   RS_Xt.pr.13C <dbl>, RS_Xt.pr.12C <dbl>, SeM_Xt.pr.13C <dbl>, ...

---

    Code
      diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .output = "augmented")
    Output
      # A tibble: 104,980 x 12
         execution type.nm   trend.nm base.nm force.nm  t.nm bl.nm  n.rw spot.nm
             <dbl> <chr>        <dbl>   <dbl>    <dbl> <int> <int> <dbl>   <int>
       1         1 symmetric      120       0      -60     1     1  3000       1
       2         1 symmetric      120       0      -60     1     1  3000       1
       3         1 symmetric      120       0      -60     1     1  3000       2
       4         1 symmetric      120       0      -60     1     1  3000       2
       5         1 symmetric      120       0      -60     1     1  3000       3
       6         1 symmetric      120       0      -60     1     1  3000       3
       7         1 symmetric      120       0      -60     2     1  3000       1
       8         1 symmetric      120       0      -60     2     1  3000       1
       9         1 symmetric      120       0      -60     2     1  3000       2
      10         1 symmetric      120       0      -60     2     1  3000       2
      # ... with 104,970 more rows, and 3 more variables: species.nm <chr>,
      #   N.pr <dbl>, Xt.pr <dbl>

---

    Code
      diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .output = "diagnostic")
    Output
      # A tibble: 27,000 x 54
         execution type.nm   spot.nm  t.nm trend.nm.13C trend.nm.12C base.nm.13C
             <dbl> <chr>       <int> <int>        <dbl>        <dbl>       <dbl>
       1         1 symmetric       1     1          120          120           0
       2         1 symmetric       1     2          120          120           0
       3         1 symmetric       1     3          120          120           0
       4         1 symmetric       1     4          120          120           0
       5         1 symmetric       1     5          120          120           0
       6         1 symmetric       1     6          120          120           0
       7         1 symmetric       1     7          120          120           0
       8         1 symmetric       1     8          120          120           0
       9         1 symmetric       1     9          120          120           0
      10         1 symmetric       1    10          120          120           0
      # ... with 26,990 more rows, and 47 more variables: base.nm.12C <dbl>,
      #   force.nm.13C <dbl>, force.nm.12C <dbl>, bl.nm.13C <int>, bl.nm.12C <int>,
      #   n.rw.13C <dbl>, n.rw.12C <dbl>, N.pr.13C <dbl>, N.pr.12C <dbl>,
      #   Xt.pr.13C <dbl>, Xt.pr.12C <dbl>, n_t.nm.13C <int>, n_t.nm.12C <int>,
      #   tot_N.pr.13C <dbl>, tot_N.pr.12C <dbl>, M_Xt.pr.13C <dbl>,
      #   M_Xt.pr.12C <dbl>, S_Xt.pr.13C <dbl>, S_Xt.pr.12C <dbl>,
      #   RS_Xt.pr.13C <dbl>, RS_Xt.pr.12C <dbl>, SeM_Xt.pr.13C <dbl>, ...

---

    Code
      diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .output = "outlier")
    Output
      # A tibble: 27,000 x 22
         execution type.nm   spot.nm  t.nm trend.nm.13C trend.nm.12C base.nm.13C
             <dbl> <chr>       <int> <int>        <dbl>        <dbl>       <dbl>
       1         1 symmetric       1     1          120          120           0
       2         1 symmetric       1     2          120          120           0
       3         1 symmetric       1     3          120          120           0
       4         1 symmetric       1     4          120          120           0
       5         1 symmetric       1     5          120          120           0
       6         1 symmetric       1     6          120          120           0
       7         1 symmetric       1     7          120          120           0
       8         1 symmetric       1     8          120          120           0
       9         1 symmetric       1     9          120          120           0
      10         1 symmetric       1    10          120          120           0
      # ... with 26,990 more rows, and 15 more variables: base.nm.12C <dbl>,
      #   force.nm.13C <dbl>, force.nm.12C <dbl>, bl.nm.13C <int>, bl.nm.12C <int>,
      #   n.rw.13C <dbl>, n.rw.12C <dbl>, N.pr.13C <dbl>, N.pr.12C <dbl>,
      #   Xt.pr.13C <dbl>, Xt.pr.12C <dbl>, hat_S_N.pr.13C <dbl>,
      #   hat_Xt.pr.13C <dbl>, CooksD <dbl>, flag <fct>

---

    Code
      diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .label = "latex")
    Output
      # A tibble: 9 x 7
        execution type.nm    spot.nm ratio.nm `$\\bar{R}$` `$F_{R}$` `$p_{R}$`
            <dbl> <chr>        <int> <chr>           <dbl>     <dbl>     <dbl>
      1         1 asymmetric       1 13C/12C        0.0111    63.2    1.28e-39
      2         1 asymmetric       2 13C/12C        0.0111    74.9    1.09e-46
      3         1 asymmetric       3 13C/12C        0.0111    71.9    6.70e-45
      4         1 ideal            1 13C/12C        0.0112     0.267  8.49e- 1
      5         1 ideal            2 13C/12C        0.0112     0.219  8.83e- 1
      6         1 ideal            3 13C/12C        0.0112     1.73   1.59e- 1
      7         1 symmetric        1 13C/12C        0.0110   105.     2.00e-64
      8         1 symmetric        2 13C/12C        0.0110   127.     2.81e-77
      9         1 symmetric        3 13C/12C        0.0110   122.     1.42e-74

---

    Code
      diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .nest = type.nm, .label = "latex")
    Output
      # A tibble: 9 x 11
        execution type.nm    spot.nm ratio.nm `$\\bar{R}$` `$F_{R}$` `$p_{R}$`
            <dbl> <chr>        <int> <chr>           <dbl>     <dbl>     <dbl>
      1         1 asymmetric       1 13C/12C        0.0111    63.2    1.28e-39
      2         1 asymmetric       2 13C/12C        0.0111    74.9    1.09e-46
      3         1 asymmetric       3 13C/12C        0.0111    71.9    6.70e-45
      4         1 ideal            1 13C/12C        0.0112     0.267  8.49e- 1
      5         1 ideal            2 13C/12C        0.0112     0.219  8.83e- 1
      6         1 ideal            3 13C/12C        0.0112     1.73   1.59e- 1
      7         1 symmetric        1 13C/12C        0.0110   105.     2.00e-64
      8         1 symmetric        2 13C/12C        0.0110   127.     2.81e-77
      9         1 symmetric        3 13C/12C        0.0110   122.     1.42e-74
      # ... with 4 more variables: `$\\hat{\\bar{R}}$` <dbl>,
      #   `$\\hat{\\epsilon}_{\\bar{R}}$ (\\text{\\textperthousand})` <dbl>,
      #   `$\\Delta AIC_{\\bar{R}}$` <dbl>, `$p_{\\bar{R}}$` <dbl>

# QQ diagnostic on synthetic data is consistent

    Code
      diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .method = "QQ")
    Output
      # A tibble: 9 x 3
        type.nm    spot.nm hyp            
        <chr>        <int> <chr>          
      1 symmetric        1 H0 (normal)    
      2 symmetric        2 H0 (normal)    
      3 symmetric        3 H0 (normal)    
      4 asymmetric       1 Ha (non-normal)
      5 asymmetric       2 Ha (non-normal)
      6 asymmetric       3 Ha (non-normal)
      7 ideal            1 H0 (normal)    
      8 ideal            2 H0 (normal)    
      9 ideal            3 H0 (normal)    

# IR diagnostic on synthetic data is consistent

    Code
      diag_R(simu_IC, "13C", "12C", type.nm, spot.nm, .method = "IR")
    Output
      # A tibble: 9 x 3
        type.nm    spot.nm hyp                           
        <chr>        <int> <chr>                         
      1 symmetric        1 Ha (dependence of residuals)  
      2 symmetric        2 Ha (dependence of residuals)  
      3 symmetric        3 Ha (dependence of residuals)  
      4 asymmetric       1 Ha (dependence of residuals)  
      5 asymmetric       2 Ha (dependence of residuals)  
      6 asymmetric       3 Ha (dependence of residuals)  
      7 ideal            1 H0 (independence of residuals)
      8 ideal            2 H0 (independence of residuals)
      9 ideal            3 H0 (independence of residuals)

# diagnostics preserve metadata

    Code
      diag_R(real_IC, "13C", "12C", file.nm)
    Output
      # A tibble: 3 x 6
        execution file.nm                ratio.nm M_R_Xt.pr F_R_Xt.pr p_R_Xt.pr
            <dbl> <chr>                  <chr>        <dbl>     <dbl>     <dbl>
      1         1 2018-01-19-GLENDON_1_1 13C/12C     0.0110      2.13   0.0940 
      2         1 2018-01-19-GLENDON_1_2 13C/12C     0.0110      1.37   0.250  
      3         1 2018-01-19-GLENDON_1_3 13C/12C     0.0110      4.43   0.00410

---

    Code
      unfold(diag_R(real_IC, "13C", "12C", file.nm, .meta = TRUE))
    Output
      # A tibble: 81,900 x 43
         execution file.nm     ratio.nm M_R_Xt.pr F_R_Xt.pr p_R_Xt.pr  t.nm species.nm
             <dbl> <chr>       <chr>        <dbl>     <dbl>     <dbl> <dbl> <chr>     
       1         1 2018-01-19~ 13C/12C     0.0110      2.13    0.0940  0.54 12C       
       2         1 2018-01-19~ 13C/12C     0.0110      2.13    0.0940  1.08 12C       
       3         1 2018-01-19~ 13C/12C     0.0110      2.13    0.0940  1.62 12C       
       4         1 2018-01-19~ 13C/12C     0.0110      2.13    0.0940  2.16 12C       
       5         1 2018-01-19~ 13C/12C     0.0110      2.13    0.0940  2.7  12C       
       6         1 2018-01-19~ 13C/12C     0.0110      2.13    0.0940  3.24 12C       
       7         1 2018-01-19~ 13C/12C     0.0110      2.13    0.0940  3.78 12C       
       8         1 2018-01-19~ 13C/12C     0.0110      2.13    0.0940  4.32 12C       
       9         1 2018-01-19~ 13C/12C     0.0110      2.13    0.0940  4.86 12C       
      10         1 2018-01-19~ 13C/12C     0.0110      2.13    0.0940  5.4  12C       
      # ... with 81,890 more rows, and 35 more variables: sample.nm <chr>,
      #   bl.nm <int>, num.mt <dbl>, bfield.mt <dbl>, rad.mt <dbl>, mass.mt <chr>,
      #   tc.mt <dbl>, coord.mt <chr>, file_raw.mt <chr>, bl_num.mt <dbl>,
      #   meas_bl.mt <dbl>, rejection.mt <dbl>, slit.mt <chr>, lens.mt <chr>,
      #   presput.mt <chr>, rast_com.mt <dbl>, frame.mt <chr>, blank_rast.mt <chr>,
      #   raster.mt <chr>, tune.mt <chr>, reg_mode.mt <chr>, chk_frm.mt <dbl>,
      #   sec_ion_cent.mt <chr>, frame_sec_ion_cent.mt <chr>, width_hor.mt <dbl>, ...

