# consistency of the evaluation of diagnostics on synthetic data

    Code
      eval_diag(tb_dia, "13C", "12C", type.nm, spot.nm, .nest = type.nm, .X = Xt.pr,
        .N = N.pr, .species = species.nm, .t = t.nm)
    Output
      # A tibble: 9 x 11
        execution type.nm    spot.nm ratio.nm M_R_Xt.pr F_R_Xt.pr p_R_Xt.pr
            <dbl> <chr>        <int> <chr>        <dbl>     <dbl>     <dbl>
      1         1 symmetric        1 13C/12C     0.0110   105.     2.00e-64
      2         1 symmetric        2 13C/12C     0.0110   127.     2.81e-77
      3         1 symmetric        3 13C/12C     0.0110   122.     1.42e-74
      4         1 asymmetric       1 13C/12C     0.0111    63.2    1.28e-39
      5         1 asymmetric       2 13C/12C     0.0111    74.9    1.09e-46
      6         1 asymmetric       3 13C/12C     0.0111    71.9    6.70e-45
      7         1 ideal            1 13C/12C     0.0112     0.267  8.49e- 1
      8         1 ideal            2 13C/12C     0.0112     0.219  8.83e- 1
      9         1 ideal            3 13C/12C     0.0112     1.73   1.59e- 1
      # ... with 4 more variables: hat_M_M_R_Xt.pr <dbl>, hat_RS_M_R_Xt.pr <dbl>,
      #   dAIC_M_R_Xt.pr <dbl>, p_M_R_Xt.pr <dbl>

