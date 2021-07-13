# Check change over time in read_IC

    Code
      read_IC(point_example("2018-01-19-GLENDON"))
    Output
      # A tibble: 81,900 x 4
         file.nm                 t.nm  N.rw species.nm
         <chr>                  <dbl> <dbl> <chr>     
       1 2018-01-19-GLENDON_1_1  0.54 12040 12C       
       2 2018-01-19-GLENDON_1_1  1.08 11950 12C       
       3 2018-01-19-GLENDON_1_1  1.62 12100 12C       
       4 2018-01-19-GLENDON_1_1  2.16 11946 12C       
       5 2018-01-19-GLENDON_1_1  2.7  12178 12C       
       6 2018-01-19-GLENDON_1_1  3.24 12114 12C       
       7 2018-01-19-GLENDON_1_1  3.78 12147 12C       
       8 2018-01-19-GLENDON_1_1  4.32 12092 12C       
       9 2018-01-19-GLENDON_1_1  4.86 12024 12C       
      10 2018-01-19-GLENDON_1_1  5.4  12045 12C       
      # ... with 81,890 more rows

---

    Code
      read_meta(point_example("2018-01-19-GLENDON"))
    Output
      [[1]]
      # A tibble: 3 x 28
        file.nm  sample.nm  coord.mt  file_raw.mt    bl_num.mt meas_bl.mt rejection.mt
        <chr>    <chr>      <chr>     <chr>              <dbl>      <dbl>        <dbl>
      1 2018-01~ Belemnite~ x=-12761~ "D:\\CamecaNa~        60         65            2
      2 2018-01~ Belemnite~ x=-12741~ "D:\\CamecaNa~        60         65            2
      3 2018-01~ Belemnite~ x=-12721~ "D:\\CamecaNa~        60         65            2
      # ... with 21 more variables: slit.mt <chr>, lens.mt <chr>, presput.mt <chr>,
      #   rast_com.mt <dbl>, frame.mt <chr>, blank_rast.mt <chr>, raster.mt <chr>,
      #   tune.mt <chr>, reg_mode.mt <chr>, chk_frm.mt <dbl>, sec_ion_cent.mt <chr>,
      #   frame_sec_ion_cent.mt <chr>, width_hor.mt <dbl>, width_ver.mt <dbl>,
      #   E0S_cent.mt <chr>, width_V.mt <dbl>, E0P_off.mt <dbl>,
      #   prim_cur_start.mt <chr>, prim_cur_after.mt <chr>, n.rw <dbl>,
      #   det_type.mt <chr>
      
      [[2]]
      # A tibble: 6 x 5
        file.nm                num.mt M_PHD.mt SD_PHD.mt EMHV.mt
        <chr>                   <dbl>    <dbl>     <dbl>   <dbl>
      1 2018-01-19-GLENDON_1_1      1     217.     23.6    2036.
      2 2018-01-19-GLENDON_1_1      7     217.     20.8       0 
      3 2018-01-19-GLENDON_1_2      1     218.      4.91   2036.
      4 2018-01-19-GLENDON_1_2      7     211.     42.8       0 
      5 2018-01-19-GLENDON_1_3      1     233.     38.8    2036.
      6 2018-01-19-GLENDON_1_3      7     218.      8.83      0 
      

