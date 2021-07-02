# Check change over time in read_IC

    Code
      read_IC(path_point)
    Output
      # A tibble: 81,900 x 5
         file.nm                 t.nm  N.rw species.nm  n.rw
         <chr>                  <dbl> <dbl> <chr>      <dbl>
       1 2018-01-19-GLENDON_1_1  0.54 12040 12C         3900
       2 2018-01-19-GLENDON_1_1  1.08 11950 12C         3900
       3 2018-01-19-GLENDON_1_1  1.62 12100 12C         3900
       4 2018-01-19-GLENDON_1_1  2.16 11946 12C         3900
       5 2018-01-19-GLENDON_1_1  2.7  12178 12C         3900
       6 2018-01-19-GLENDON_1_1  3.24 12114 12C         3900
       7 2018-01-19-GLENDON_1_1  3.78 12147 12C         3900
       8 2018-01-19-GLENDON_1_1  4.32 12092 12C         3900
       9 2018-01-19-GLENDON_1_1  4.86 12024 12C         3900
      10 2018-01-19-GLENDON_1_1  5.4  12045 12C         3900
      # ... with 81,890 more rows

---

    Code
      read_meta(path_point)
    Output
      # A tibble: 6 x 31
        file.nm  num.mt mean_PHD.mt SD_PHD.mt EMHV.mt coord.mt  file_raw.mt  bl_num.mt
        <chr>     <dbl>       <dbl>     <dbl>   <dbl> <chr>     <chr>            <dbl>
      1 2018-01~      1        217.     23.6    2036. x=-12761~ "D:\\Cameca~        60
      2 2018-01~      7        217.     20.8       0  x=-12761~ "D:\\Cameca~        60
      3 2018-01~      1        218.      4.91   2036. x=-12741~ "D:\\Cameca~        60
      4 2018-01~      7        211.     42.8       0  x=-12741~ "D:\\Cameca~        60
      5 2018-01~      1        233.     38.8    2036. x=-12721~ "D:\\Cameca~        60
      6 2018-01~      7        218.      8.83      0  x=-12721~ "D:\\Cameca~        60
      # ... with 23 more variables: meas_bl.mt <dbl>, rejection.mt <dbl>,
      #   slit.mt <chr>, lens.mt <chr>, presput.mt <chr>, rast_com.mt <dbl>,
      #   frame.mt <chr>, blank_rast.mt <chr>, raster.mt <chr>, tune.mt <chr>,
      #   reg_mode.mt <chr>, chk_frm.mt <dbl>, sec_ion_cent.mt <chr>,
      #   frame_sec_ion_cent.mt <chr>, width_hor.mt <dbl>, width_ver.mt <dbl>,
      #   E0S_cent.mt <chr>, width_V.mt <dbl>, E0P_off.mt <dbl>,
      #   prim_cur_start.mt <chr>, prim_cur_after.mt <chr>, n.rw <dbl>,
      #   det_type.mt <chr>

# directory check

    Code
      ICdir_chk(point_example("2018-01-19-GLENDON"))
    Output
      $chk_is
                                                                                                2018-01-19-GLENDON_1_1 
      "/home/amandus/Documents/work/projects/code/point/inst/extdata/2018-01-19-GLENDON/2018-01-19-GLENDON_1_1.chk_is" 
                                                                                                2018-01-19-GLENDON_1_2 
      "/home/amandus/Documents/work/projects/code/point/inst/extdata/2018-01-19-GLENDON/2018-01-19-GLENDON_1_2.chk_is" 
                                                                                                2018-01-19-GLENDON_1_3 
      "/home/amandus/Documents/work/projects/code/point/inst/extdata/2018-01-19-GLENDON/2018-01-19-GLENDON_1_3.chk_is" 
      
      $is_txt
                                                                                                2018-01-19-GLENDON_1_1 
      "/home/amandus/Documents/work/projects/code/point/inst/extdata/2018-01-19-GLENDON/2018-01-19-GLENDON_1_1.is_txt" 
                                                                                                2018-01-19-GLENDON_1_2 
      "/home/amandus/Documents/work/projects/code/point/inst/extdata/2018-01-19-GLENDON/2018-01-19-GLENDON_1_2.is_txt" 
                                                                                                2018-01-19-GLENDON_1_3 
      "/home/amandus/Documents/work/projects/code/point/inst/extdata/2018-01-19-GLENDON/2018-01-19-GLENDON_1_3.is_txt" 
      
      $stat
                                                                                              2018-01-19-GLENDON_1_1 
      "/home/amandus/Documents/work/projects/code/point/inst/extdata/2018-01-19-GLENDON/2018-01-19-GLENDON_1_1.stat" 
                                                                                              2018-01-19-GLENDON_1_2 
      "/home/amandus/Documents/work/projects/code/point/inst/extdata/2018-01-19-GLENDON/2018-01-19-GLENDON_1_2.stat" 
                                                                                              2018-01-19-GLENDON_1_3 
      "/home/amandus/Documents/work/projects/code/point/inst/extdata/2018-01-19-GLENDON/2018-01-19-GLENDON_1_3.stat" 
      

