# Check change over time in read_IC

    Code
      read_IC(path_point)
    Output
      # A tibble: 81,900 x 7
         file.nm                 t.nm  N.rw species.nm sample.nm          n.rw bl.nm
         <chr>                  <dbl> <dbl> <chr>      <chr>             <dbl> <int>
       1 2018-01-19-GLENDON_1_1  0.54 12040 12C        Belemnite, Indium  3900     1
       2 2018-01-19-GLENDON_1_1  1.08 11950 12C        Belemnite, Indium  3900     1
       3 2018-01-19-GLENDON_1_1  1.62 12100 12C        Belemnite, Indium  3900     1
       4 2018-01-19-GLENDON_1_1  2.16 11946 12C        Belemnite, Indium  3900     1
       5 2018-01-19-GLENDON_1_1  2.7  12178 12C        Belemnite, Indium  3900     1
       6 2018-01-19-GLENDON_1_1  3.24 12114 12C        Belemnite, Indium  3900     1
       7 2018-01-19-GLENDON_1_1  3.78 12147 12C        Belemnite, Indium  3900     1
       8 2018-01-19-GLENDON_1_1  4.32 12092 12C        Belemnite, Indium  3900     1
       9 2018-01-19-GLENDON_1_1  4.86 12024 12C        Belemnite, Indium  3900     1
      10 2018-01-19-GLENDON_1_1  5.4  12045 12C        Belemnite, Indium  3900     1
      # ... with 81,890 more rows

# directory check

    Code
      ICdir_chk(point_example("2018-01-19-GLENDON"))
    Output
      $ion
               2018-01-19-GLENDON_1_1          2018-01-19-GLENDON_1_2 
      "2018-01-19-GLENDON_1_1.is_txt" "2018-01-19-GLENDON_1_2.is_txt" 
               2018-01-19-GLENDON_1_3 
      "2018-01-19-GLENDON_1_3.is_txt" 
      
      $optic
               2018-01-19-GLENDON_1_1          2018-01-19-GLENDON_1_2 
      "2018-01-19-GLENDON_1_1.chk_is" "2018-01-19-GLENDON_1_2.chk_is" 
               2018-01-19-GLENDON_1_3 
      "2018-01-19-GLENDON_1_3.chk_is" 
      
      $stat
             2018-01-19-GLENDON_1_1        2018-01-19-GLENDON_1_2 
      "2018-01-19-GLENDON_1_1.stat" "2018-01-19-GLENDON_1_2.stat" 
             2018-01-19-GLENDON_1_3 
      "2018-01-19-GLENDON_1_3.stat" 
      

