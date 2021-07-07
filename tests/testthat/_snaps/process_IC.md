# consistency of systematic corrections on internal dataset with custom settings

    Code
      cor_IC(tb_rw)
    Output
      # A tibble: 81,900 x 5
         file.nm                 t.nm species.nm  Xt.pr  N.pr
         <chr>                  <dbl> <chr>       <dbl> <dbl>
       1 2018-01-19-GLENDON_1_1  0.54 12C        34460. 12040
       2 2018-01-19-GLENDON_1_1  1.08 12C        34202. 11950
       3 2018-01-19-GLENDON_1_1  1.62 12C        34632. 12100
       4 2018-01-19-GLENDON_1_1  2.16 12C        34191. 11946
       5 2018-01-19-GLENDON_1_1  2.7  12C        34855. 12178
       6 2018-01-19-GLENDON_1_1  3.24 12C        34672. 12114
       7 2018-01-19-GLENDON_1_1  3.78 12C        34766. 12147
       8 2018-01-19-GLENDON_1_1  4.32 12C        34609. 12092
       9 2018-01-19-GLENDON_1_1  4.86 12C        34414. 12024
      10 2018-01-19-GLENDON_1_1  5.4  12C        34474. 12045
      # ... with 81,890 more rows

---

    Code
      cor_IC(tb_rw, .bl_t = 200)
    Output
      # A tibble: 81,900 x 5
         file.nm                 t.nm species.nm  Xt.pr  N.pr
         <chr>                  <dbl> <chr>       <dbl> <dbl>
       1 2018-01-19-GLENDON_1_1  0.54 12C        35412. 12040
       2 2018-01-19-GLENDON_1_1  1.08 12C        35147. 11950
       3 2018-01-19-GLENDON_1_1  1.62 12C        35588. 12100
       4 2018-01-19-GLENDON_1_1  2.16 12C        35135. 11946
       5 2018-01-19-GLENDON_1_1  2.7  12C        35818. 12178
       6 2018-01-19-GLENDON_1_1  3.24 12C        35629. 12114
       7 2018-01-19-GLENDON_1_1  3.78 12C        35726. 12147
       8 2018-01-19-GLENDON_1_1  4.32 12C        35565. 12092
       9 2018-01-19-GLENDON_1_1  4.86 12C        35365. 12024
      10 2018-01-19-GLENDON_1_1  5.4  12C        35426. 12045
      # ... with 81,890 more rows

---

    Code
      cor_IC(tb_rw, .bl_t = 200, .deadtime = 44, .det = "EM")
    Output
      # A tibble: 81,900 x 5
         file.nm                 t.nm species.nm  Xt.pr   N.pr
         <chr>                  <dbl> <chr>       <dbl>  <dbl>
       1 2018-01-19-GLENDON_1_1  0.54 12C        35467. 12059.
       2 2018-01-19-GLENDON_1_1  1.08 12C        35201. 11969.
       3 2018-01-19-GLENDON_1_1  1.62 12C        35644. 12119.
       4 2018-01-19-GLENDON_1_1  2.16 12C        35190. 11964.
       5 2018-01-19-GLENDON_1_1  2.7  12C        35874. 12197.
       6 2018-01-19-GLENDON_1_1  3.24 12C        35685. 12133.
       7 2018-01-19-GLENDON_1_1  3.78 12C        35783. 12166.
       8 2018-01-19-GLENDON_1_1  4.32 12C        35620. 12111.
       9 2018-01-19-GLENDON_1_1  4.86 12C        35420. 12043.
      10 2018-01-19-GLENDON_1_1  5.4  12C        35482. 12064.
      # ... with 81,890 more rows

---

    Code
      cor_IC(tb_rw, .bl_t = 200, .thr_PHD = 150, .M_PHD = 210, .SD_PHD = 60, .det = "EM")
    Output
      # A tibble: 81,900 x 5
         file.nm                 t.nm species.nm  Xt.pr   N.pr
         <chr>                  <dbl> <chr>       <dbl>  <dbl>
       1 2018-01-19-GLENDON_1_1  0.54 12C        42135. 14326.
       2 2018-01-19-GLENDON_1_1  1.08 12C        41820. 14219.
       3 2018-01-19-GLENDON_1_1  1.62 12C        42345. 14397.
       4 2018-01-19-GLENDON_1_1  2.16 12C        41806. 14214.
       5 2018-01-19-GLENDON_1_1  2.7  12C        42618. 14490.
       6 2018-01-19-GLENDON_1_1  3.24 12C        42394. 14414.
       7 2018-01-19-GLENDON_1_1  3.78 12C        42510. 14453.
       8 2018-01-19-GLENDON_1_1  4.32 12C        42317. 14388.
       9 2018-01-19-GLENDON_1_1  4.86 12C        42079. 14307.
      10 2018-01-19-GLENDON_1_1  5.4  12C        42153. 14332.
      # ... with 81,890 more rows

# consistency of systematic corrections on internal dataset with internal settings

    Code
      cor_IC(tb_rw)
    Output
      # A tibble: 81,900 x 7
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
      # ... with 81,890 more rows

---

    Code
      cor_IC(tb_rw, .bl_t = 200)
    Output
      # A tibble: 81,900 x 7
         file.nm                 t.nm species.nm sample.nm        bl.nm  Xt.pr  N.pr
         <chr>                  <dbl> <chr>      <chr>            <int>  <dbl> <dbl>
       1 2018-01-19-GLENDON_1_1  0.54 12C        Belemnite,Indium     1 35412. 12040
       2 2018-01-19-GLENDON_1_1  1.08 12C        Belemnite,Indium     1 35147. 11950
       3 2018-01-19-GLENDON_1_1  1.62 12C        Belemnite,Indium     1 35588. 12100
       4 2018-01-19-GLENDON_1_1  2.16 12C        Belemnite,Indium     1 35135. 11946
       5 2018-01-19-GLENDON_1_1  2.7  12C        Belemnite,Indium     1 35818. 12178
       6 2018-01-19-GLENDON_1_1  3.24 12C        Belemnite,Indium     1 35629. 12114
       7 2018-01-19-GLENDON_1_1  3.78 12C        Belemnite,Indium     1 35726. 12147
       8 2018-01-19-GLENDON_1_1  4.32 12C        Belemnite,Indium     1 35565. 12092
       9 2018-01-19-GLENDON_1_1  4.86 12C        Belemnite,Indium     1 35365. 12024
      10 2018-01-19-GLENDON_1_1  5.4  12C        Belemnite,Indium     1 35426. 12045
      # ... with 81,890 more rows

---

    Code
      cor_IC(tb_rw, .bl_t = 200, .deadtime = 44)
    Output
      # A tibble: 81,900 x 7
         file.nm                 t.nm species.nm sample.nm        bl.nm  Xt.pr   N.pr
         <chr>                  <dbl> <chr>      <chr>            <int>  <dbl>  <dbl>
       1 2018-01-19-GLENDON_1_1  0.54 12C        Belemnite,Indium     1 35467. 12059.
       2 2018-01-19-GLENDON_1_1  1.08 12C        Belemnite,Indium     1 35201. 11969.
       3 2018-01-19-GLENDON_1_1  1.62 12C        Belemnite,Indium     1 35644. 12119.
       4 2018-01-19-GLENDON_1_1  2.16 12C        Belemnite,Indium     1 35190. 11964.
       5 2018-01-19-GLENDON_1_1  2.7  12C        Belemnite,Indium     1 35874. 12197.
       6 2018-01-19-GLENDON_1_1  3.24 12C        Belemnite,Indium     1 35685. 12133.
       7 2018-01-19-GLENDON_1_1  3.78 12C        Belemnite,Indium     1 35783. 12166.
       8 2018-01-19-GLENDON_1_1  4.32 12C        Belemnite,Indium     1 35620. 12111.
       9 2018-01-19-GLENDON_1_1  4.86 12C        Belemnite,Indium     1 35420. 12043.
      10 2018-01-19-GLENDON_1_1  5.4  12C        Belemnite,Indium     1 35482. 12064.
      # ... with 81,890 more rows

---

    Code
      cor_IC(tb_rw, .bl_t = 200, .thr_PHD = 150)
    Output
      # A tibble: 81,900 x 7
         file.nm                 t.nm species.nm sample.nm        bl.nm  Xt.pr   N.pr
         <chr>                  <dbl> <chr>      <chr>            <int>  <dbl>  <dbl>
       1 2018-01-19-GLENDON_1_1  0.54 12C        Belemnite,Indium     1 35457. 12055.
       2 2018-01-19-GLENDON_1_1  1.08 12C        Belemnite,Indium     1 35192. 11965.
       3 2018-01-19-GLENDON_1_1  1.62 12C        Belemnite,Indium     1 35634. 12115.
       4 2018-01-19-GLENDON_1_1  2.16 12C        Belemnite,Indium     1 35180. 11961.
       5 2018-01-19-GLENDON_1_1  2.7  12C        Belemnite,Indium     1 35863. 12194.
       6 2018-01-19-GLENDON_1_1  3.24 12C        Belemnite,Indium     1 35675. 12129.
       7 2018-01-19-GLENDON_1_1  3.78 12C        Belemnite,Indium     1 35772. 12162.
       8 2018-01-19-GLENDON_1_1  4.32 12C        Belemnite,Indium     1 35610. 12107.
       9 2018-01-19-GLENDON_1_1  4.86 12C        Belemnite,Indium     1 35410. 12039.
      10 2018-01-19-GLENDON_1_1  5.4  12C        Belemnite,Indium     1 35472. 12060.
      # ... with 81,890 more rows

